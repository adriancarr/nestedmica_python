#!/usr/bin/env python3
"""
NestedMICA Auto-Discovery Tool
Automatically determines the optimal number of motifs using Bayesian Model Selection.
"""

import argparse
import time
import numpy as np
from Bio import SeqIO
from nestedmica.trainer.fast import FastTrainer
from nestedmica.utils.console import AsciiPlotter

def run_model(sequences, num_motifs, args):
    """Run a single model configuration and return (log_evidence, best_model)."""
    print(f"\n{'='*60}")
    print(f"  Testing Hypothesis: N = {num_motifs} Motifs")
    print(f"{'='*60}")
    
    # Initialize Trainer
    # We use a default starting length, let Indel/Zap handle the rest.
    start_length = args.startLength 
    trainer = FastTrainer(sequences, num_motifs, start_length, 
                          args.ensembleSize, n_jobs=args.threads)
    
    plotter_iqr = AsciiPlotter(height=3, width=40)
    
    start_time = time.time()
    converged = False
    
    # Run loop
    for cycle in range(args.maxCycles):
        worst, hood = trainer.step()
        
        # Monitor convergence
        if cycle % 100 == 0 and cycle > 0:
            likelihoods = np.array(trainer.model_likelihoods)
            q75, q25 = np.percentile(likelihoods, [75, 25])
            iqr = q75 - q25
            
            # Print brief status
            elapsed = time.time() - start_time
            print(f"Cycle {cycle}: LogZ={trainer.log_evidence:.2f} IQR={iqr:.4f} ({elapsed:.1f}s)")
            plotter_iqr.add(iqr)
            print(plotter_iqr.plot())
            
            # Convergence Check
            if iqr < args.stopIqr:
                print(f"--> Converged at Cycle {cycle} (IQR < {args.stopIqr})")
                converged = True
                break
    
    if not converged:
        print(f"--> Stopped at max cycles ({args.maxCycles})")
        
    best_model, best_hood = trainer.get_best_model()
    return trainer.log_evidence, best_model

def save_result(filename, model, log_evidence, num_motifs):
    """Save the winning model."""
    with open(filename, 'w') as f:
        f.write("<motifs>\n")
        f.write(f"<!-- Auto-Discovered Model -->\n")
        f.write(f"<!-- Optimal Motifs: {num_motifs} -->\n")
        f.write(f"<!-- GlobalLogEvidence: {log_evidence:.4f} -->\n")
        
        for i, wm in enumerate(model['motifs']):
             f.write(f"  <motif id='{i}' weight='1.0'>\n")
             f.write(f"    <weightMatrix length='{wm.get_length()}'>\n")
             cols = wm.get_columns()
             length = cols.shape[1]
             for pos in range(length):
                 probs = np.power(2.0, cols[:, pos])
                 probs = probs / probs.sum()
                 f.write(f"      {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}\n")
             f.write("    </weightMatrix>\n")
             f.write("  </motif>\n")
        f.write("</motifs>\n")
    print(f"Saved optimal solution to {filename}")

def main():
    parser = argparse.ArgumentParser(description='NestedMICA Auto-Discovery')
    parser.add_argument('-seqs', required=True, help='Input FASTA sequences')
    parser.add_argument('-out', required=True, help='Output XMS file')
    parser.add_argument('-maxMotifs', type=int, default=10, help='Maximum N to test')
    parser.add_argument('-minMotifs', type=int, default=1, help='Minimum N to start search from')
    parser.add_argument('-bayesThreshold', type=float, default=5.0, help='Log Bayes Factor threshold to accept new motif')
    
    # Tuning params
    parser.add_argument('-ensembleSize', type=int, default=None, help='Ensemble size (default: Auto-calculated)')
    parser.add_argument('-maxCycles', type=int, default=5000, help='Max cycles per run')
    parser.add_argument('-stopIqr', type=float, default=0.5, help='Convergence threshold')
    parser.add_argument('-startLength', type=int, default=10, help='Initial motif length')
    parser.add_argument('-threads', type=int, default=-1, help='Number of threads')
    
    args = parser.parse_args()
    
    print(f"Loading sequences from {args.seqs}...")
    sequences = list(SeqIO.parse(args.seqs, "fasta"))
    print(f"Loaded {len(sequences)} sequences.")
    
    # Smart Ensemble Sizing
    if args.ensembleSize is None:
        total_bases = sum(len(s.seq) for s in sequences)
        # Formula: 50 * log10(bases). E.g. 1MB -> 50 * 6 = 300. 10KB -> 50 * 4 = 200.
        # Minimum of 50 to ensure coverage.
        calc_size = int(50 * np.log10(max(total_bases, 100)))
        args.ensembleSize = max(50, calc_size)
        print(f"Auto-configured Ensemble Size: {args.ensembleSize} (based on {total_bases} bases)")
    else:
        print(f"Using manual Ensemble Size: {args.ensembleSize}")
    
    history = [] # List of (N, LogZ)    best_model_overall = None
    best_N = 0
    max_log_z = -np.float64('inf')
    
    print("\nStarting Bayesian Model Selection Scan...")
    print(f"Range: N={args.minMotifs} to {args.maxMotifs}")
    print("Criterion: Log Bayes Factor >", args.bayesThreshold)
    
    for n in range(args.minMotifs, args.maxMotifs + 1):
        log_z, model = run_model(sequences, n, args)
        history.append((n, log_z))
        
        # Compare with previous
        if n == args.minMotifs:
            print(f"Base Evidence (N={n}): {log_z:.2f}")
            max_log_z = log_z
            best_model_overall = model
            best_N = n
        else:
            prev_n, prev_z = history[-2]
            bayes_factor = log_z - prev_z
            
            print(f"\nComparator: N={prev_n} -> N={n}")
            print(f"  Delta LogZ (Bayes Factor): {bayes_factor:.2f}")
            
            if bayes_factor > args.bayesThreshold:
                print(f"  Decision: ACCEPT (Evidence is strong)")
                max_log_z = log_z
                best_model_overall = model
                best_N = n
            else:
                print(f"  Decision: REJECT (Insufficient evidence)")
                print(f"Stopping scan. Optimal complexity reached at N={prev_n}.")
                break
    
    print("\n" + "="*60)
    print(f"FINAL RESULT: {best_N} Motifs found.")
    print("="*60)
    
    if best_model_overall:
        save_result(args.out, best_model_overall, max_log_z, best_N)

if __name__ == '__main__':
    main()
