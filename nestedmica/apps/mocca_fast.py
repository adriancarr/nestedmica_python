#!/usr/bin/env python3
"""
mocca_fast.py - Cython-optimized motif finder with Data-Oriented Threaded Parallelism.
Uses ThreadPoolExecutor + GIL-releasing Cython batch function. (No OpenMP dependency).
"""

import argparse
import numpy as np
from Bio import SeqIO
import time
import os

from nestedmica.trainer.fast import FastTrainer
from nestedmica.utils.console import AsciiPlotter
from nestedmica.utils.checkpoint import save_checkpoint, load_checkpoint

def main():
    parser = argparse.ArgumentParser(description='Fast Cython-optimized Mocca')
    parser.add_argument('-seqs', required=True, help='FASTA file')
    parser.add_argument('-numMotifs', type=int, default=1)
    parser.add_argument('-motifLength', type=int, default=10)
    parser.add_argument('-maxCycles', type=int, default=500)
    parser.add_argument('-ensembleSize', type=int, default=20)
    parser.add_argument('-threads', type=int, default=-1, help='Number of threads (-1=all)')
    parser.add_argument('-out', required=True, help='Output XMS file')
    parser.add_argument('-stopIqr', type=float, default=0.01)
    parser.add_argument('-checkpoint', help='File to save checkpoints to')
    parser.add_argument('-checkpointInterval', type=int, default=100)
    parser.add_argument('-restart', help='Checkpoint file to restart from')
    args = parser.parse_args()
    
    print(f"Loading sequences from {args.seqs}...")
    sequences = list(SeqIO.parse(args.seqs, "fasta"))
    print(f"Loaded {len(sequences)} sequences.")
    
    start_cycle = 0
    if args.restart:
        try:
            # Pass FastTrainer class to loader to avoid circular dependency
            trainer, start_cycle = load_checkpoint(args.restart, sequences, args.threads, FastTrainer)
            print(f"Resuming from Cycle {start_cycle}")
        except Exception as e:
            print(f"Failed to load checkpoint: {e}")
            return
    else:
        print("Initializing Cython Data-Oriented Threaded trainer...")
        trainer = FastTrainer(sequences, args.numMotifs, args.motifLength, 
                             args.ensembleSize, n_jobs=args.threads)
    
    print("Starting sampling...")
    start_time = time.time()
    
    plotter_l = AsciiPlotter(height=5, width=60)
    plotter_iqr = AsciiPlotter(height=5, width=60)
    
    for cycle in range(start_cycle, args.maxCycles):
        worst, hood = trainer.step()
        if cycle % 100 == 0:
            best_model, best_hood = trainer.get_best_model()

            # Calculate IQR
            likelihoods = np.array(trainer.model_likelihoods)
            q75, q25 = np.percentile(likelihoods, [75, 25])
            iqr = q75 - q25

            # Plotting
            plotter_l.add(best_hood)
            plotter_iqr.add(iqr)
            
            print(f"\nCycle {cycle}: L={best_hood:.2f} IQR={iqr:.4f}")
            print("Likelihood History:")
            print(plotter_l.plot())
            
            # Checkpoint
            if args.checkpoint and cycle % args.checkpointInterval == 0 and cycle > start_cycle:
                save_checkpoint(trainer, args.checkpoint, cycle)
            
            if iqr < args.stopIqr and cycle > start_cycle + 500:
                print(f"Converged at Cycle {cycle} (IQR {iqr:.4f} < {args.stopIqr})")
                break
    
    end_time = time.time()
    print(f"Finished in {end_time - start_time:.2f}s")
    
    # Save final checkpoint if requested
    if args.checkpoint:
         save_checkpoint(trainer, args.checkpoint, args.maxCycles)
    
    best_model, best_hood = trainer.get_best_model()
    print(f"Best Likelihood: {best_hood}")
    
    trainer.cleanup()
    
    with open(args.out, 'w') as f:
        f.write("<motifs>\n")
        f.write(f"<!-- Final Cycle: {args.maxCycles}, Best Likelihood: {best_hood} -->\n")
        for i, wm in enumerate(best_model['motifs']):
            f.write(f"  <motif id='{i}' weight='1.0'>\n")
            f.write(f"    <weightMatrix length='{wm.get_length()}'>\n")
            columns = wm.get_columns()
            for pos in range(wm.get_length()):
                probs = np.power(2.0, columns[:, pos])
                probs = probs / probs.sum()
                f.write(f"      {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}\n")
            f.write("    </weightMatrix>\n")
            f.write("  </motif>\n")
        f.write("</motifs>\n")
    
    print(f"Saved motifs to {args.out}")


if __name__ == "__main__":
    main()
