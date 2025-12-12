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
from nestedmica.utils.export import export_xms, export_meme, export_pfm, export_transfac
from nestedmica.utils.kmer_seeds import generate_seed_pwms

def main():
    parser = argparse.ArgumentParser(description='Fast Cython-optimized Mocca')
    parser.add_argument('-seqs', required=True, help='FASTA file')
    parser.add_argument('-numMotifs', type=int, default=1)
    parser.add_argument('-motifLength', type=int, default=10)
    parser.add_argument('-maxCycles', type=int, default=500)
    parser.add_argument('-ensembleSize', type=int, default=20)
    parser.add_argument('-threads', type=int, default=-1, help='Number of threads (-1=all)')
    parser.add_argument('-out', required=True, help='Output XMS file')
    parser.add_argument('-stopIqr', type=float, default=0.01, help='IQR convergence threshold (used with -convergenceMode iqr)')
    parser.add_argument('-stopH', type=float, default=0.1, help='H-hat convergence threshold (default mode)')
    parser.add_argument('-convergenceMode', choices=['skilling', 'iqr'], default='skilling', 
                        help='Convergence criterion: skilling (H-hat) or iqr (default: skilling)')
    parser.add_argument('-checkpoint', help='File to save checkpoints to')
    parser.add_argument('-checkpointInterval', type=int, default=100)
    parser.add_argument('-restart', help='Checkpoint file to restart from')
    parser.add_argument('-bgOrder', type=int, default=3, help='Background Markov order (0-5, default=3)')
    parser.add_argument('-adaptiveMCMC', action='store_true', default=False, help='[EXPERIMENTAL] Use adaptive proposals')
    parser.add_argument('-noAdaptiveMCMC', action='store_false', dest='adaptiveMCMC', help='Use fixed proposals (default)')
    parser.add_argument('-format', choices=['xms', 'meme', 'pfm', 'transfac'], default='xms',
                        help='Output format (default: xms)')
    parser.add_argument('-kmerSeeds', action='store_true', default=False,
                        help='Use k-mer enrichment for seed initialization (MEME-style)')
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
        # Generate k-mer seeds if requested
        seed_pwms = None
        if args.kmerSeeds:
            print("Generating k-mer enrichment seeds...")
            seed_pwms = generate_seed_pwms(sequences, args.numMotifs, args.motifLength)
            if seed_pwms:
                print(f"  Generated {len([s for s in seed_pwms if s is not None])} k-mer seeds")
            else:
                print("  No enriched k-mers found, using random initialization")
        
        print("Initializing Cython Data-Oriented Threaded trainer...")
        trainer = FastTrainer(sequences, args.numMotifs, args.motifLength, 
                             args.ensembleSize, n_jobs=args.threads, bg_order=args.bgOrder,
                             adaptive_mcmc=args.adaptiveMCMC, seed_pwms=seed_pwms)
    
    print("Starting sampling...")
    start_time = time.time()
    
    plotter_l = AsciiPlotter(height=5, width=60)
    plotter_iqr = AsciiPlotter(height=5, width=60)
    
    for cycle in range(start_cycle, args.maxCycles):
        worst, hood = trainer.step()
        if cycle % 100 == 0:
            best_model, best_hood = trainer.get_best_model()

            # Calculate convergence metrics
            likelihoods = np.array(trainer.model_likelihoods)
            q75, q25 = np.percentile(likelihoods, [75, 25])
            iqr = q75 - q25
            h_hat = trainer.get_skilling_h()
            remaining = trainer.get_remaining_info()

            # Plotting
            plotter_l.add(best_hood)
            plotter_iqr.add(iqr)
            
            print(f"\nCycle {cycle}: L={best_hood:.2f} LogZ={trainer.log_evidence:.2f} IQR={iqr:.4f} H={remaining:.2f}")
            print("Likelihood History:")
            print(plotter_l.plot())
            
            # Checkpoint
            if args.checkpoint and cycle % args.checkpointInterval == 0 and cycle > start_cycle:
                save_checkpoint(trainer, args.checkpoint, cycle)
            
            # Convergence check
            if cycle > start_cycle + 500:
                if args.convergenceMode == 'skilling' and remaining < args.stopH:
                    print(f"Converged at Cycle {cycle} (remaining info {remaining:.4f} < {args.stopH})")
                    break
                elif args.convergenceMode == 'iqr' and iqr < args.stopIqr:
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
    
    # Export motifs in selected format
    motifs = best_model['motifs']
    if args.format == 'xms':
        export_xms(motifs, args.out, log_evidence=trainer.log_evidence)
    elif args.format == 'meme':
        export_meme(motifs, args.out)
    elif args.format == 'pfm':
        export_pfm(motifs, args.out)
    elif args.format == 'transfac':
        export_transfac(motifs, args.out)
    
    print(f"Saved motifs to {args.out} ({args.format.upper()} format)")


if __name__ == "__main__":
    main()
