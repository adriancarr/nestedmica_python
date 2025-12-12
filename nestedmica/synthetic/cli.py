#!/usr/bin/env python3
"""
Command-line interface for synthetic data generation.

Usage:
    python -m nestedmica.synthetic generate -o output.fa --motif ACGTACGT
    python -m nestedmica.synthetic generate -o output.fa --motif "ACGT[5-10]TGCA" --learn-bg-from promoters.fa
"""

import argparse
import sys
from typing import Optional

from nestedmica.synthetic.background import (
    UniformBackground, LearnedBackground, 
    NormalLength, LogNormalLength, LearnedLength, FixedLength
)
from nestedmica.synthetic.motifs import (
    ConsensusMotif, GappedSyntheticMotif, PWMMotif
)
from nestedmica.synthetic.planting import UniformPlanting
from nestedmica.synthetic.datasets import (
    SimpleBenchmark, GappedMotifBenchmark, SyntheticDataset
)


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic DNA sequences with planted motifs',
        prog='python -m nestedmica.synthetic'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command')
    
    # Generate command
    gen = subparsers.add_parser('generate', help='Generate synthetic dataset')
    
    # Required
    gen.add_argument('-o', '--output', required=True,
                     help='Output FASTA file path')
    
    # Sequence parameters
    gen.add_argument('-n', '--num-sequences', type=int, default=100,
                     help='Number of sequences (default: 100)')
    gen.add_argument('-L', '--length', type=int, default=200,
                     help='Sequence length (default: 200)')
    gen.add_argument('--length-std', type=float, default=0,
                     help='Std dev for variable length (0=fixed, default: 0)')
    gen.add_argument('--length-min', type=int, default=50,
                     help='Minimum sequence length (default: 50)')
    gen.add_argument('--length-max', type=int, default=1000,
                     help='Maximum sequence length (default: 1000)')
    
    # Background parameters
    gen.add_argument('--gc', type=float, default=0.5,
                     help='GC content (0-1, default: 0.5)')
    gen.add_argument('--bg-order', type=int, default=0,
                     help='Background Markov order (default: 0 = IID)')
    gen.add_argument('--learn-bg-from', type=str, default=None,
                     help='FASTA file to learn background from')
    gen.add_argument('--bg-method', type=str, choices=['markov', 'shuffle'], default='markov',
                     help='Method for learned background: markov (generative) or shuffle (exact dinucl).')
    
    # Motif parameters
    gen.add_argument('--motif', type=str, action='append',
                     help='Motif pattern (can specify multiple). '
                          'Use IUPAC codes or [min-max] for gaps.')
    gen.add_argument('--motif-file', type=str,
                     help='MEME or JASPAR file with motifs')
    
    # Planting parameters
    gen.add_argument('--density', type=float, default=0.8,
                     help='Fraction of sequences with motif (default: 0.8)')
    gen.add_argument('--strand-bias', type=float, default=0.5,
                     help='Forward strand probability (default: 0.5)')
    gen.add_argument('--occurrences', type=int, default=1,
                     help='Motif occurrences per sequence (default: 1)')
    gen.add_argument('--planting', type=str, choices=['uniform', 'center'], default='uniform',
                     help='Planting strategy: uniform or center (Gaussian bias).')
    
    # Other
    gen.add_argument('--seed', type=int, default=42,
                     help='Random seed (default: 42)')
    gen.add_argument('--truth', type=str, default=None,
                     help='Output JSON file for ground truth')
                     
    # Benchmark command
    bench = subparsers.add_parser('benchmark', help='Generate full train/test benchmark suite')
    bench.add_argument('-o', '--output-dir', required=True, help='Output directory')
    bench.add_argument('--train', type=int, default=1000, help='Train size')
    bench.add_argument('--test', type=int, default=200, help='Test size')
    bench.add_argument('--motif', type=str, required=True, help='Motif pattern')
    bench.add_argument('--bg-fasta', type=str, help='Learned background FASTA')
    
    args = parser.parse_args()
    
    if args.command == 'benchmark':
        from nestedmica.synthetic.datasets import BenchmarkGenerator
        
        gen = BenchmarkGenerator(args.output_dir)
        config = {'motif': args.motif}
        
        if args.bg_fasta:
            # For benchmark, we assume simple config for now or expand this
            # Passing raw config dict for flexibility
            config['background_fasta'] = args.bg_fasta
            
        print(f"Generating benchmark in {args.output_dir}...")
        gen.create_standard_suite(config, train_size=args.train, test_size=args.test)
        print("Done.")
        return

    if args.command != 'generate':
        parser.print_help()
        return
    
    # Build length distribution
    if args.length_std > 0:
        length_dist = NormalLength(
            mean=args.length, 
            std=args.length_std,
            min_len=args.length_min,
            max_len=args.length_max
        )
    else:
        length_dist = FixedLength(args.length)
    
    # Build background
    if args.learn_bg_from:
        if args.bg_method == 'shuffle':
            from nestedmica.synthetic.background import ShuffledBackground
            print(f"Creating shuffled background from {args.learn_bg_from}...")
            background = ShuffledBackground.from_fasta(args.learn_bg_from)
        else:
            print(f"Learning background from {args.learn_bg_from}...")
            background = LearnedBackground.from_fasta(args.learn_bg_from, order=args.bg_order)
            print(f"  Learned {args.bg_order}-order Markov model, GC={background.gc_content:.1%}")
    else:
        background = UniformBackground(gc_content=args.gc)
    
    # Create dataset
    dataset = SyntheticDataset(
        n_sequences=args.num_sequences,
        seq_length=length_dist,
        background=background,
        seed=args.seed
    )
    
    # Add motifs
    motifs = args.motif or []
    
    if args.motif_file:
        if args.motif_file.endswith('.meme'):
            pwm_motifs = PWMMotif.from_meme(args.motif_file)
            for m in pwm_motifs:
                dataset.add_motif(m)
                print(f"  Loaded motif: {m.name} ({m.length}bp)")
        else:
            pwm = PWMMotif.from_jaspar_file(args.motif_file)
            dataset.add_motif(pwm)
            print(f"  Loaded motif: {pwm.name} ({pwm.length}bp)")
    
    for pattern in motifs:
        dataset.add_motif(pattern)
        print(f"  Added motif: {pattern}")
    
    if not dataset.motifs:
        print("Warning: No motifs specified. Generating background-only sequences.")
    
    # Set planting strategy
    if args.planting == 'center':
        from nestedmica.synthetic.planting import PositionalPlanting
        strategy = PositionalPlanting(
            density=args.density,
            strand_bias=args.strand_bias,
            center_fraction=0.5,
            std_fraction=0.1
        )
    else:
        strategy = UniformPlanting(
            density=args.density,
            strand_bias=args.strand_bias,
            occurrences_per_seq=args.occurrences
        )
        
    dataset.set_planting(strategy)
    
    # Generate
    print(f"Generating {args.num_sequences} sequences...")
    dataset.generate(args.output)
    print(f"Wrote {args.num_sequences} sequences to {args.output}")
    
    # Save truth if requested
    if args.truth:
        dataset.save_truth(args.truth)
        print(f"Wrote ground truth to {args.truth}")
    else:
        # Auto-generate truth filename
        truth_path = args.output.replace('.fa', '_truth.json').replace('.fasta', '_truth.json')
        dataset.save_truth(truth_path)
        print(f"Wrote ground truth to {truth_path}")


if __name__ == '__main__':
    main()
