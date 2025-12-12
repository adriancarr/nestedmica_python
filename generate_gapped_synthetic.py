#!/usr/bin/env python3
"""
Generate synthetic sequences with planted gapped motifs.

Usage:
    python generate_gapped_synthetic.py output.fa -n 100 -L 200 \
        -motif "ACGTACGT[5-8]TGCATGCA" -density 0.8

The motif format uses [min-max] to specify gap length ranges.
"""

import argparse
import numpy as np
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_gapped_motif(motif_str: str):
    """
    Parse a gapped motif string like "ACGT[5-8]TGCA".
    
    Returns:
        List of (block_consensus, gap_range) tuples.
        gap_range is (min, max) or None for no gap after this block.
    """
    # Split on gap pattern [min-max]
    pattern = r'\[(\d+)-(\d+)\]'
    
    blocks = []
    remaining = motif_str
    
    while remaining:
        match = re.search(pattern, remaining)
        if match:
            block_seq = remaining[:match.start()]
            gap_min = int(match.group(1))
            gap_max = int(match.group(2))
            blocks.append((block_seq, (gap_min, gap_max)))
            remaining = remaining[match.end():]
        else:
            # Last block (no gap after)
            blocks.append((remaining, None))
            break
    
    return blocks


def consensus_to_pwm(consensus: str) -> np.ndarray:
    """Convert consensus string to PWM (log2 space)."""
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    length = len(consensus)
    pwm = np.full((4, length), -10.0)  # Low probability for non-consensus
    
    for i, base in enumerate(consensus):
        if base in base_map:
            pwm[base_map[base], i] = 0.0  # High probability for consensus
    
    return pwm


def sample_from_pwm(pwm: np.ndarray, rng: np.random.Generator) -> str:
    """Sample a sequence from a PWM."""
    bases = ['A', 'C', 'G', 'T']
    length = pwm.shape[1]
    result = []
    
    for i in range(length):
        # Convert log2 to probabilities
        probs = np.power(2, pwm[:, i])
        probs = probs / probs.sum()
        result.append(bases[rng.choice(4, p=probs)])
    
    return ''.join(result)


def generate_background(length: int, rng: np.random.Generator) -> str:
    """Generate random background sequence."""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(rng.choice(bases, size=length))


def plant_gapped_motif(seq_length: int, blocks: list, rng: np.random.Generator) -> tuple:
    """
    Plant a gapped motif into a random sequence.
    
    Returns:
        Tuple of (sequence, motif_start, gap_lengths)
    """
    # Calculate total min/max span
    total_block_len = sum(len(b[0]) for b in blocks)
    min_gaps = sum(g[0] for b, g in blocks if g is not None)
    max_gaps = sum(g[1] for b, g in blocks if g is not None)
    
    min_span = total_block_len + min_gaps
    max_span = total_block_len + max_gaps
    
    if max_span > seq_length - 10:
        raise ValueError(f"Motif span {max_span} too large for sequence length {seq_length}")
    
    # Choose motif start position
    start = rng.integers(5, seq_length - max_span - 5)
    
    # Generate the motif instance
    motif_parts = []
    gap_lengths = []
    
    for block_seq, gap_range in blocks:
        # Sample block with some variation
        pwm = consensus_to_pwm(block_seq)
        motif_parts.append(sample_from_pwm(pwm, rng))
        
        if gap_range is not None:
            gap_len = rng.integers(gap_range[0], gap_range[1] + 1)
            gap_lengths.append(gap_len)
            # Gap is random background
            motif_parts.append(generate_background(gap_len, rng))
    
    motif_instance = ''.join(motif_parts)
    
    # Build full sequence
    prefix = generate_background(start, rng)
    suffix = generate_background(seq_length - start - len(motif_instance), rng)
    
    sequence = prefix + motif_instance + suffix
    
    return sequence, start, gap_lengths


def main():
    parser = argparse.ArgumentParser(description='Generate synthetic gapped motif data')
    parser.add_argument('output', help='Output FASTA file')
    parser.add_argument('-n', type=int, default=100, help='Number of sequences')
    parser.add_argument('-L', type=int, default=200, help='Sequence length')
    parser.add_argument('-motif', default='ACGTACGT[5-10]TGCATGCA', 
                        help='Gapped motif pattern, e.g. ACGT[5-8]TGCA')
    parser.add_argument('-density', type=float, default=0.8,
                        help='Fraction of sequences with motif (default: 0.8)')
    parser.add_argument('-seed', type=int, default=42, help='Random seed')
    args = parser.parse_args()
    
    rng = np.random.default_rng(args.seed)
    
    # Parse motif
    blocks = parse_gapped_motif(args.motif)
    print(f"Gapped motif structure:")
    for i, (block_seq, gap_range) in enumerate(blocks):
        print(f"  Block {i}: {block_seq} ({len(block_seq)}bp)")
        if gap_range:
            print(f"  Gap {i}: {gap_range[0]}-{gap_range[1]}bp")
    
    # Generate sequences
    records = []
    truth = []
    
    for i in range(args.n):
        if rng.random() < args.density:
            seq, start, gaps = plant_gapped_motif(args.L, blocks, rng)
            truth.append({
                'seq_id': f'seq_{i}',
                'has_motif': True,
                'start': start,
                'gap_lengths': gaps
            })
        else:
            seq = generate_background(args.L, rng)
            truth.append({
                'seq_id': f'seq_{i}',
                'has_motif': False
            })
        
        record = SeqRecord(Seq(seq), id=f'seq_{i}', description='')
        records.append(record)
    
    # Write FASTA
    SeqIO.write(records, args.output, 'fasta')
    print(f"Wrote {len(records)} sequences to {args.output}")
    
    # Write truth file
    truth_file = args.output.replace('.fa', '_truth.txt').replace('.fasta', '_truth.txt')
    with open(truth_file, 'w') as f:
        f.write(f"# Motif: {args.motif}\n")
        f.write(f"# Density: {args.density}\n")
        planted = sum(1 for t in truth if t['has_motif'])
        f.write(f"# Planted: {planted}/{args.n}\n")
        for t in truth:
            if t['has_motif']:
                f.write(f"{t['seq_id']}\t{t['start']}\t{t['gap_lengths']}\n")
    print(f"Wrote truth to {truth_file}")


if __name__ == '__main__':
    main()
