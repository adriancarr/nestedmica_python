#!/usr/bin/env python3
"""
Generate synthetic FASTA sequences with a known planted motif.

This script creates sequences that contain a specific motif pattern,
allowing us to verify that both Java and Python implementations
correctly recover the planted motif.
"""

import random
import numpy as np

# Configuration
NUM_SEQUENCES = 100
SEQ_LENGTH = 100
MOTIF_LENGTH = 10

# The planted motif PWM (probability distribution at each position)
# This is a strong, unambiguous motif for testing
PLANTED_MOTIF = np.array([
    [0.9, 0.03, 0.03, 0.04],  # A-rich
    [0.03, 0.9, 0.04, 0.03],  # C-rich
    [0.03, 0.04, 0.9, 0.03],  # G-rich
    [0.04, 0.03, 0.03, 0.9],  # T-rich
    [0.9, 0.03, 0.03, 0.04],  # A-rich
    [0.03, 0.9, 0.04, 0.03],  # C-rich
    [0.03, 0.04, 0.9, 0.03],  # G-rich
    [0.04, 0.03, 0.03, 0.9],  # T-rich
    [0.5, 0.5, 0.0, 0.0],     # A or C
    [0.0, 0.0, 0.5, 0.5],     # G or T
])

BASES = ['A', 'C', 'G', 'T']

def sample_from_pwm(pwm):
    """Sample a sequence from a PWM."""
    seq = ""
    for row in pwm:
        base = np.random.choice(BASES, p=row)
        seq += base
    return seq

def generate_background(length):
    """Generate uniform random background sequence."""
    return ''.join(random.choice(BASES) for _ in range(length))

def generate_sequence_with_motif(seq_length, motif_pwm, num_occurrences=1):
    """Generate a sequence with planted motif occurrences."""
    motif_len = len(motif_pwm)
    
    # Start with background
    seq = list(generate_background(seq_length))
    
    # Plant motifs at random positions
    positions = []
    for _ in range(num_occurrences):
        # Find valid position (not overlapping with existing)
        max_attempts = 100
        for _ in range(max_attempts):
            pos = random.randint(0, seq_length - motif_len)
            # Check overlap
            overlap = False
            for existing_pos in positions:
                if abs(pos - existing_pos) < motif_len:
                    overlap = True
                    break
            if not overlap:
                positions.append(pos)
                break
    
    # Insert motif instances
    for pos in positions:
        motif_instance = sample_from_pwm(motif_pwm)
        for i, base in enumerate(motif_instance):
            seq[pos + i] = base
    
    return ''.join(seq), positions

def main():
    output_file = "synthetic_test.fa"
    ground_truth_file = "synthetic_ground_truth.txt"
    
    with open(output_file, 'w') as f, open(ground_truth_file, 'w') as gt:
        gt.write("# Ground truth motif positions\n")
        gt.write("# Planted motif consensus: ACGTACGT[AC][GT]\n")
        gt.write(f"# Motif length: {MOTIF_LENGTH}\n")
        gt.write("# PWM (A C G T):\n")
        for i, row in enumerate(PLANTED_MOTIF):
            gt.write(f"#   pos {i}: {row[0]:.2f} {row[1]:.2f} {row[2]:.2f} {row[3]:.2f}\n")
        gt.write("\n")
        
        for i in range(NUM_SEQUENCES):
            seq, positions = generate_sequence_with_motif(SEQ_LENGTH, PLANTED_MOTIF, num_occurrences=1)
            f.write(f">seq{i}\n")
            f.write(f"{seq}\n")
            gt.write(f"seq{i}: {positions}\n")
    
    print(f"Generated {NUM_SEQUENCES} sequences to {output_file}")
    print(f"Ground truth written to {ground_truth_file}")
    print(f"Expected motif consensus: ACGTACGT[AC][GT]")

if __name__ == "__main__":
    main()
