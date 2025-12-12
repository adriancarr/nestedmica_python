#!/usr/bin/env python3
"""
Generate a challenging synthetic FASTA dataset with:
- 10 different motifs with varying signal strengths
- 0-3 motif occurrences per sequence
- 500 sequences
- Variable background
"""

import random
import numpy as np

# Configuration
NUM_SEQUENCES = 100
SEQ_LENGTH = 200
MOTIF_LENGTH = 10
NUM_MOTIFS = 10

BASES = ['A', 'C', 'G', 'T']

def generate_random_motif(strength):
    """
    Generate a random motif PWM with given signal strength.
    strength: 0.0 (uniform/weak) to 1.0 (deterministic/strong)
    """
    pwm = np.zeros((MOTIF_LENGTH, 4))
    
    for pos in range(MOTIF_LENGTH):
        # Pick a "preferred" base for this position
        preferred = np.random.randint(4)
        
        # Set probabilities based on strength
        # strength=1.0 -> 100% preferred, 0% others
        # strength=0.5 -> 62.5% preferred, 12.5% each other
        # strength=0.25 -> ~44% preferred, ~19% each other
        for base in range(4):
            if base == preferred:
                pwm[pos, base] = 0.25 + 0.75 * strength
            else:
                pwm[pos, base] = (1.0 - (0.25 + 0.75 * strength)) / 3.0
    
    # Normalize rows
    pwm = pwm / pwm.sum(axis=1, keepdims=True)
    return pwm

def get_consensus(pwm):
    """Get consensus sequence from PWM."""
    return ''.join(BASES[np.argmax(row)] for row in pwm)

def sample_from_pwm(pwm):
    """Sample a sequence from a PWM."""
    seq = ""
    for row in pwm:
        base = np.random.choice(BASES, p=row)
        seq += base
    return seq

def generate_background(length, gc_bias=0.5):
    """Generate random background with optional GC bias."""
    probs = [(1-gc_bias)/2, gc_bias/2, gc_bias/2, (1-gc_bias)/2]  # A, C, G, T
    return ''.join(np.random.choice(BASES, p=probs) for _ in range(length))

def generate_sequence_with_motifs(seq_length, motif_pwms, num_occurrences):
    """Generate a sequence with multiple planted motif occurrences."""
    motif_len = MOTIF_LENGTH
    
    # Start with background
    seq = list(generate_background(seq_length, gc_bias=0.5))
    
    # Track planted motifs
    planted = []
    positions = []
    
    for _ in range(num_occurrences):
        # Pick a random motif
        motif_idx = np.random.randint(len(motif_pwms))
        pwm = motif_pwms[motif_idx]
        
        # Find valid position (not overlapping with existing)
        max_attempts = 50
        for _ in range(max_attempts):
            pos = random.randint(0, seq_length - motif_len)
            # Check overlap
            overlap = False
            for existing_pos, _, _ in planted:
                if abs(pos - existing_pos) < motif_len + 2:  # Buffer of 2bp
                    overlap = True
                    break
            if not overlap:
                break
        else:
            continue  # Couldn't find non-overlapping position
        
        # Sample and insert motif instance
        motif_instance = sample_from_pwm(pwm)
        for i, base in enumerate(motif_instance):
            seq[pos + i] = base
        
        planted.append((pos, motif_idx, motif_instance))
        positions.append(pos)
    
    return ''.join(seq), planted

def main():
    np.random.seed(42)  # Reproducible
    random.seed(42)
    
    output_file = "synthetic_hard_test.fa"
    ground_truth_file = "synthetic_hard_ground_truth.txt"
    
    # Generate 10 motifs with varying strengths
    # Strengths from very strong (0.95) to weak (0.40)
    strengths = np.linspace(0.95, 0.40, NUM_MOTIFS)
    motif_pwms = []
    
    print("Generating motifs with varying signal strengths:")
    for i, strength in enumerate(strengths):
        pwm = generate_random_motif(strength)
        motif_pwms.append(pwm)
        consensus = get_consensus(pwm)
        print(f"  Motif {i}: {consensus} (strength={strength:.2f})")
    
    # Statistics
    total_occurrences = {i: 0 for i in range(NUM_MOTIFS)}
    seqs_with_n_motifs = {0: 0, 1: 0, 2: 0, 3: 0}
    
    with open(output_file, 'w') as f, open(ground_truth_file, 'w') as gt:
        gt.write("# Challenging Synthetic Test - Ground Truth\n")
        gt.write(f"# {NUM_SEQUENCES} sequences, {SEQ_LENGTH}bp each\n")
        gt.write(f"# {NUM_MOTIFS} motifs, length {MOTIF_LENGTH}\n")
        gt.write("# 0-3 motif occurrences per sequence\n\n")
        
        gt.write("# Motif PWMs (consensus, strength):\n")
        for i, (pwm, strength) in enumerate(zip(motif_pwms, strengths)):
            gt.write(f"# Motif {i}: {get_consensus(pwm)} (strength={strength:.2f})\n")
        gt.write("\n")
        
        # Save full PWMs
        gt.write("# Full PWMs (A C G T):\n")
        for i, pwm in enumerate(motif_pwms):
            gt.write(f"MOTIF_{i}_PWM:\n")
            for pos, row in enumerate(pwm):
                gt.write(f"  {pos}: {row[0]:.3f} {row[1]:.3f} {row[2]:.3f} {row[3]:.3f}\n")
        gt.write("\n# Sequence annotations:\n")
        
        for seq_num in range(NUM_SEQUENCES):
            # Random number of motifs (0-3)
            num_motifs = np.random.choice([0, 1, 2, 3], p=[0.10, 0.40, 0.35, 0.15])
            seqs_with_n_motifs[num_motifs] += 1
            
            seq, planted = generate_sequence_with_motifs(SEQ_LENGTH, motif_pwms, num_motifs)
            
            f.write(f">seq{seq_num}\n")
            f.write(f"{seq}\n")
            
            # Record ground truth
            if planted:
                for pos, motif_idx, instance in planted:
                    gt.write(f"seq{seq_num}: pos={pos}, motif={motif_idx}, seq={instance}\n")
                    total_occurrences[motif_idx] += 1
            else:
                gt.write(f"seq{seq_num}: (no motifs)\n")
    
    print(f"\nGenerated {NUM_SEQUENCES} sequences to {output_file}")
    print(f"Ground truth written to {ground_truth_file}")
    
    print("\n### Statistics ###")
    print(f"Sequences with 0 motifs: {seqs_with_n_motifs[0]}")
    print(f"Sequences with 1 motif:  {seqs_with_n_motifs[1]}")
    print(f"Sequences with 2 motifs: {seqs_with_n_motifs[2]}")
    print(f"Sequences with 3 motifs: {seqs_with_n_motifs[3]}")
    print("\nMotif occurrence counts:")
    for i in range(NUM_MOTIFS):
        print(f"  Motif {i} (strength={strengths[i]:.2f}): {total_occurrences[i]} occurrences")

if __name__ == "__main__":
    main()
