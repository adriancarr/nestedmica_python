#!/usr/bin/env python3
"""
Generate medium-difficulty synthetic test:
- 50 sequences
- 2 motifs: one very strong (0.95), one medium (0.70)
- 0-2 motif occurrences per sequence
"""

import random
import numpy as np

# Configuration
NUM_SEQUENCES = 50
SEQ_LENGTH = 100
MOTIF_LENGTH = 10
NUM_MOTIFS = 2

BASES = ['A', 'C', 'G', 'T']

# Define specific motifs with different strengths
MOTIF_CONFIGS = [
    {"consensus": "ACGTACGTAC", "strength": 0.95, "name": "strong"},  # Very strong
    {"consensus": "GCTAGCTAGT", "strength": 0.70, "name": "medium"},  # Medium
]

def consensus_to_pwm(consensus, strength):
    """Convert consensus to PWM with given strength."""
    pwm = np.zeros((len(consensus), 4))
    for i, base in enumerate(consensus):
        base_idx = BASES.index(base)
        for j in range(4):
            if j == base_idx:
                pwm[i, j] = 0.25 + 0.75 * strength
            else:
                pwm[i, j] = (1.0 - (0.25 + 0.75 * strength)) / 3.0
    return pwm

def sample_from_pwm(pwm):
    """Sample a sequence from a PWM."""
    seq = ""
    for row in pwm:
        base = np.random.choice(BASES, p=row)
        seq += base
    return seq

def generate_background(length):
    """Generate uniform random background."""
    return ''.join(random.choice(BASES) for _ in range(length))

def generate_sequence_with_motifs(seq_length, motif_pwms, num_occurrences):
    """Generate a sequence with planted motif occurrences."""
    motif_len = MOTIF_LENGTH
    seq = list(generate_background(seq_length))
    planted = []
    
    for _ in range(num_occurrences):
        motif_idx = np.random.randint(len(motif_pwms))
        pwm = motif_pwms[motif_idx]
        
        max_attempts = 50
        for _ in range(max_attempts):
            pos = random.randint(0, seq_length - motif_len)
            overlap = False
            for existing_pos, _, _ in planted:
                if abs(pos - existing_pos) < motif_len + 2:
                    overlap = True
                    break
            if not overlap:
                break
        else:
            continue
        
        motif_instance = sample_from_pwm(pwm)
        for i, base in enumerate(motif_instance):
            seq[pos + i] = base
        planted.append((pos, motif_idx, motif_instance))
    
    return ''.join(seq), planted

def main():
    np.random.seed(123)
    random.seed(123)
    
    output_file = "synthetic_medium_test.fa"
    ground_truth_file = "synthetic_medium_ground_truth.txt"
    
    # Generate PWMs
    motif_pwms = []
    print("Generating motifs:")
    for i, config in enumerate(MOTIF_CONFIGS):
        pwm = consensus_to_pwm(config["consensus"], config["strength"])
        motif_pwms.append(pwm)
        print(f"  Motif {i} ({config['name']}): {config['consensus']} (strength={config['strength']:.2f})")
    
    # Statistics
    total_occurrences = {i: 0 for i in range(NUM_MOTIFS)}
    seqs_with_n_motifs = {0: 0, 1: 0, 2: 0}
    
    with open(output_file, 'w') as f, open(ground_truth_file, 'w') as gt:
        gt.write("# Medium Difficulty Synthetic Test - Ground Truth\n")
        gt.write(f"# {NUM_SEQUENCES} sequences, {SEQ_LENGTH}bp each\n")
        gt.write(f"# {NUM_MOTIFS} motifs, length {MOTIF_LENGTH}\n\n")
        
        for i, config in enumerate(MOTIF_CONFIGS):
            gt.write(f"# Motif {i} ({config['name']}): {config['consensus']} (strength={config['strength']:.2f})\n")
        gt.write("\n")
        
        for seq_num in range(NUM_SEQUENCES):
            # 0-2 motifs per sequence
            num_motifs = np.random.choice([0, 1, 2], p=[0.15, 0.50, 0.35])
            seqs_with_n_motifs[num_motifs] += 1
            
            seq, planted = generate_sequence_with_motifs(SEQ_LENGTH, motif_pwms, num_motifs)
            
            f.write(f">seq{seq_num}\n")
            f.write(f"{seq}\n")
            
            if planted:
                for pos, motif_idx, instance in planted:
                    gt.write(f"seq{seq_num}: pos={pos}, motif={motif_idx}, seq={instance}\n")
                    total_occurrences[motif_idx] += 1
            else:
                gt.write(f"seq{seq_num}: (no motifs)\n")
    
    print(f"\nGenerated {NUM_SEQUENCES} sequences to {output_file}")
    print(f"\n### Statistics ###")
    print(f"Sequences with 0 motifs: {seqs_with_n_motifs[0]}")
    print(f"Sequences with 1 motif:  {seqs_with_n_motifs[1]}")
    print(f"Sequences with 2 motifs: {seqs_with_n_motifs[2]}")
    print(f"\nMotif occurrence counts:")
    for i, config in enumerate(MOTIF_CONFIGS):
        print(f"  Motif {i} ({config['name']}): {total_occurrences[i]} occurrences")

if __name__ == "__main__":
    main()
