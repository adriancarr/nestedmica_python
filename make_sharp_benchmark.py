#!/usr/bin/env python3
"""
Create a benchmark dataset with HIGHLY CONSERVED (sharp) motifs.
Adaptive MCMC should excel here because it needs high α for fine-tuning.
"""

import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Settings for ~5 min run
NUM_SEQS = 500
SEQ_LEN = 200
NUM_MOTIFS = 2

# Two HIGHLY conserved motifs (95-100% at each position)
SHARP_MOTIFS = [
    "TATAAATA",    # TATA-like (very sharp)
    "GCGGCGGC",    # GC-rich (very sharp)
]

def sample_from_pwm(consensus, conservation=0.95):
    """Sample a sequence from a sharp PWM."""
    result = []
    bases = ['A', 'C', 'G', 'T']
    for base in consensus:
        if random.random() < conservation:
            result.append(base)
        else:
            # Pick random other base
            result.append(random.choice([b for b in bases if b != base]))
    return ''.join(result)

def generate_sequence(length, motifs, gc_bias=0.4):
    """Generate random sequence with planted sharp motifs."""
    bases = ['A', 'C', 'G', 'T']
    weights = [(1-gc_bias)/2, gc_bias/2, gc_bias/2, (1-gc_bias)/2]
    
    seq = []
    for _ in range(length):
        seq.append(random.choices(bases, weights)[0])
    seq = list(''.join(seq))
    
    # Plant each motif once with 95% conservation
    for motif in motifs:
        motif_len = len(motif)
        max_pos = length - motif_len - 10
        if max_pos > 10:
            pos = random.randint(10, max_pos)
            sampled = sample_from_pwm(motif, conservation=0.95)
            for i, base in enumerate(sampled):
                seq[pos + i] = base
    
    return ''.join(seq)

# Generate sequences
records = []
for i in range(NUM_SEQS):
    seq = generate_sequence(SEQ_LEN, SHARP_MOTIFS)
    record = SeqRecord(Seq(seq), id=f"seq_{i}", description="")
    records.append(record)

# Write FASTA
SeqIO.write(records, "sharp_benchmark.fa", "fasta")
print(f"Created sharp_benchmark.fa")
print(f"  Sequences: {NUM_SEQS}")
print(f"  Length: {SEQ_LEN}")
print(f"  Motifs: {SHARP_MOTIFS}")
print(f"  Conservation: 95%")
print(f"  Expected run time: ~5 min each")
