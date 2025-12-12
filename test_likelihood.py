#!/usr/bin/env python3
"""
Unit test to verify likelihood calculation against manual calculation.
"""

import numpy as np
from nestedmica.maths import log2, addLog2
from nestedmica.model.motif import WeightMatrix, WeightedWeightMatrix
from nestedmica.model.core import SimpleContributionItem

# Simple test: 1 short sequence, 1 simple motif

def manual_likelihood(seq_indices, motif_log_probs, bg_score=-2.0, uncounted_exp=1.0):
    """
    Manually compute likelihood using the DP recurrence.
    
    seq_indices: array of base indices (0=A, 1=C, 2=G, 3=T)
    motif_log_probs: 2D array (L, 4) of log2 probabilities
    """
    seq_len = len(seq_indices)
    motif_len = len(motif_log_probs)
    
    # Transition probability
    trans_prob = 1.0 * uncounted_exp / seq_len
    trans_log = log2(trans_prob)
    base_penalty = log2(1.0 - trans_prob)
    
    print(f"trans_prob = {trans_prob:.6f}")
    print(f"trans_log = {trans_log:.2f}")
    print(f"base_penalty = {base_penalty:.6f}")
    
    # Pre-compute motif scores at each position
    motif_scores = []
    for pos in range(seq_len - motif_len + 1):
        score = 0.0
        for col in range(motif_len):
            base_idx = seq_indices[pos + col]
            score += motif_log_probs[col, base_idx]
        motif_scores.append(score)
    motif_scores = np.array(motif_scores)
    
    print(f"\nMotif scores at each position:")
    for i, s in enumerate(motif_scores):
        print(f"  pos {i}: {s:.2f}")
    
    # DP
    matrix = np.zeros(seq_len + 1)
    matrix[0] = 0.0
    
    print(f"\nDP iteration:")
    for i in range(1, seq_len + 1):
        # Background path
        bg_path = matrix[i-1] + bg_score + base_penalty
        
        # Motif path
        if i >= motif_len:
            motif_emit = motif_scores[i - motif_len]
            motif_path = matrix[i - motif_len] + motif_emit + trans_log
            score = addLog2(bg_path, motif_path)
            print(f"  i={i}: bg_path={bg_path:.2f}, motif_path={motif_path:.2f}, score={score:.2f}")
        else:
            score = bg_path
            print(f"  i={i}: bg_path={bg_path:.2f}, score={score:.2f}")
        
        matrix[i] = score
    
    return matrix[seq_len]

def test_simple():
    # Simple sequence: ACGTACGT (length 8, no N)
    # The motif should match perfectly at position 0 and 4
    seq = "ACGTACGT"
    seq_indices = []
    for c in seq:
        if c == 'A': seq_indices.append(0)
        elif c == 'C': seq_indices.append(1)
        elif c == 'G': seq_indices.append(2)
        elif c == 'T': seq_indices.append(3)
        else: seq_indices.append(-1)  # Invalid
    seq_indices = np.array(seq_indices)
    
    print(f"Sequence: {seq}")
    print(f"Indices: {seq_indices}")
    
    # Simple motif: perfect ACGT pattern
    # High probability for A at pos 0, C at pos 1, etc.
    motif_probs = np.array([
        [0.9, 0.03, 0.03, 0.04],  # A
        [0.03, 0.9, 0.04, 0.03],  # C
        [0.03, 0.04, 0.9, 0.03],  # G
        [0.04, 0.03, 0.03, 0.9],  # T
    ])
    motif_log_probs = log2(motif_probs)
    
    print(f"\nMotif log probabilities:")
    for i, row in enumerate(motif_log_probs):
        print(f"  pos {i}: {row}")
    
    result = manual_likelihood(seq_indices, motif_log_probs)
    print(f"\nFinal likelihood (manual): {result:.2f}")
    
    # Compare with Python implementation
    print("\n--- Python Implementation ---")
    from nestedmica.model.motif_generative import MotifFacette, MotifUncountedLikelihood
    from nestedmica.model.core import SimpleContributionItem
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    facette = MotifFacette(uncounted_expectation=1.0)
    
    # Create datum with SAME sequence
    datum = SeqRecord(Seq(seq))
    
    # Create calculator
    calc = MotifUncountedLikelihood(facette, datum)
    
    # Create contribution
    wm = WeightMatrix(motif_log_probs)
    wwm = WeightedWeightMatrix(wm, 1.0)
    contrib = SimpleContributionItem(wwm)
    
    weights = np.array([1.0])  # Motif is active
    
    py_result = calc.likelihood([contrib], weights)
    print(f"Python likelihood: {py_result:.2f}")
    
    print(f"\nDifference: {abs(result - py_result):.2f}")

if __name__ == "__main__":
    test_simple()
