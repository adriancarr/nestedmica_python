"""
Higher-order Markov background model for motif discovery.
"""

import numpy as np
from typing import List, Any


def learn_markov_background(sequences: List[Any], order: int = 3) -> np.ndarray:
    """
    Learn k-th order Markov model from sequences.
    
    Args:
        sequences: List of BioPython sequence objects.
        order: Markov order (0-5). Default 3.
        
    Returns:
        np.ndarray: (4^order, 4) log2 probability matrix.
                    Row = context index (e.g., 64 for order-3)
                    Col = next base probability (A, C, G, T)
                    
    For order=0, returns (1, 4) matrix with base frequencies.
    """
    if order < 0 or order > 5:
        raise ValueError("Background order must be 0-5")
    
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num_contexts = 4 ** order if order > 0 else 1
    
    # Count matrix: counts[context_idx][next_base]
    counts = np.zeros((num_contexts, 4), dtype=np.float64)
    
    # Add pseudocounts to avoid zero probabilities
    counts += 0.01
    
    for seq_obj in sequences:
        seq = str(seq_obj.seq).upper()
        
        if order == 0:
            # Order-0: just count base frequencies
            for base in seq:
                if base in base_map:
                    counts[0, base_map[base]] += 1
        else:
            # Order-k: count (context, next_base) pairs
            for i in range(order, len(seq)):
                context = seq[i-order:i]
                next_base = seq[i]
                
                # Skip if any base is unknown
                if next_base not in base_map:
                    continue
                if not all(c in base_map for c in context):
                    continue
                
                # Compute context index
                ctx_idx = context_to_index(context, base_map)
                counts[ctx_idx, base_map[next_base]] += 1
    
    # Normalize rows to probabilities
    row_sums = counts.sum(axis=1, keepdims=True)
    probs = counts / row_sums
    
    # Convert to log2
    log_probs = np.log2(probs + 1e-10)
    
    return log_probs.astype(np.float64)


def context_to_index(context: str, base_map: dict) -> int:
    """
    Convert a context string (e.g., "CGT") to an integer index.
    
    Uses base-4 encoding: A=0, C=1, G=2, T=3
    "CGT" = 1*16 + 2*4 + 3*1 = 27
    """
    idx = 0
    for i, base in enumerate(context):
        idx = idx * 4 + base_map[base]
    return idx


def get_background_score_0th_order() -> np.ndarray:
    """
    Return uniform 0th-order background (for compatibility).
    Returns (1, 4) array with log2(0.25) = -2.0 for each base.
    """
    return np.full((1, 4), -2.0, dtype=np.float64)


def print_background_stats(bg_model: np.ndarray, order: int):
    """Print summary statistics of learned background model."""
    if order == 0:
        probs = np.power(2.0, bg_model[0])
        print(f"Background Model (Order-0):")
        print(f"  A: {probs[0]:.3f}, C: {probs[1]:.3f}, G: {probs[2]:.3f}, T: {probs[3]:.3f}")
        gc = probs[1] + probs[2]
        print(f"  GC Content: {gc*100:.1f}%")
    else:
        probs = np.power(2.0, bg_model)
        mean_probs = probs.mean(axis=0)
        print(f"Background Model (Order-{order}):")
        print(f"  Mean: A={mean_probs[0]:.3f}, C={mean_probs[1]:.3f}, G={mean_probs[2]:.3f}, T={mean_probs[3]:.3f}")
        gc = mean_probs[1] + mean_probs[2]
        print(f"  Mean GC: {gc*100:.1f}%")
        print(f"  Contexts: {bg_model.shape[0]}")
