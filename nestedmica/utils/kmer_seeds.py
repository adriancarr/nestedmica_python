"""
K-mer enrichment analysis for motif seed initialization.

This module provides MEME-style k-mer counting and enrichment scoring
to generate better initial PWM seeds instead of random Dirichlet samples.
"""

import numpy as np
from collections import Counter
from typing import List, Any, Tuple, Optional


def count_kmers(sequences: List[Any], k: int = 8) -> Counter:
    """
    Count all k-mers in sequences.
    
    Args:
        sequences: List of BioPython sequence objects.
        k: K-mer length.
        
    Returns:
        Counter: K-mer counts.
    """
    counts = Counter()
    
    for seq_obj in sequences:
        seq = str(seq_obj.seq).upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            # Skip k-mers containing N
            if 'N' not in kmer:
                counts[kmer] += 1
                # Also count reverse complement
                rc = reverse_complement(kmer)
                counts[rc] += 1
    
    return counts


def reverse_complement(kmer: str) -> str:
    """Return reverse complement of a k-mer."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, 'N') for base in reversed(kmer))


def calculate_expected_frequency(k: int, gc_content: float = 0.5) -> dict:
    """
    Calculate expected k-mer frequency under background model.
    
    Args:
        k: K-mer length.
        gc_content: GC content (0-1).
        
    Returns:
        dict: Expected frequency for each k-mer.
    """
    at_prob = (1 - gc_content) / 2  # A and T probability
    gc_prob = gc_content / 2        # G and C probability
    
    base_probs = {'A': at_prob, 'T': at_prob, 'C': gc_prob, 'G': gc_prob}
    
    def kmer_prob(kmer):
        prob = 1.0
        for base in kmer:
            prob *= base_probs.get(base, 0.25)
        return prob
    
    # Generate all k-mers
    from itertools import product
    bases = 'ACGT'
    expected = {}
    for kmer_tuple in product(bases, repeat=k):
        kmer = ''.join(kmer_tuple)
        expected[kmer] = kmer_prob(kmer)
    
    return expected


def find_enriched_kmers(sequences: List[Any], k: int = 8, 
                        top_n: int = 10, min_count: int = 5) -> List[Tuple[str, float]]:
    """
    Find k-mers enriched above background expectation.
    
    Args:
        sequences: List of BioPython sequence objects.
        k: K-mer length.
        top_n: Number of top enriched k-mers to return.
        min_count: Minimum count to consider.
        
    Returns:
        List of (kmer, enrichment_ratio) tuples, sorted by enrichment.
    """
    # Count k-mers
    counts = count_kmers(sequences, k)
    
    # Calculate total k-mer positions
    total_positions = sum(len(str(s.seq)) - k + 1 for s in sequences) * 2  # Both strands
    
    # Estimate GC content
    total_gc = 0
    total_bases = 0
    for seq_obj in sequences:
        seq = str(seq_obj.seq).upper()
        total_gc += seq.count('G') + seq.count('C')
        total_bases += len(seq.replace('N', ''))
    gc_content = total_gc / total_bases if total_bases > 0 else 0.5
    
    # Calculate expected frequencies
    expected = calculate_expected_frequency(k, gc_content)
    
    # Calculate enrichment ratios
    enrichment = []
    for kmer, count in counts.items():
        if count >= min_count:
            expected_count = expected.get(kmer, 1e-10) * total_positions
            ratio = count / max(expected_count, 1)
            enrichment.append((kmer, ratio, count))
    
    # Sort by enrichment ratio
    enrichment.sort(key=lambda x: -x[1])
    
    # Return top N (de-duplicate reverse complements)
    seen = set()
    result = []
    for kmer, ratio, count in enrichment:
        rc = reverse_complement(kmer)
        if kmer not in seen and rc not in seen:
            result.append((kmer, ratio))
            seen.add(kmer)
            seen.add(rc)
            if len(result) >= top_n:
                break
    
    return result


def kmer_to_pwm(kmer: str, pseudocount: float = 0.1) -> np.ndarray:
    """
    Convert a k-mer to a log2-space PWM.
    
    Creates a "soft" PWM where the consensus base has high probability
    and other bases have low probability.
    
    Args:
        kmer: K-mer string.
        pseudocount: Pseudocount for non-consensus bases.
        
    Returns:
        np.ndarray: Log2-space PWM (4, length).
    """
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    length = len(kmer)
    
    # Start with pseudocounts
    pwm = np.full((4, length), pseudocount)
    
    # Set high probability for consensus bases
    for pos, base in enumerate(kmer):
        if base in base_map:
            pwm[base_map[base], pos] = 1.0 - 3 * pseudocount
    
    # Normalize columns and convert to log2
    pwm = pwm / pwm.sum(axis=0, keepdims=True)
    log_pwm = np.log2(pwm + 1e-10)
    
    return log_pwm.astype(np.float64)


def generate_seed_pwms(sequences: List[Any], num_motifs: int, 
                       motif_length: int) -> List[np.ndarray]:
    """
    Generate seed PWMs from enriched k-mers.
    
    Args:
        sequences: List of BioPython sequence objects.
        num_motifs: Number of motifs to find.
        motif_length: Target motif length.
        
    Returns:
        List of log2-space PWM arrays (4, motif_length).
    """
    # Find enriched k-mers
    enriched = find_enriched_kmers(sequences, k=motif_length, top_n=num_motifs * 2)
    
    if not enriched:
        # Fallback: return None to signal random initialization
        return []
    
    # Convert top k-mers to PWMs
    seeds = []
    for kmer, ratio in enriched[:num_motifs]:
        pwm = kmer_to_pwm(kmer)
        seeds.append(pwm)
        print(f"  Seed {len(seeds)}: {kmer} (enrichment={ratio:.1f}x)")
    
    # If not enough enriched k-mers, pad with None (random init)
    while len(seeds) < num_motifs:
        seeds.append(None)
    
    return seeds
