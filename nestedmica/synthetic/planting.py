"""
Motif planting strategies for synthetic data.

This module provides strategies for placing motif instances into sequences,
including position selection, strand choice, and overlap avoidance.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class PlantedMotif:
    """Record of a planted motif instance."""
    motif_id: int
    start: int
    end: int
    strand: str  # '+' or '-'
    sequence: str  # The actual planted sequence
    gap_lengths: Optional[List[int]] = None  # For gapped motifs


class PlantingStrategy(ABC):
    """Abstract base class for motif planting strategies."""
    
    @abstractmethod
    def plant(self, sequence: str, motif_length: int, motif_id: int,
              rng: np.random.Generator,
              occupied: Optional[List[Tuple[int, int]]] = None) -> Tuple[str, Optional[PlantedMotif]]:
        """
        Plant a motif into a sequence.
        
        Args:
            sequence: Background sequence.
            motif_length: Length of motif to plant.
            motif_id: ID of this motif.
            rng: Random generator.
            
        Returns:
            Tuple of (modified_sequence, planted_info_or_None).
        """
        pass


class UniformPlanting(PlantingStrategy):
    """
    Plant motifs at uniformly random positions.
    
    Args:
        density: Probability that a sequence receives a motif (0-1).
        occurrences_per_seq: Number of motif instances per sequence.
        strand_bias: Probability of forward strand (0.5 = equal).
        min_edge_distance: Minimum distance from sequence edges.
        avoid_overlap: If True, prevent overlapping motif instances.
    """
    
    def __init__(self, density: float = 0.8, occurrences_per_seq: int = 1,
                 strand_bias: float = 0.5, min_edge_distance: int = 5,
                 avoid_overlap: bool = True):
        self.density = density
        self.occurrences_per_seq = occurrences_per_seq
        self.strand_bias = strand_bias
        self.min_edge_distance = min_edge_distance
        self.avoid_overlap = avoid_overlap
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq))
    
    def _find_valid_position(self, seq_length: int, motif_length: int,
                             occupied: List[Tuple[int, int]],
                             rng: np.random.Generator) -> Optional[int]:
        """Find a valid position that doesn't overlap with existing placements."""
        min_pos = self.min_edge_distance
        max_pos = seq_length - motif_length - self.min_edge_distance
        
        if max_pos < min_pos:
            return None
        
        # Try random positions
        for _ in range(50):  # Max attempts
            pos = rng.integers(min_pos, max_pos + 1)
            
            if not self.avoid_overlap:
                return pos
            
            # Check overlap
            overlap = False
            for occ_start, occ_end in occupied:
                if not (pos + motif_length <= occ_start or pos >= occ_end):
                    overlap = True
                    break
            
            if not overlap:
                return pos
        
        return None
    
    def plant(self, sequence: str, motif_length: int, motif_id: int,
              rng: np.random.Generator,
              occupied: Optional[List[Tuple[int, int]]] = None) -> Tuple[str, Optional[PlantedMotif]]:
        """Plant a single motif instance."""
        # Check if this sequence should receive a motif
        if rng.random() > self.density:
            return sequence, None
        
        if occupied is None:
            occupied = []
            
        pos = self._find_valid_position(len(sequence), motif_length, occupied, rng)
        if pos is None:
            return sequence, None
        
        # Choose strand
        strand = '+' if rng.random() < self.strand_bias else '-'
        
        # Create planted record (actual motif insertion done by caller)
        planted = PlantedMotif(
            motif_id=motif_id,
            start=pos,
            end=pos + motif_length,
            strand=strand,
            sequence=""  # Filled by caller
        )
        
        return sequence, planted
    
    def plant_multiple(self, sequence: str, motif_length: int, motif_id: int,
                       rng: np.random.Generator,
                       occupied: Optional[List[Tuple[int, int]]] = None) -> Tuple[str, List[PlantedMotif]]:
        """Plant multiple motif instances into one sequence."""
        # Check if this sequence should receive motifs
        if rng.random() > self.density:
            return sequence, []
        
        if occupied is None:
            occupied = []
        
        # We modify the occupied list in-place if provided, or use local one?
        # Ideally we append to the provided list so caller sees updates if they shared it.
        # But here occupied is local variable.
        # If caller passed occupied, we should append to it.
        # But this function returns planted_list and doesn't explicitly return new occupied state 
        # (unless it modifies the passed list).
        # We'll assume modification of passed list or tracking internally.
        # UniformPlanting.plant_multiple appends to `occupied` internally but that loop variable 
        # might shadow the argument if we aren't careful.
        
        # Wait, the method signature had `occupied`.
        # Existing logic:
        # occupied = []
        # We should use the argument.
        
        planted_list = []
        
        for _ in range(self.occurrences_per_seq):
            pos = self._find_valid_position(len(sequence), motif_length, occupied, rng)
            if pos is None:
                break
            
            # Choose strand
            strand = '+' if rng.random() < self.strand_bias else '-'
            
            planted = PlantedMotif(
                motif_id=motif_id,
                start=pos,
                end=pos + motif_length,
                strand=strand,
                sequence=""
            )
            planted_list.append(planted)
            occupied.append((pos, pos + motif_length))
        
        return sequence, planted_list


class ClusteredPlanting(PlantingStrategy):
    """
    Plant motifs in clusters (e.g., for testing combinatorial regulation).
    """
    
    def __init__(self, cluster_size: int = 3, cluster_spread: int = 50,
                 density: float = 0.5):
        self.cluster_size = cluster_size
        self.cluster_spread = cluster_spread
        self.density = density
    
    def plant(self, sequence: str, motif_length: int, motif_id: int,
              rng: np.random.Generator,
              occupied: Optional[List[Tuple[int, int]]] = None) -> Tuple[str, Optional[PlantedMotif]]:
        """Plant first motif of a cluster."""
        if rng.random() > self.density:
            return sequence, None
        
        # Choose cluster center
        min_pos = self.cluster_spread
        max_pos = len(sequence) - self.cluster_spread - motif_length
        
        if max_pos < min_pos:
            return sequence, None
        
        pos = rng.integers(min_pos, max_pos + 1)
        
        planted = PlantedMotif(
            motif_id=motif_id,
            start=pos,
            end=pos + motif_length,
            strand='+',
            sequence=""
        )
        
        return sequence, planted


class PositionalPlanting(PlantingStrategy):
    """
    Plant motifs with a positional bias (e.g., center of sequence).
    Typical for ChIP-seq peaks where the motif is central.
    """
    
    def __init__(self, distribution: str = 'normal', center_fraction: float = 0.5,
                 std_fraction: float = 0.2, density: float = 0.8,
                 strand_bias: float = 0.5):
        """
        Args:
            distribution: 'normal' or 'uniform' (localized).
            center_fraction: Center position as fraction of sequence length (0.0-1.0).
            std_fraction: Standard deviation as fraction of sequence length.
            density: Planting density.
            strand_bias: Forward strand bias.
        """
        self.distribution = distribution
        self.center_fraction = center_fraction
        self.std_fraction = std_fraction
        self.density = density
        self.strand_bias = strand_bias
        
    def plant(self, sequence: str, motif_length: int, motif_id: int,
              rng: np.random.Generator,
              occupied: Optional[List[Tuple[int, int]]] = None) -> Tuple[str, Optional[PlantedMotif]]:
        if rng.random() > self.density:
            return sequence, None
            
        seq_len = len(sequence)
        center = seq_len * self.center_fraction
        std = seq_len * self.std_fraction
        
        # Determine valid range
        min_pos = 0
        max_pos = seq_len - motif_length
        if max_pos < 0:
            return sequence, None
            
        # Sample position
        for _ in range(50): # Rejection sampling
            if self.distribution == 'normal':
                pos = rng.normal(center, std)
            else:
                pos = rng.uniform(center - std, center + std) # Approximate width
                
            pos = int(pos)
            if min_pos <= pos <= max_pos:
                # Check overlap
                overlap = False
                if occupied:
                    for occ_start, occ_end in occupied:
                        if not (pos + motif_length <= occ_start or pos >= occ_end):
                            overlap = True
                            break
                            
                if not overlap:
                    break
        else:
            return sequence, None # Fallback if standard deviation is too wide/narrow
            
        # Choose strand
        strand = '+' if rng.random() < self.strand_bias else '-'
        
        planted = PlantedMotif(
            motif_id=motif_id,
            start=pos,
            end=pos + motif_length,
            strand=strand,
            sequence=""
        )
        
        return sequence, planted
