"""
Gapped Motif Data Structures for NestedMICA.

This module provides classes for representing motifs with variable-length gaps
(hard gaps scored under background model). Supports multiple gaps per motif
with configurable constraints and prior distributions.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Literal
from enum import Enum


class GapPriorType(Enum):
    """Supported gap length prior distributions."""
    GEOMETRIC = "geometric"   # Favors shorter gaps (default)
    UNIFORM = "uniform"       # Equal probability across range
    POISSON = "poisson"       # Centered around expected length


@dataclass
class GapConfig:
    """
    User-configurable gap constraints and prior settings.
    
    Attributes:
        allow_gaps: Master switch for gapped motif discovery.
        min_gap_length: Minimum gap size (inclusive).
        max_gap_length: Maximum gap size (inclusive).
        max_num_gaps: Maximum number of gaps per motif.
        gap_prior: Prior distribution type for gap lengths.
        gap_prior_param: Parameter for the prior (p for geometric, λ for poisson).
    """
    allow_gaps: bool = False
    min_gap_length: int = 0
    max_gap_length: int = 20
    max_num_gaps: int = 2
    gap_prior: GapPriorType = GapPriorType.GEOMETRIC
    gap_prior_param: float = 0.1  # p for geometric, λ for poisson
    
    def __post_init__(self):
        if self.min_gap_length < 0:
            raise ValueError("min_gap_length must be >= 0")
        if self.max_gap_length < self.min_gap_length:
            raise ValueError("max_gap_length must be >= min_gap_length")
        if self.max_num_gaps < 0:
            raise ValueError("max_num_gaps must be >= 0")
        if not 0 < self.gap_prior_param < 1 and self.gap_prior == GapPriorType.GEOMETRIC:
            raise ValueError("geometric prior param p must be in (0, 1)")
    
    def log_prior_gap_length(self, gap_length: int) -> float:
        """
        Compute log prior probability for a gap length.
        
        Args:
            gap_length: The gap length to evaluate.
            
        Returns:
            Log probability (natural log).
        """
        if gap_length < self.min_gap_length or gap_length > self.max_gap_length:
            return -np.inf
        
        if self.gap_prior == GapPriorType.UNIFORM:
            # Uniform over [min, max]
            return -np.log(self.max_gap_length - self.min_gap_length + 1)
        
        elif self.gap_prior == GapPriorType.GEOMETRIC:
            # Geometric: P(g) ∝ (1-p)^g, normalized over [min, max]
            p = self.gap_prior_param
            # Unnormalized
            log_unnorm = gap_length * np.log(1 - p)
            # Normalization constant
            norm = sum((1 - p) ** g for g in range(self.min_gap_length, self.max_gap_length + 1))
            return log_unnorm - np.log(norm)
        
        elif self.gap_prior == GapPriorType.POISSON:
            # Poisson: P(g) ∝ λ^g / g!, normalized over [min, max]
            lam = self.gap_prior_param
            from math import lgamma
            log_unnorm = gap_length * np.log(lam) - lgamma(gap_length + 1)
            # Normalization
            norm = sum(np.exp(g * np.log(lam) - lgamma(g + 1)) 
                      for g in range(self.min_gap_length, self.max_gap_length + 1))
            return log_unnorm - np.log(norm)
        
        return 0.0
    
    def sample_gap_length(self, rng: np.random.Generator) -> int:
        """
        Sample a gap length from the prior distribution.
        
        Args:
            rng: NumPy random generator.
            
        Returns:
            Sampled gap length.
        """
        lengths = np.arange(self.min_gap_length, self.max_gap_length + 1)
        log_probs = np.array([self.log_prior_gap_length(g) for g in lengths])
        probs = np.exp(log_probs - np.max(log_probs))  # Numerical stability
        probs /= probs.sum()
        return int(rng.choice(lengths, p=probs))


@dataclass
class GapSpec:
    """
    Specifies a single gap within a gapped motif.
    
    Attributes:
        length: Gap length in base pairs.
    """
    length: int
    
    def __post_init__(self):
        if self.length < 0:
            raise ValueError("Gap length must be >= 0")


class GappedMotif:
    """
    Represents a gapped motif as alternating conserved blocks and gaps.
    
    Structure: [Block0][Gap0][Block1][Gap1]...[BlockN]
    
    Blocks are PWM arrays in log2 space.
    Gaps are hard gaps (scored under background model).
    
    Attributes:
        blocks: List of PWM arrays, each shape (4, L_i) in log2 space.
        gaps: List of GapSpec objects, len = len(blocks) - 1.
        config: GapConfig for constraints and priors.
    """
    
    def __init__(self, 
                 blocks: List[np.ndarray], 
                 gaps: Optional[List[GapSpec]] = None,
                 config: Optional[GapConfig] = None):
        """
        Initialize a gapped motif.
        
        Args:
            blocks: List of PWM arrays (4, L_i) in log2 space.
            gaps: List of GapSpec objects. If None, creates ungapped motif.
            config: Gap configuration. Uses defaults if None.
        """
        if len(blocks) == 0:
            raise ValueError("Must have at least one block")
        
        self.blocks = [b.astype(np.float64) for b in blocks]
        self.gaps = gaps if gaps is not None else []
        self.config = config if config is not None else GapConfig()
        
        # Validate structure
        if len(self.gaps) != len(self.blocks) - 1 and len(self.gaps) != 0:
            if len(self.blocks) > 1:
                raise ValueError(f"Need {len(self.blocks)-1} gaps for {len(self.blocks)} blocks, got {len(self.gaps)}")
        
        # Fill in missing gaps with zero-length gaps
        while len(self.gaps) < len(self.blocks) - 1:
            self.gaps.append(GapSpec(length=0))
    
    @property
    def num_blocks(self) -> int:
        """Number of conserved blocks."""
        return len(self.blocks)
    
    @property
    def num_gaps(self) -> int:
        """Number of gaps (non-zero length)."""
        return sum(1 for g in self.gaps if g.length > 0)
    
    @property
    def total_block_length(self) -> int:
        """Total length of all conserved blocks."""
        return sum(b.shape[1] for b in self.blocks)
    
    @property
    def total_gap_length(self) -> int:
        """Total length of all gaps."""
        return sum(g.length for g in self.gaps)
    
    @property
    def total_span(self) -> int:
        """Total span in bp (blocks + gaps)."""
        return self.total_block_length + self.total_gap_length
    
    @property
    def is_gapped(self) -> bool:
        """True if motif has any non-zero gaps."""
        return self.num_gaps > 0
    
    def get_block_lengths(self) -> np.ndarray:
        """Return array of block lengths."""
        return np.array([b.shape[1] for b in self.blocks], dtype=np.int64)
    
    def get_gap_lengths(self) -> np.ndarray:
        """Return array of gap lengths."""
        return np.array([g.length for g in self.gaps], dtype=np.int64)
    
    def to_flat_columns(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Flatten blocks for Cython processing.
        
        Returns:
            Tuple of (all_columns, block_offsets, block_lengths)
            - all_columns: (4, total_block_length) concatenated PWM columns
            - block_offsets: (num_blocks,) start offset for each block
            - block_lengths: (num_blocks,) length of each block
        """
        all_columns = np.concatenate(self.blocks, axis=1)
        block_lengths = self.get_block_lengths()
        block_offsets = np.zeros(len(self.blocks), dtype=np.int64)
        offset = 0
        for i, length in enumerate(block_lengths):
            block_offsets[i] = offset
            offset += length
        return all_columns, block_offsets, block_lengths
    
    def log_prior(self) -> float:
        """
        Compute log prior probability of this gapped motif structure.
        
        Includes prior on:
        - Gap lengths
        - Number of gaps (uniform over [0, max_num_gaps])
        
        Returns:
            Log prior probability.
        """
        log_p = 0.0
        
        # Prior on number of gaps
        actual_gaps = self.num_gaps
        if actual_gaps > self.config.max_num_gaps:
            return -np.inf
        log_p += -np.log(self.config.max_num_gaps + 1)  # Uniform over 0..max
        
        # Prior on each gap length
        for gap in self.gaps:
            if gap.length > 0:
                log_p += self.config.log_prior_gap_length(gap.length)
        
        return log_p
    
    def get_consensus(self) -> str:
        """
        Get consensus string with gap notation.
        
        Returns:
            String like "ACGT[5]TGCA" where [5] indicates 5bp gap.
        """
        bases = ['A', 'C', 'G', 'T']
        result = []
        
        for i, block in enumerate(self.blocks):
            # Consensus for this block
            for pos in range(block.shape[1]):
                max_idx = np.argmax(block[:, pos])
                result.append(bases[max_idx])
            
            # Add gap notation if not last block
            if i < len(self.gaps) and self.gaps[i].length > 0:
                result.append(f"[{self.gaps[i].length}]")
        
        return "".join(result)
    
    def copy(self) -> 'GappedMotif':
        """Create a deep copy of this motif."""
        return GappedMotif(
            blocks=[b.copy() for b in self.blocks],
            gaps=[GapSpec(length=g.length) for g in self.gaps],
            config=self.config
        )
    
    @classmethod
    def from_contiguous(cls, columns: np.ndarray, 
                        config: Optional[GapConfig] = None) -> 'GappedMotif':
        """
        Create a GappedMotif from a contiguous PWM (single block, no gaps).
        
        Args:
            columns: (4, L) PWM in log2 space.
            config: Gap configuration.
            
        Returns:
            GappedMotif with single block.
        """
        return cls(blocks=[columns], gaps=[], config=config)
    
    def __repr__(self) -> str:
        return (f"GappedMotif(blocks={self.num_blocks}, "
                f"gaps={self.num_gaps}, span={self.total_span}bp, "
                f"consensus='{self.get_consensus()}')")


# =============================================================================
# RJMCMC Proposal Functions for Gapped Motifs
# =============================================================================

def propose_insert_gap(motif: GappedMotif, rng: np.random.Generator,
                       min_block_size: int = 4) -> Tuple[Optional[GappedMotif], float]:
    """
    Propose inserting a new gap by splitting an existing block.
    
    This is a trans-dimensional move: increases number of blocks by 1.
    
    Args:
        motif: Current gapped motif.
        rng: Random generator.
        min_block_size: Minimum size for resulting blocks.
        
    Returns:
        Tuple of (proposed_motif, log_hastings_ratio).
        Returns (None, 0) if proposal is invalid.
    """
    config = motif.config
    
    # Check if we can add more gaps
    if motif.num_gaps >= config.max_num_gaps:
        return None, 0.0
    
    # Find blocks large enough to split
    splittable = []
    for i, block in enumerate(motif.blocks):
        block_len = block.shape[1]
        if block_len >= 2 * min_block_size:
            splittable.append(i)
    
    if not splittable:
        return None, 0.0
    
    # Choose block to split
    block_idx = rng.choice(splittable)
    block = motif.blocks[block_idx]
    block_len = block.shape[1]
    
    # Choose split position (must leave min_block_size on each side)
    min_split = min_block_size
    max_split = block_len - min_block_size
    if max_split < min_split:
        return None, 0.0
    split_pos = rng.integers(min_split, max_split + 1)
    
    # Sample initial gap length
    gap_length = config.sample_gap_length(rng)
    
    # Create new blocks
    block_a = block[:, :split_pos].copy()
    block_b = block[:, split_pos:].copy()
    
    # Build new motif structure
    new_blocks = motif.blocks[:block_idx] + [block_a, block_b] + motif.blocks[block_idx + 1:]
    new_gaps = list(motif.gaps)
    new_gaps.insert(block_idx, GapSpec(length=gap_length))
    
    proposed = GappedMotif(blocks=new_blocks, gaps=new_gaps, config=config)
    
    # Compute Hastings ratio for RJMCMC
    # Forward: P(choose block) * P(choose split) * P(gap length)
    # Reverse: P(choose gap to delete) in the proposed model
    num_splittable = len(splittable)
    num_split_positions = max_split - min_split + 1
    
    # Forward proposal probability
    log_q_forward = -np.log(num_splittable) - np.log(num_split_positions)
    log_q_forward += config.log_prior_gap_length(gap_length)
    
    # Reverse proposal: choose which gap to delete
    num_deletable_gaps = proposed.num_gaps
    log_q_reverse = -np.log(num_deletable_gaps) if num_deletable_gaps > 0 else 0.0
    
    log_hastings = log_q_reverse - log_q_forward
    
    return proposed, log_hastings


def propose_delete_gap(motif: GappedMotif, rng: np.random.Generator,
                       min_block_size: int = 4) -> Tuple[Optional[GappedMotif], float]:
    """
    Propose deleting a gap by merging adjacent blocks.
    
    This is a trans-dimensional move: decreases number of blocks by 1.
    
    Args:
        motif: Current gapped motif.
        rng: Random generator.
        min_block_size: Used for computing reverse Hastings ratio.
        
    Returns:
        Tuple of (proposed_motif, log_hastings_ratio).
        Returns (None, 0) if proposal is invalid.
    """
    if motif.num_gaps == 0:
        return None, 0.0
    
    config = motif.config
    
    # Find gaps with non-zero length that can be deleted
    deletable = [i for i, gap in enumerate(motif.gaps) if gap.length > 0]
    if not deletable:
        return None, 0.0
    
    # Choose gap to delete
    gap_idx = rng.choice(deletable)
    deleted_gap_length = motif.gaps[gap_idx].length
    
    # Merge the two adjacent blocks
    block_a = motif.blocks[gap_idx]
    block_b = motif.blocks[gap_idx + 1]
    merged = np.concatenate([block_a, block_b], axis=1)
    
    # Build new structure
    new_blocks = motif.blocks[:gap_idx] + [merged] + motif.blocks[gap_idx + 2:]
    new_gaps = motif.gaps[:gap_idx] + motif.gaps[gap_idx + 1:]
    
    proposed = GappedMotif(blocks=new_blocks, gaps=new_gaps, config=config)
    
    # Compute Hastings ratio (reverse of insertion)
    num_deletable = len(deletable)
    log_q_forward = -np.log(num_deletable)
    
    # Reverse: need to compute probability of inserting at this position
    # This requires knowing how many splittable blocks there would be
    merged_len = merged.shape[1]
    can_split = merged_len >= 2 * min_block_size
    if not can_split:
        # Can't reverse this move - invalid
        return None, 0.0
    
    split_positions = merged_len - 2 * min_block_size + 1
    log_q_reverse = -np.log(proposed.num_blocks)  # Approximate
    log_q_reverse += -np.log(split_positions)
    log_q_reverse += config.log_prior_gap_length(deleted_gap_length)
    
    log_hastings = log_q_reverse - log_q_forward
    
    return proposed, log_hastings


def propose_perturb_gap_length(motif: GappedMotif, rng: np.random.Generator,
                                delta: int = 3) -> Tuple[Optional[GappedMotif], float]:
    """
    Propose changing the length of an existing gap.
    
    This is a dimension-preserving move with symmetric proposal.
    
    Args:
        motif: Current gapped motif.
        rng: Random generator.
        delta: Maximum change in gap length.
        
    Returns:
        Tuple of (proposed_motif, log_hastings_ratio).
        Returns (None, 0) if proposal is invalid.
    """
    if motif.num_gaps == 0:
        return None, 0.0
    
    config = motif.config
    
    # Find gaps to perturb
    gap_indices = [i for i, gap in enumerate(motif.gaps) if gap.length > 0]
    if not gap_indices:
        return None, 0.0
    
    # Choose gap
    gap_idx = rng.choice(gap_indices)
    old_length = motif.gaps[gap_idx].length
    
    # Propose new length (symmetric random walk)
    change = rng.integers(-delta, delta + 1)
    new_length = old_length + change
    
    # Check bounds
    if new_length < config.min_gap_length or new_length > config.max_gap_length:
        return None, 0.0
    
    # Create new motif
    new_gaps = [GapSpec(length=g.length) for g in motif.gaps]
    new_gaps[gap_idx] = GapSpec(length=new_length)
    
    proposed = GappedMotif(
        blocks=[b.copy() for b in motif.blocks],
        gaps=new_gaps,
        config=config
    )
    
    # Symmetric proposal -> Hastings ratio = 1 (log = 0)
    return proposed, 0.0


def propose_shift_block_boundary(motif: GappedMotif, rng: np.random.Generator,
                                  delta: int = 2, min_block_size: int = 4) -> Tuple[Optional[GappedMotif], float]:
    """
    Propose shifting the boundary between a block and an adjacent gap.
    
    This moves columns from block to gap or vice versa.
    Note: Since gaps are "hard" (background model), this effectively
    changes what positions are considered conserved vs random.
    
    Args:
        motif: Current gapped motif.
        rng: Random generator.
        delta: Maximum shift amount.
        min_block_size: Minimum block size to maintain.
        
    Returns:
        Tuple of (proposed_motif, log_hastings_ratio).
        Returns (None, 0) if proposal is invalid.
    """
    if motif.num_gaps == 0:
        return None, 0.0
    
    config = motif.config
    
    # Find valid block-gap interfaces
    # Each gap has a left boundary (with preceding block) and right boundary
    interfaces = []
    for gap_idx, gap in enumerate(motif.gaps):
        if gap.length > 0:
            interfaces.append(('left', gap_idx))   # Left block's right edge
            interfaces.append(('right', gap_idx))  # Right block's left edge
    
    if not interfaces:
        return None, 0.0
    
    # Choose interface
    side, gap_idx = interfaces[rng.integers(len(interfaces))]
    
    # Propose shift direction and magnitude
    shift = rng.integers(-delta, delta + 1)
    if shift == 0:
        return None, 0.0  # No change
    
    # Copy current state
    new_blocks = [b.copy() for b in motif.blocks]
    new_gap_length = motif.gaps[gap_idx].length
    
    if side == 'left':
        # Shift affects block[gap_idx] and gap[gap_idx]
        block_idx = gap_idx
        block = new_blocks[block_idx]
        
        if shift > 0:
            # Move columns from block to gap (shrink block, expand gap)
            if block.shape[1] - shift < min_block_size:
                return None, 0.0
            new_blocks[block_idx] = block[:, :-shift]
            new_gap_length += shift
        else:
            # Move gap space to block (expand block, shrink gap)
            if new_gap_length + shift < config.min_gap_length:
                return None, 0.0
            # Need to generate new PWM columns (sample from prior)
            new_cols = np.log2(np.full((4, -shift), 0.25))  # Uniform prior
            new_blocks[block_idx] = np.concatenate([block, new_cols], axis=1)
            new_gap_length += shift  # shift is negative
    else:
        # side == 'right': affects gap[gap_idx] and block[gap_idx + 1]
        block_idx = gap_idx + 1
        block = new_blocks[block_idx]
        
        if shift > 0:
            # Move columns from block to gap
            if block.shape[1] - shift < min_block_size:
                return None, 0.0
            new_blocks[block_idx] = block[:, shift:]
            new_gap_length += shift
        else:
            # Move gap to block
            if new_gap_length + shift < config.min_gap_length:
                return None, 0.0
            new_cols = np.log2(np.full((4, -shift), 0.25))
            new_blocks[block_idx] = np.concatenate([new_cols, block], axis=1)
            new_gap_length += shift
    
    # Check gap length bounds
    if new_gap_length > config.max_gap_length:
        return None, 0.0
    
    # Create new gaps list
    new_gaps = [GapSpec(length=g.length) for g in motif.gaps]
    new_gaps[gap_idx] = GapSpec(length=new_gap_length)
    
    proposed = GappedMotif(blocks=new_blocks, gaps=new_gaps, config=config)
    
    # Symmetric proposal -> Hastings ratio = 1
    return proposed, 0.0

