"""
Motif definitions for synthetic data generation.

This module provides classes for defining motifs to plant in synthetic
sequences, including support for importing from JASPAR, MEME, and other formats.
"""

import numpy as np
import re
from abc import ABC, abstractmethod
from typing import Optional, List, Tuple, Union
from dataclasses import dataclass


# IUPAC ambiguity codes
IUPAC_CODES = {
    'A': [1, 0, 0, 0],
    'C': [0, 1, 0, 0],
    'G': [0, 0, 1, 0],
    'T': [0, 0, 0, 1],
    'U': [0, 0, 0, 1],
    'R': [0.5, 0, 0.5, 0],      # A or G
    'Y': [0, 0.5, 0, 0.5],      # C or T
    'S': [0, 0.5, 0.5, 0],      # G or C
    'W': [0.5, 0, 0, 0.5],      # A or T
    'K': [0, 0, 0.5, 0.5],      # G or T
    'M': [0.5, 0.5, 0, 0],      # A or C
    'B': [0, 1/3, 1/3, 1/3],    # C, G, or T
    'D': [1/3, 0, 1/3, 1/3],    # A, G, or T
    'H': [1/3, 1/3, 0, 1/3],    # A, C, or T
    'V': [1/3, 1/3, 1/3, 0],    # A, C, or G
    'N': [0.25, 0.25, 0.25, 0.25],  # Any
}


class SyntheticMotif(ABC):
    """Abstract base class for motif definitions."""
    
    @abstractmethod
    def sample(self, rng: np.random.Generator) -> str:
        """Sample a motif instance."""
        pass
    
    @abstractmethod
    def get_pwm(self) -> np.ndarray:
        """Get PWM representation (4, L) in probability space."""
        pass
    
    @property
    @abstractmethod
    def length(self) -> int:
        """Motif length in base pairs."""
        pass
    
    @property
    def consensus(self) -> str:
        """Get consensus sequence."""
        bases = ['A', 'C', 'G', 'T']
        pwm = self.get_pwm()
        return ''.join(bases[np.argmax(pwm[:, i])] for i in range(pwm.shape[1]))


class ConsensusMotif(SyntheticMotif):
    """
    Motif defined by IUPAC consensus string.
    
    Examples:
        "ACGT" - exact match
        "ACGTNNN" - with wildcards
        "ACGTWSY" - with ambiguity codes
    """
    
    def __init__(self, consensus: str, noise: float = 0.05):
        """
        Create motif from consensus.
        
        Args:
            consensus: IUPAC consensus string.
            noise: Background probability for non-consensus bases.
        """
        self.consensus_str = consensus.upper()
        self.noise = noise
        self._pwm = self._consensus_to_pwm()
    
    def _consensus_to_pwm(self) -> np.ndarray:
        """Convert consensus to PWM."""
        pwm = np.zeros((4, len(self.consensus_str)))
        
        for i, char in enumerate(self.consensus_str):
            if char in IUPAC_CODES:
                probs = np.array(IUPAC_CODES[char])
                # Add noise
                probs = probs * (1 - self.noise) + self.noise * 0.25
                pwm[:, i] = probs
            else:
                # Unknown character - treat as N
                pwm[:, i] = 0.25
        
        return pwm
    
    def sample(self, rng: np.random.Generator) -> str:
        bases = ['A', 'C', 'G', 'T']
        result = []
        for i in range(self._pwm.shape[1]):
            idx = rng.choice(4, p=self._pwm[:, i])
            result.append(bases[idx])
        return ''.join(result)
    
    def get_pwm(self) -> np.ndarray:
        return self._pwm.copy()
    
    @property
    def length(self) -> int:
        return len(self.consensus_str)


class PWMMotif(SyntheticMotif):
    """
    Motif defined by a Position Weight Matrix (PWM).
    
    Supports import from JASPAR, MEME, and other formats.
    """
    
    def __init__(self, pwm: np.ndarray, name: str = "motif"):
        """
        Create motif from PWM.
        
        Args:
            pwm: (4, L) probability matrix (A, C, G, T order).
            name: Motif name.
        """
        if pwm.shape[0] != 4:
            raise ValueError("PWM must have shape (4, length)")
        
        self._pwm = pwm.astype(np.float64)
        self.name = name
        
        # Normalize columns to sum to 1
        col_sums = self._pwm.sum(axis=0, keepdims=True)
        self._pwm = self._pwm / col_sums
    
    def sample(self, rng: np.random.Generator) -> str:
        bases = ['A', 'C', 'G', 'T']
        result = []
        for i in range(self._pwm.shape[1]):
            idx = rng.choice(4, p=self._pwm[:, i])
            result.append(bases[idx])
        return ''.join(result)
    
    def get_pwm(self) -> np.ndarray:
        return self._pwm.copy()
    
    @property
    def length(self) -> int:
        return self._pwm.shape[1]
    
    @classmethod
    def from_counts(cls, counts: np.ndarray, pseudocount: float = 0.1,
                    name: str = "motif") -> 'PWMMotif':
        """
        Create PWM from count matrix.
        
        Args:
            counts: (4, L) count matrix.
            pseudocount: Pseudocount for regularization.
            name: Motif name.
        """
        pwm = counts + pseudocount
        pwm = pwm / pwm.sum(axis=0, keepdims=True)
        return cls(pwm, name)
    
    @classmethod
    def from_meme(cls, filepath: str) -> List['PWMMotif']:
        """
        Parse MEME minimal format file.
        
        Args:
            filepath: Path to MEME file.
            
        Returns:
            List of PWMMotif objects.
        """
        motifs = []
        current_name = None
        current_pwm = []
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('MOTIF'):
                    # Save previous motif
                    if current_name and current_pwm:
                        pwm = np.array(current_pwm).T  # (L, 4) -> (4, L)
                        motifs.append(cls(pwm, current_name))
                    
                    # Start new motif
                    parts = line.split()
                    current_name = parts[1] if len(parts) > 1 else f"motif_{len(motifs)}"
                    current_pwm = []
                
                elif line and not line.startswith('#') and not line.startswith('letter-probability'):
                    # Try to parse as PWM row
                    try:
                        values = [float(x) for x in line.split()[:4]]
                        if len(values) == 4:
                            current_pwm.append(values)
                    except ValueError:
                        pass
        
        # Save last motif
        if current_name and current_pwm:
            pwm = np.array(current_pwm).T
            motifs.append(cls(pwm, current_name))
        
        return motifs
    
    @classmethod
    def from_jaspar_file(cls, filepath: str) -> 'PWMMotif':
        """
        Parse JASPAR format file.
        
        Args:
            filepath: Path to JASPAR .pfm file.
            
        Returns:
            PWMMotif object.
        """
        counts = []
        name = "motif"
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name = line[1:].split()[0]
                elif line:
                    # Parse row: A [ 1 2 3 4 ] or just numbers
                    match = re.search(r'\[(.*?)\]', line)
                    if match:
                        values = [float(x) for x in match.group(1).split()]
                    else:
                        # Try direct number parsing
                        parts = line.split()
                        if parts[0] in 'ACGT':
                            values = [float(x) for x in parts[1:]]
                        else:
                            values = [float(x) for x in parts]
                    counts.append(values)
        
        if len(counts) == 4:
            return cls.from_counts(np.array(counts), name=name)
        else:
            raise ValueError(f"Expected 4 rows (A,C,G,T), got {len(counts)}")


class GappedSyntheticMotif(SyntheticMotif):
    """
    Gapped motif with blocks and variable-length gaps.
    
    Structure: Block1 - Gap1 - Block2 - Gap2 - ...
    """
    
    def __init__(self, blocks: List[Union[SyntheticMotif, str]],
                 gap_ranges: Optional[List[Tuple[int, int]]] = None):
        """
        Create gapped motif.
        
        Args:
            blocks: List of motifs or consensus strings.
            gap_ranges: List of (min, max) gap lengths between blocks.
        """
        self.blocks = []
        for block in blocks:
            if isinstance(block, str):
                self.blocks.append(ConsensusMotif(block))
            else:
                self.blocks.append(block)
        
        if gap_ranges is None:
            gap_ranges = [(0, 0)] * (len(self.blocks) - 1)
        
        if len(gap_ranges) != len(self.blocks) - 1:
            raise ValueError(f"Need {len(self.blocks)-1} gap ranges, got {len(gap_ranges)}")
        
        self.gap_ranges = gap_ranges
    
    def sample(self, rng: np.random.Generator) -> Tuple[str, List[int]]:
        """
        Sample a gapped motif instance.
        
        Returns:
            Tuple of (sequence, gap_lengths)
        """
        parts = []
        gap_lengths = []
        
        for i, block in enumerate(self.blocks):
            parts.append(block.sample(rng))
            if i < len(self.gap_ranges):
                gap_min, gap_max = self.gap_ranges[i]
                gap_len = rng.integers(gap_min, gap_max + 1)
                gap_lengths.append(gap_len)
                # Gap is filled with random bases
                gap_seq = ''.join(rng.choice(['A', 'C', 'G', 'T'], size=gap_len))
                parts.append(gap_seq)
        
        return ''.join(parts), gap_lengths
    
    def get_pwm(self) -> np.ndarray:
        """Get concatenated block PWMs (gaps as uniform)."""
        parts = []
        for i, block in enumerate(self.blocks):
            parts.append(block.get_pwm())
            if i < len(self.gap_ranges):
                gap_min, gap_max = self.gap_ranges[i]
                avg_gap = (gap_min + gap_max) // 2
                if avg_gap > 0:
                    gap_pwm = np.full((4, avg_gap), 0.25)
                    parts.append(gap_pwm)
        return np.concatenate(parts, axis=1)
    
    @property
    def length(self) -> int:
        """Minimum total length."""
        block_len = sum(b.length for b in self.blocks)
        gap_len = sum(g[0] for g in self.gap_ranges)
        return block_len + gap_len
    
    @classmethod
    def from_pattern(cls, pattern: str) -> 'GappedSyntheticMotif':
        """
        Parse pattern string like "ACGT[5-10]TGCA".
        
        Args:
            pattern: Pattern with [min-max] gap notation.
            
        Returns:
            GappedSyntheticMotif instance.
        """
        # Split on gap patterns
        gap_pattern = r'\[(\d+)-(\d+)\]'
        
        blocks = []
        gap_ranges = []
        remaining = pattern
        
        while remaining:
            match = re.search(gap_pattern, remaining)
            if match:
                block_seq = remaining[:match.start()]
                if block_seq:
                    blocks.append(block_seq)
                gap_min = int(match.group(1))
                gap_max = int(match.group(2))
                gap_ranges.append((gap_min, gap_max))
                remaining = remaining[match.end():]
            else:
                if remaining:
                    blocks.append(remaining)
                break
        
        return cls(blocks, gap_ranges)
