"""
Background sequence generators for synthetic data.

This module provides classes for generating random DNA background sequences
with configurable properties like GC content, Markov order, and learned
distributions from real sequence data.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Optional, List, Any, Union
from dataclasses import dataclass
from Bio import SeqIO


# =============================================================================
# Length Distribution Classes
# =============================================================================

class LengthDistribution(ABC):
    """Abstract base class for sequence length distributions."""
    
    @abstractmethod
    def sample(self, rng: np.random.Generator) -> int:
        """Sample a sequence length."""
        pass
    
    @abstractmethod
    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Sample n sequence lengths."""
        pass


@dataclass
class FixedLength(LengthDistribution):
    """Fixed length for all sequences."""
    length: int
    
    def sample(self, rng: np.random.Generator) -> int:
        return self.length
    
    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return np.full(n, self.length, dtype=np.int64)


@dataclass
class NormalLength(LengthDistribution):
    """Normal distribution for sequence lengths (typical for ChIP-seq)."""
    mean: float
    std: float
    min_len: int = 50
    max_len: int = 1000
    
    def sample(self, rng: np.random.Generator) -> int:
        length = int(rng.normal(self.mean, self.std))
        return max(self.min_len, min(self.max_len, length))
    
    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        lengths = rng.normal(self.mean, self.std, size=n).astype(np.int64)
        return np.clip(lengths, self.min_len, self.max_len)


@dataclass
class LogNormalLength(LengthDistribution):
    """Log-normal distribution (common for fragment sizes)."""
    mean: float  # Mean of log(length)
    sigma: float  # Std of log(length)
    min_len: int = 50
    max_len: int = 2000
    
    def sample(self, rng: np.random.Generator) -> int:
        length = int(rng.lognormal(self.mean, self.sigma))
        return max(self.min_len, min(self.max_len, length))
    
    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        lengths = rng.lognormal(self.mean, self.sigma, size=n).astype(np.int64)
        return np.clip(lengths, self.min_len, self.max_len)


class LearnedLength(LengthDistribution):
    """Learn length distribution from real FASTA data."""
    
    def __init__(self, lengths: np.ndarray):
        """
        Initialize with observed lengths.
        
        Args:
            lengths: Array of observed sequence lengths.
        """
        self.lengths = lengths
        self.min_len = int(lengths.min())
        self.max_len = int(lengths.max())
        self.mean = float(lengths.mean())
        self.std = float(lengths.std())
    
    @classmethod
    def from_fasta(cls, fasta_path: str) -> 'LearnedLength':
        """
        Learn length distribution from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file.
            
        Returns:
            LearnedLength instance.
        """
        lengths = []
        for record in SeqIO.parse(fasta_path, 'fasta'):
            lengths.append(len(record.seq))
        return cls(np.array(lengths, dtype=np.int64))
    
    def sample(self, rng: np.random.Generator) -> int:
        return int(rng.choice(self.lengths))
    
    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        return rng.choice(self.lengths, size=n, replace=True)


# =============================================================================
# Background Generator Classes
# =============================================================================

class BackgroundGenerator(ABC):
    """Abstract base class for background sequence generators."""
    
    @abstractmethod
    def generate(self, length: int, rng: np.random.Generator) -> str:
        """Generate a background sequence of given length."""
        pass
    
    def generate_n(self, lengths: Union[int, np.ndarray, LengthDistribution], 
                   n: int, rng: np.random.Generator) -> List[str]:
        """
        Generate n background sequences.
        
        Args:
            lengths: Fixed length, array of lengths, or LengthDistribution.
            n: Number of sequences.
            rng: Random generator.
            
        Returns:
            List of sequences.
        """
        if isinstance(lengths, int):
            lens = np.full(n, lengths)
        elif isinstance(lengths, LengthDistribution):
            lens = lengths.sample_n(n, rng)
        else:
            lens = lengths
        
        return [self.generate(int(l), rng) for l in lens]


class UniformBackground(BackgroundGenerator):
    """
    Generate IID sequences with configurable base composition.
    
    Args:
        gc_content: GC content as fraction (0.0 to 1.0).
    """
    
    def __init__(self, gc_content: float = 0.5):
        if not 0.0 <= gc_content <= 1.0:
            raise ValueError("gc_content must be between 0 and 1")
        
        self.gc_content = gc_content
        at = (1.0 - gc_content) / 2.0
        gc = gc_content / 2.0
        self.probs = np.array([at, gc, gc, at])  # A, C, G, T
        self.bases = np.array(['A', 'C', 'G', 'T'])
    
    def generate(self, length: int, rng: np.random.Generator) -> str:
        indices = rng.choice(4, size=length, p=self.probs)
        return ''.join(self.bases[indices])


class MarkovBackground(BackgroundGenerator):
    """
    Generate sequences from a k-th order Markov model.
    
    Args:
        order: Markov order (1-5).
        gc_content: Target GC content for transition matrix.
        transition_matrix: Optional pre-computed transition matrix.
    """
    
    def __init__(self, order: int = 3, gc_content: float = 0.5,
                 transition_matrix: Optional[np.ndarray] = None):
        if order < 1 or order > 5:
            raise ValueError("order must be 1-5")
        
        self.order = order
        self.gc_content = gc_content
        self.bases = np.array(['A', 'C', 'G', 'T'])
        self.base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        if transition_matrix is not None:
            self.transition = transition_matrix
        else:
            self.transition = self._create_gc_biased_transition(gc_content)
    
    def _create_gc_biased_transition(self, gc: float) -> np.ndarray:
        """Create transition matrix with target GC content."""
        num_contexts = 4 ** self.order
        at = (1.0 - gc) / 2.0
        gc_prob = gc / 2.0
        
        # Simple: same distribution for all contexts
        base_probs = np.array([at, gc_prob, gc_prob, at])
        transition = np.tile(base_probs, (num_contexts, 1))
        return transition
    
    def _context_to_idx(self, context: str) -> int:
        """Convert context string to integer index."""
        idx = 0
        for base in context:
            idx = idx * 4 + self.base_map[base]
        return idx
    
    def generate(self, length: int, rng: np.random.Generator) -> str:
        if length == 0:
            return ""
            
        # Try using fast Cython implementation
        try:
            from nestedmica.model.cython_model import generate_markov_sequence
            
            # Pre-generate random numbers
            rand_vals = rng.random(size=length).astype(np.float64)
            
            return generate_markov_sequence(
                length, 
                self.transition.astype(np.float64), 
                self._get_init_probs().astype(np.float64), 
                self.order, 
                rand_vals
            )
        except (ImportError, AttributeError):
            # Fallback to slow Python
            pass
        
        # Generate initial context with uniform IID
        at = (1.0 - self.gc_content) / 2.0
        gc = self.gc_content / 2.0
        init_probs = np.array([at, gc, gc, at])
        
        result = []
        for _ in range(min(self.order, length)):
            idx = rng.choice(4, p=init_probs)
            result.append(self.bases[idx])
        
        # Generate rest using Markov model
        for i in range(self.order, length):
            context = ''.join(result[i - self.order:i])
            ctx_idx = self._context_to_idx(context)
            probs = self.transition[ctx_idx]
            idx = rng.choice(4, p=probs)
            result.append(self.bases[idx])
        
        return ''.join(result)

    def _get_init_probs(self) -> np.ndarray:
        at = (1.0 - self.gc_content) / 2.0
        gc = self.gc_content / 2.0
        return np.array([at, gc, gc, at])


class LearnedBackground(BackgroundGenerator):
    """
    Generate sequences matching properties learned from real data.
    
    Learn Markov model from a FASTA file of real sequences
    (e.g., Drosophila promoters, intergenic regions).
    """
    
    def __init__(self, transition_matrix: np.ndarray, order: int,
                 init_probs: np.ndarray):
        """
        Initialize with learned parameters.
        
        Args:
            transition_matrix: (4^order, 4) probability matrix.
            order: Markov order.
            init_probs: (4,) initial base probabilities.
        """
        self.transition = transition_matrix
        self.order = order
        self.init_probs = init_probs
        self.bases = np.array(['A', 'C', 'G', 'T'])
        self.base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    @classmethod
    def from_fasta(cls, fasta_path: str, order: int = 3) -> 'LearnedBackground':
        """
        Learn background model from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file.
            order: Markov order to learn (1-5).
            
        Returns:
            LearnedBackground instance.
        """
        base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_contexts = 4 ** order if order > 0 else 1
        
        # Count transition frequencies
        counts = np.zeros((num_contexts, 4), dtype=np.float64)
        init_counts = np.zeros(4, dtype=np.float64)
        counts += 0.01  # Pseudocounts
        init_counts += 0.01
        
        for record in SeqIO.parse(fasta_path, 'fasta'):
            seq = str(record.seq).upper()
            
            # Count initial bases
            for i in range(min(order, len(seq))):
                if seq[i] in base_map:
                    init_counts[base_map[seq[i]]] += 1
            
            # Count transitions
            for i in range(order, len(seq)):
                context = seq[i-order:i]
                next_base = seq[i]
                
                if next_base not in base_map:
                    continue
                if not all(c in base_map for c in context):
                    continue
                
                # Compute context index
                ctx_idx = 0
                for c in context:
                    ctx_idx = ctx_idx * 4 + base_map[c]
                
                counts[ctx_idx, base_map[next_base]] += 1
        
        # Normalize to probabilities
        transition = counts / counts.sum(axis=1, keepdims=True)
        init_probs = init_counts / init_counts.sum()
        
        return cls(transition, order, init_probs)
    
    def _context_to_idx(self, context: str) -> int:
        """Convert context string to integer index."""
        idx = 0
        for base in context:
            idx = idx * 4 + self.base_map[base]
        return idx
    
    def generate(self, length: int, rng: np.random.Generator) -> str:
        if length == 0:
            return ""
            
        # Try using fast Cython implementation
        try:
            from nestedmica.model.cython_model import generate_markov_sequence
            
            # Pre-generate random numbers
            rand_vals = rng.random(size=length).astype(np.float64)
            
            return generate_markov_sequence(
                length, 
                self.transition.astype(np.float64), 
                self.init_probs.astype(np.float64), 
                self.order, 
                rand_vals
            )
        except (ImportError, AttributeError):
            pass
        
        result = []
        
        # Generate initial context
        for _ in range(min(self.order, length)):
            idx = rng.choice(4, p=self.init_probs)
            result.append(self.bases[idx])
        
        # Generate rest using learned Markov model
        for i in range(self.order, length):
            context = ''.join(result[i - self.order:i])
            ctx_idx = self._context_to_idx(context)
            probs = self.transition[ctx_idx]
            idx = rng.choice(4, p=probs)
            result.append(self.bases[idx])
        
        return ''.join(result)
    
    @property
    def gc_content(self) -> float:
        """Estimated GC content from learned model."""
        # Average over all contexts
        avg_probs = self.transition.mean(axis=0)
        return float(avg_probs[1] + avg_probs[2])  # C + G


# =============================================================================
# Shuffling Utilities
# =============================================================================

def dinucl_shuffle(seq: str, rng: np.random.Generator) -> str:
    """
    Perform a dinucleotide shuffle of the sequence using the Altschul-Erickson algorithm.
    Preserves exact dinucleotide counts.
    """
    if len(seq) < 2:
        return seq
    
    seq = seq.upper()
    bases = ['A', 'C', 'G', 'T']
    base_map = {b: i for i, b in enumerate(bases)}
    
    # transform seq to indices
    s = [base_map.get(c, 0) for c in seq]  # separate unknown as 0 or handle?
    # For now assume ACGT. If others, map to A
    
    n = len(s)
    last_char = s[-1]
    
    # 1. Count edges (transitions)
    counts = np.zeros((4, 4), dtype=np.int64)
    for i in range(n - 1):
        counts[s[i], s[i+1]] += 1
        
    # 2. For each vertex, create a list of outgoing edges
    lists = [[] for _ in range(4)]
    for i in range(n - 1):
        lists[s[i]].append(s[i+1])
    
    # 3. Shuffle the lists
    for u in range(4):
        rng.shuffle(lists[u])
        
    # 4. Construct the graph for Eulerian path
    degrees = counts.sum(axis=1)
    
    # Vertices with non-zero degree are the ones we care about
    nodes = [i for i in range(4) if degrees[i] > 0]
    if not nodes:
        return seq
        
    # Altschul-Erickson algorithm:
    # Ensure the "last edges" (one per node except last_char) form a tree rooted at last_char.
    root = last_char
    present_nodes = set(s)
    if len(present_nodes) == 1:
        return seq
    
    # Rejection sampling for the tree condition (efficient for N=4)
    while True:
        valid_tree = True
        
        # We need to check connectivity for all nodes that have outgoing edges.
        active_nodes = [u for u in range(4) if len(lists[u]) > 0]
        
        # Temporarily look at the last item of each list as the "tree edge"
        tree_edges = {}
        for u in active_nodes:
            if u == root:
                continue
            if not lists[u]:
                valid_tree = False
                break
            tree_edges[u] = lists[u][-1]
        
        if valid_tree:
            # Check if all active nodes reach root without cycles
            for u in active_nodes:
                if u == root:
                    continue
                curr = u
                path = set()
                while curr != root:
                    if curr in path: # Cycle detected
                        valid_tree = False
                        break
                    path.add(curr)
                    if curr not in tree_edges:
                        valid_tree = False
                        break
                    curr = tree_edges[curr]
                if not valid_tree:
                    break
        
        if valid_tree:
            break
            
        # If not valid, permute and try again
        for u in range(4):
            rng.shuffle(lists[u])

    # Construct the path
    curr = s[0]
    res_indices = [curr]
    consumed = [0] * 4
    
    for _ in range(n - 1):
        # Determine which edge to take
        u = curr
        # We consume from index 0 upwards. The tree edges are at the END of the lists.
        # We implicitly reserve them for last.
        next_v = lists[u][consumed[u]]
        consumed[u] += 1
        
        res_indices.append(next_v)
        curr = next_v
        
    return ''.join(bases[i] for i in res_indices)


class ShuffledBackground(BackgroundGenerator):
    """
    Background generator that produces dinucleotide-shuffled versions
    of provided source sequences.
    """
    
    def __init__(self, sequences: List[str]):
        """
        Initialize with source sequences.
        
        Args:
            sequences: List of DNA sequences to use as shuffle templates.
        """
        self.sequences = sequences
        self.lengths = {}
        for s in sequences:
            l = len(s)
            if l not in self.lengths:
                self.lengths[l] = []
            self.lengths[l].append(s)
            
    @classmethod
    def from_fasta(cls, fasta_path: str) -> 'ShuffledBackground':
        """Load sequences from FASTA."""
        seqs = []
        for record in SeqIO.parse(fasta_path, 'fasta'):
            seqs.append(str(record.seq).upper())
        return cls(seqs)
            
    def generate(self, length: int, rng: np.random.Generator) -> str:
        """
        Generate a shuffled sequence.
        
        If a source sequence of the requested length exists, use it.
        Otherwise, pick a random source sequence and trim/pad.
        """
        if length in self.lengths:
            template = rng.choice(self.lengths[length])
        else:
            # Fallback: pick random
            template = rng.choice(self.sequences)
            
        shuffled = dinucl_shuffle(template, rng)
        
        # Adjust length if needed
        if len(shuffled) == length:
            return shuffled
        elif len(shuffled) > length:
            return shuffled[:length]
        else:
            # Random padding
            pad_len = length - len(shuffled)
            pad = ''.join(rng.choice(['A', 'C', 'G', 'T'], size=pad_len))
            return shuffled + pad
