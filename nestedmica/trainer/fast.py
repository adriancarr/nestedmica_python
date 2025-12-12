"""
Core Nested MICA Trainer implementation (Cython-Accelerated).
"""

import numpy as np
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from typing import List, Tuple, Dict, Any, Optional

# Import Cython-optimized modules
from nestedmica.model.cython_model import CythonWeightMatrix, batch_likelihood

def sample_dirichlet_pwm(length: int, alpha: float = 1.0) -> np.ndarray:
    """
    Sample a random PWM from Dirichlet prior.
    
    Args:
        length (int): Length of the PWM.
        alpha (float): Dirichlet concentration parameter.
        
    Returns:
        np.ndarray: Log2-space probability matrix (4, length).
    """
    columns = np.random.dirichlet([alpha] * 4, size=length).T
    return np.log2(columns + 1e-10).astype(np.float64)


class FastTrainer:
    """
    Optimized Nested Sampling trainer with Data-Oriented Threaded Parallelism.
    
    Attributes:
        sequences: List of input sequences.
        num_motifs (int): Number of motifs to discover.
        motif_length (int): Initial motif length (can change with Indel/Zap).
        ensemble_size (int): Number of particles in Nested Sampling ensemble.
        n_jobs (int): Number of threads.
    """
    
    def __init__(self, sequences: List[Any], num_motifs: int, motif_length: int, 
                 ensemble_size: int, n_jobs: int = -1):
        """
        Initialize the trainer.

        Args:
            sequences (List): BioPython sequence objects.
            num_motifs (int): Number of motifs.
            motif_length (int): Initial target length.
            ensemble_size (int): Population size.
            n_jobs (int): Threads (-1 for all cores).
        """
        self.sequences = sequences
        self.num_motifs = num_motifs
        self.motif_length = motif_length
        self.ensemble_size = ensemble_size
        self.num_sequences = len(sequences)
        
        # Number of parallel threads
        self.n_jobs = n_jobs if n_jobs > 0 else multiprocessing.cpu_count()
        
        # Pre-process sequences into flat int array
        max_len = 0
        seq_indices = []
        char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        for seq in sequences:
            s_str = str(seq.seq).upper()
            indices = [char_to_idx.get(c, -1) for c in s_str]
            seq_indices.append(indices)
            if len(indices) > max_len:
                max_len = len(indices)
        
        # Build buffer (padding with -1)
        self.all_indices = np.full((self.num_sequences, max_len), -1, dtype=np.int64)
        self.seq_lengths = np.zeros(self.num_sequences, dtype=np.int64)
        
        for i, indices in enumerate(seq_indices):
            l = len(indices)
            self.seq_lengths[i] = l
            self.all_indices[i, :l] = indices
            
        # Create thread pool and chunks
        if self.n_jobs > 1:
            self.executor = ThreadPoolExecutor(max_workers=self.n_jobs)
            print(f"Using {self.n_jobs} threads for batch processing")
            
            # Create chunks (start, end)
            chunk_size = int(np.ceil(self.num_sequences / self.n_jobs))
            self.chunks = []
            for i in range(0, self.num_sequences, chunk_size):
                self.chunks.append((i, min(i + chunk_size, self.num_sequences)))
            # print(f"Sequence chunks: {len(self.chunks)} x ~{chunk_size}") # Verbose
            
        else:
            self.executor = None
            print("Using sequential processing")
        
        # Pre-allocate for Dirichlet
        self._dirichlet_alpha = np.zeros(4, dtype=np.float64)
        
        # Initialize ensemble
        self.models = []
        self.model_likelihoods = []
        
        for _ in range(ensemble_size):
            model = self._sample_model()
            self.models.append(model)
            self.model_likelihoods.append(self._total_likelihood(model['motifs'], model['weights']))
    
    def _sample_model(self) -> Dict[str, Any]:
        """Sample a new model from prior."""
        motifs = []
        for _ in range(self.num_motifs):
            columns = sample_dirichlet_pwm(self.motif_length)
            motifs.append(CythonWeightMatrix(columns))
        weights = np.ones(self.num_motifs, dtype=np.float64)
        return {'motifs': motifs, 'weights': weights}
    
    
    def _total_likelihood(self, motifs: List[CythonWeightMatrix], weights: np.ndarray) -> float:
        """
        Compute total likelihood for a set of motifs.
        Uses Cython and Threading for speed.
        """
        
        # 1. Identify active motifs and build flattened arrays
        active_indices = [i for i, w in enumerate(weights) if w > 0.5]
        if not active_indices:
            return -2.0 * np.sum(self.seq_lengths)
            
        num_active = len(active_indices)
        
        # Parameters
        mean_len = np.mean(self.seq_lengths)
        trans = 1.0 / mean_len
        sum_trans = num_active * trans
        if sum_trans >= 1.0: sum_trans = 0.9999
        base_penalty = np.log2(1.0 - sum_trans)
        
        motif_pen = np.log2(trans)
        motif_penalties = np.full(num_active, motif_pen, dtype=np.float64)
        
        # Build flattened columns
        active_cols_list = [motifs[i].columns for i in active_indices] 
        all_motif_columns = np.concatenate(active_cols_list, axis=1)
        
        # Offsets and lengths
        lengths = np.array([cols.shape[1] for cols in active_cols_list], dtype=np.int64)
        offsets = np.zeros(num_active, dtype=np.int64)
        curr = 0
        for i in range(num_active):
            offsets[i] = curr
            curr += lengths[i]
            
        # 2. Parallel Execution
        if self.executor and self.num_sequences >= len(self.chunks):
            futures = []
            for start, end in self.chunks:
                # Slicing creates lightweight views passed to Cython
                # batch_likelihood releases GIL, so these run in parallel!
                futures.append(self.executor.submit(
                    batch_likelihood,
                    self.all_indices[start:end],
                    self.seq_lengths[start:end],
                    all_motif_columns,
                    offsets,
                    lengths,
                    motif_penalties,
                    base_penalty
                ))
            return sum(f.result() for f in futures)
        else:
            # Sequential call
            return batch_likelihood(
                self.all_indices,
                self.seq_lengths,
                all_motif_columns,
                offsets,
                lengths,
                motif_penalties,
                base_penalty
            )

    def _total_likelihood_from_columns(self, columns_list: List[np.ndarray], weights: np.ndarray) -> float:
        """Compute likelihood from raw column arrays."""
        temp_motifs = [CythonWeightMatrix(cols) for cols in columns_list]
        return self._total_likelihood(temp_motifs, weights)
    
    def _perturb_column_fast(self, log_probs: np.ndarray, alpha: float = 10.0) -> np.ndarray:
        """Fast column perturbation using Dirichlet proposal."""
        probs = np.power(2.0, log_probs)
        probs = np.maximum(probs, 1e-10)
        probs /= probs.sum()
        self._dirichlet_alpha[:] = probs * alpha + 0.1
        new_probs = np.random.dirichlet(self._dirichlet_alpha)
        return np.log2(new_probs + 1e-10)

    def _zap_column(self, cols: np.ndarray) -> np.ndarray:
        """Remove a random column (Zap)."""
        length = cols.shape[1]
        if length <= 5: return cols # Min length constraint
        idx = np.random.randint(length)
        return np.delete(cols, idx, axis=1)

    def _indel_column(self, cols: np.ndarray) -> np.ndarray:
        """Insert a random column (Indel)."""
        length = cols.shape[1]
        if length >= 20: return cols # Max length constraint
        idx = np.random.randint(length + 1)
        new_col = sample_dirichlet_pwm(1) # Shape (4, 1)
        return np.insert(cols, idx, new_col[:, 0], axis=1)
    
    def _decorrelate_adaptive(self, motifs, weights, min_hood, max_steps=100, early_stop=20) -> Tuple[List[Any], float]:
        """Adaptive M-H decorrelation with Variable Length moves."""
        # Convert to list of numpy arrays for mutation
        columns_list = [m.get_columns().copy() for m in motifs]
        
        # Initial score
        temp_motifs = [CythonWeightMatrix(cols) for cols in columns_list]
        current_hood = self._total_likelihood(temp_motifs, weights)
        
        no_improve_count = 0
        best_hood = current_hood
        
        for step in range(max_steps):
            # Choose move type: 0=Perturb(70%), 1=Zap(15%), 2=Indel(15%)
            move_type = np.random.choice([0, 1, 2], p=[0.7, 0.15, 0.15])
            m_idx = np.random.randint(self.num_motifs)
            
            old_cols = columns_list[m_idx].copy()
            
            if move_type == 0: # Perturb
                length = old_cols.shape[1]
                if length > 0:
                    c_idx = np.random.randint(length)
                    columns_list[m_idx][:, c_idx] = self._perturb_column_fast(old_cols[:, c_idx])
            elif move_type == 1: # Zap
                columns_list[m_idx] = self._zap_column(old_cols)
            else: # Indel
                columns_list[m_idx] = self._indel_column(old_cols)
            
            # Check acceptance
            new_hood = self._total_likelihood_from_columns(columns_list, weights)
            
            if new_hood > min_hood:
                current_hood = new_hood
                if new_hood > best_hood:
                    best_hood = new_hood
                    no_improve_count = 0
                else:
                    no_improve_count += 1
            else:
                # Revert
                columns_list[m_idx] = old_cols
                no_improve_count += 1
            
            if no_improve_count >= early_stop and step > 30:
                break
        
        new_motifs = [CythonWeightMatrix(cols.astype(np.float64)) for cols in columns_list]
        return new_motifs, current_hood
    
    def step(self) -> Tuple[Dict[str, Any], float]:
        """
        Perform one nested sampling step (Kill worst, Decorrelate survivor).
        
        Returns:
            Tuple[Dict, float]: The killed model and the likelihood threshold (worst).
        """
        min_idx = int(np.argmin(self.model_likelihoods))
        min_hood = self.model_likelihoods[min_idx]
        worst = self.models.pop(min_idx)
        self.model_likelihoods.pop(min_idx)
        
        survivor_idx = np.random.randint(len(self.models))
        survivor = self.models[survivor_idx]
        
        new_motifs, new_hood = self._decorrelate_adaptive(
            survivor['motifs'], survivor['weights'], min_hood
        )
        
        new_model = {'motifs': new_motifs, 'weights': survivor['weights'].copy()}
        self.models.append(new_model)
        self.model_likelihoods.append(new_hood)
        
        return worst, min_hood
    
    def get_best_model(self) -> Tuple[Dict[str, Any], float]:
        """Return the best model in ensemble."""
        best_idx = int(np.argmax(self.model_likelihoods))
        return self.models[best_idx], self.model_likelihoods[best_idx]
    
    def cleanup(self) -> None:
        """Shutdown thread pool."""
        if self.executor:
            self.executor.shutdown(wait=False)
