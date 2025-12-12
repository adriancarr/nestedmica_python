"""
Core Nested MICA Trainer implementation (Cython-Accelerated).
"""

import numpy as np
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from typing import List, Tuple, Dict, Any, Optional

# Import Cython-optimized modules
from nestedmica.model.cython_model import (
    CythonWeightMatrix, 
    batch_likelihood, 
    batch_likelihood_dual,
    reverse_complement_columns
)
from nestedmica.model.background import learn_markov_background, print_background_stats

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
        both_strands (bool): If True, scan both forward and reverse complement.
    """
    
    def __init__(self, sequences: List[Any], num_motifs: int, motif_length: int, 
                 ensemble_size: int, n_jobs: int = -1, both_strands: bool = True,
                 bg_order: int = 3, adaptive_mcmc: bool = False):
        """
        Initialize the trainer.

        Args:
            sequences (List): BioPython sequence objects.
            num_motifs (int): Number of motifs.
            motif_length (int): Initial target length.
            ensemble_size (int): Population size.
            n_jobs (int): Threads (-1 for all cores).
            both_strands (bool): Scan both strands (default: True).
            bg_order (int): Background Markov order (default: 3).
            adaptive_mcmc (bool): Use adaptive proposal distribution (default: True).
        """
        self.sequences = sequences
        self.num_motifs = num_motifs
        self.motif_length = motif_length
        self.ensemble_size = ensemble_size
        self.num_sequences = len(sequences)
        self.both_strands = both_strands
        self.bg_order = bg_order
        self.adaptive_mcmc = adaptive_mcmc
        
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
        
        # Learn higher-order background model
        self.bg_model = learn_markov_background(sequences, order=bg_order)
        print_background_stats(self.bg_model, bg_order)
        
        # Adaptive MCMC: per-motif proposal concentration and acceptance tracking
        self.proposal_alphas = np.full(num_motifs, 10.0, dtype=np.float64)
        self.accept_counts = np.zeros(num_motifs, dtype=np.float64)
        self.attempt_counts = np.zeros(num_motifs, dtype=np.float64)
        if adaptive_mcmc:
            print("Adaptive MCMC proposals: ENABLED")
        
        # Thread-safe random number generator
        self._rng = np.random.default_rng()
        
        # Pre-allocate for Dirichlet
        self._dirichlet_alpha = np.zeros(4, dtype=np.float64)
        
        # Evidence tracking
        self.log_evidence = -np.inf # Log2 scale
        self.step_count = 0
        
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
        Supports both-strand scanning if self.both_strands is True.
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
        
        # Build flattened columns (forward strand)
        active_cols_list = [motifs[i].columns for i in active_indices] 
        all_motif_columns_fwd = np.concatenate(active_cols_list, axis=1)
        
        # Offsets and lengths
        lengths = np.array([cols.shape[1] for cols in active_cols_list], dtype=np.int64)
        offsets = np.zeros(num_active, dtype=np.int64)
        curr = 0
        for i in range(num_active):
            offsets[i] = curr
            curr += lengths[i]
        
        # 2. Compute reverse complement columns if both_strands enabled
        if self.both_strands:
            rc_cols_list = [reverse_complement_columns(cols) for cols in active_cols_list]
            all_motif_columns_rc = np.concatenate(rc_cols_list, axis=1)
        
        # 3. Execute likelihood calculation
        if self.both_strands:
            # Dual-strand mode
            if self.executor and self.num_sequences >= len(self.chunks):
                futures = []
                for start, end in self.chunks:
                    futures.append(self.executor.submit(
                        batch_likelihood_dual,
                        self.all_indices[start:end],
                        self.seq_lengths[start:end],
                        all_motif_columns_fwd,
                        all_motif_columns_rc,
                        offsets,
                        lengths,
                        motif_penalties,
                        base_penalty,
                        self.bg_model,
                        self.bg_order
                    ))
                return sum(f.result() for f in futures)
            else:
                return batch_likelihood_dual(
                    self.all_indices,
                    self.seq_lengths,
                    all_motif_columns_fwd,
                    all_motif_columns_rc,
                    offsets,
                    lengths,
                    motif_penalties,
                    base_penalty,
                    self.bg_model,
                    self.bg_order
                )
        else:
            # Forward-strand only mode
            if self.executor and self.num_sequences >= len(self.chunks):
                futures = []
                for start, end in self.chunks:
                    futures.append(self.executor.submit(
                        batch_likelihood,
                        self.all_indices[start:end],
                        self.seq_lengths[start:end],
                        all_motif_columns_fwd,
                        offsets,
                        lengths,
                        motif_penalties,
                        base_penalty,
                        self.bg_model,
                        self.bg_order
                    ))
                return sum(f.result() for f in futures)
            else:
                return batch_likelihood(
                    self.all_indices,
                    self.seq_lengths,
                    all_motif_columns_fwd,
                    offsets,
                    lengths,
                    motif_penalties,
                    base_penalty,
                    self.bg_model,
                    self.bg_order
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
        new_probs = self._rng.dirichlet(self._dirichlet_alpha)
        return np.log2(new_probs + 1e-10)

    def _delete_column(self, cols: np.ndarray) -> np.ndarray:
        """Remove a random column (Delete)."""
        length = cols.shape[1]
        if length <= 5: return cols # Min length constraint
        idx = self._rng.integers(length)
        return np.delete(cols, idx, axis=1)

    def _insert_column(self, cols: np.ndarray) -> np.ndarray:
        """Insert a random column (Insert)."""
        length = cols.shape[1]
        if length >= 20: return cols # Max length constraint
        idx = self._rng.integers(length + 1)
        new_col = sample_dirichlet_pwm(1) # Shape (4, 1)
        return np.insert(cols, idx, new_col[:, 0], axis=1)
    
    def _decorrelate_adaptive(self, motifs, weights, min_hood, max_steps=100, early_stop=20) -> Tuple[List[Any], float]:
        """Adaptive M-H decorrelation with Variable Length moves and adaptive proposals."""
        # Convert to list of numpy arrays for mutation
        columns_list = [m.get_columns().copy() for m in motifs]
        
        # Initial score
        temp_motifs = [CythonWeightMatrix(cols) for cols in columns_list]
        current_hood = self._total_likelihood(temp_motifs, weights)
        
        no_improve_count = 0
        best_hood = current_hood
        
        # Per-motif tracking for this decorrelation cycle
        local_accepts = np.zeros(self.num_motifs, dtype=np.float64)
        local_attempts = np.zeros(self.num_motifs, dtype=np.float64)
        
        for step in range(max_steps):
            # Choose move type: 0=Perturb(70%), 1=Delete(15%), 2=Insert(15%)
            move_type = self._rng.choice([0, 1, 2], p=[0.7, 0.15, 0.15])
            m_idx = self._rng.integers(self.num_motifs)
            
            old_cols = columns_list[m_idx].copy()
            
            if move_type == 0: # Perturb
                length = old_cols.shape[1]
                if length > 0:
                    c_idx = self._rng.integers(length)
                    # Use per-motif adaptive alpha
                    alpha = self.proposal_alphas[m_idx] if self.adaptive_mcmc else 10.0
                    columns_list[m_idx][:, c_idx] = self._perturb_column_fast(old_cols[:, c_idx], alpha=alpha)
                    local_attempts[m_idx] += 1
            elif move_type == 1: # Delete
                columns_list[m_idx] = self._delete_column(old_cols)
            else: # Insert
                columns_list[m_idx] = self._insert_column(old_cols)            
            
            # Check acceptance
            new_hood = self._total_likelihood_from_columns(columns_list, weights)
            
            if new_hood > min_hood:
                current_hood = new_hood
                if move_type == 0:
                    local_accepts[m_idx] += 1
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
        
        # Adapt proposal alphas based on acceptance rate
        if self.adaptive_mcmc:
            for m in range(self.num_motifs):
                if local_attempts[m] >= 5:
                    rate = local_accepts[m] / local_attempts[m]
                    # Target ~23% acceptance rate
                    if rate > 0.35:
                        self.proposal_alphas[m] *= 1.3  # Too many accepts: explore more
                    elif rate < 0.15:
                        self.proposal_alphas[m] /= 1.3  # Too many rejects: fine-tune
                    # Clamp to [1.0, 1000.0]
                    self.proposal_alphas[m] = np.clip(self.proposal_alphas[m], 1.0, 1000.0)
            
            # Update global tracking
            self.accept_counts += local_accepts
            self.attempt_counts += local_attempts
        
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
        
        survivor_idx = self._rng.integers(len(self.models))
        survivor = self.models[survivor_idx]
        
        new_motifs, new_hood = self._decorrelate_adaptive(
            survivor['motifs'], survivor['weights'], min_hood
        )
        
        new_model = {'motifs': new_motifs, 'weights': survivor['weights'].copy()}
        self.models.append(new_model)
        self.model_likelihoods.append(new_hood)
        
        # --- Evidence Calculation (Log Space) ---
        # Prior volume shrinkage factor: X_i = exp(-i / N)
        # Weight w_i = 0.5 * (X_{i-1} - X_{i+1})
        # log(w_i) = log(0.5) + log(X_{i-1} - X_{i+1})
        
        # For numerical stability with X close to 0:
        # X_i = exp(-i/N)
        # X_{i-1} = exp(-(i-1)/N) = X_i * exp(1/N)
        # X_{i+1} = exp(-(i+1)/N) = X_i * exp(-1/N)
        # w_i = 0.5 * X_i * (exp(1/N) - exp(-1/N))
        # log(w_i) = log(0.5) + log(X_i) + log(exp(1/N) - exp(-1/N))
        
        # We track cumulative LogZ.
        # Delta Z = L_min * w_i
        # log(Delta Z) = log(L_min) + log(w_i)
        
        # Count is actually 'iterations done', which we need to track if we restart.
        # But for new runs: i = 1, 2, ...
        # We need a counter.
        self.step_count += 1
        i = self.step_count
        N = self.ensemble_size
        
        # Calculate log(w_i)
        # log(X_i) = -i / N
        log_Xi = -float(i) / N
        term = np.log(np.exp(1.0/N) - np.exp(-1.0/N))
        log_wi = np.log(0.5) + log_Xi + term
        
        # Update LogZ
        # min_hood is Log Likelihood (ln L or log2 L? Code uses log2 everywhere)
        # Wait, Evidence integral requires natural log usually, or consistent base.
        # Code uses log2 for PWMs and likelihoods.
        # Let's stick to log2 for Z as well to differ only by scale factor.
        # Z_2 = Integral (2^L * w)
        # log2(dZ) = L_min + log2(w_i)
        
        log2_wi = log_wi / np.log(2) # Convert natural log weight to log2
        log2_dZ = min_hood + log2_wi
        
        if self.log_evidence == -np.inf:
            self.log_evidence = log2_dZ
        else:
            # logaddexp2 is not standard numpy, use logaddexp with conversion
            # log2(a + b) = log2(2^a + 2^b) = log2(e^(a ln2) + e^(b ln2)) 
            #             = np.logaddexp(a*ln2, b*ln2) / ln2
            ln2 = np.log(2)
            self.log_evidence = np.logaddexp(self.log_evidence * ln2, log2_dZ * ln2) / ln2

        return worst, min_hood
    
    def get_best_model(self) -> Tuple[Dict[str, Any], float]:
        """Return the best model in ensemble."""
        best_idx = int(np.argmax(self.model_likelihoods))
        return self.models[best_idx], self.model_likelihoods[best_idx]
    
    def cleanup(self) -> None:
        """Shutdown thread pool."""
        if self.executor:
            self.executor.shutdown(wait=False)
