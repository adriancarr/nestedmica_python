#!/usr/bin/env python3
"""
mocca_fast.py - Cython-optimized motif finder with Data-Oriented Threaded Parallelism.
Uses ThreadPoolExecutor + GIL-releasing Cython batch function. (No OpenMP dependency).
"""

import argparse
import numpy as np
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

# Import Cython-optimized modules
from nestedmica.model.cython_model import CythonWeightMatrix, batch_likelihood


def sample_dirichlet_pwm(length, alpha=1.0):
    """Sample a random PWM from Dirichlet prior."""
    columns = np.random.dirichlet([alpha] * 4, size=length).T
    return np.log2(columns + 1e-10).astype(np.float64)


class FastTrainer:
    """Optimized trainer with Threaded Batch Parallelism."""
    
    def __init__(self, sequences, num_motifs, motif_length, ensemble_size, n_jobs=-1):
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
            print(f"Sequence chunks: {len(self.chunks)} x ~{chunk_size}")
            
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
    
    def _sample_model(self):
        """Sample a new model from prior."""
        motifs = []
        for _ in range(self.num_motifs):
            columns = sample_dirichlet_pwm(self.motif_length)
            motifs.append(CythonWeightMatrix(columns))
        weights = np.ones(self.num_motifs, dtype=np.float64)
        return {'motifs': motifs, 'weights': weights}
    
    
    def _total_likelihood(self, motifs, weights):
        """Compute total likelihood."""
        
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

    def _total_likelihood_from_columns(self, columns_list, weights):
        """Compute likelihood from raw column arrays."""
        temp_motifs = [CythonWeightMatrix(cols) for cols in columns_list]
        return self._total_likelihood(temp_motifs, weights)
    
    def _perturb_column_fast(self, log_probs, alpha=10.0):
        """Fast column perturbation."""
        probs = np.power(2.0, log_probs)
        probs = np.maximum(probs, 1e-10)
        probs /= probs.sum()
        self._dirichlet_alpha[:] = probs * alpha + 0.1
        new_probs = np.random.dirichlet(self._dirichlet_alpha)
        return np.log2(new_probs + 1e-10)

    def _zap_column(self, cols):
        """Remove a random column (Zap)."""
        length = cols.shape[1]
        if length <= 5: return cols # Min length constraint
        idx = np.random.randint(length)
        return np.delete(cols, idx, axis=1)

    def _indel_column(self, cols):
        """Insert a random column (Indel)."""
        length = cols.shape[1]
        if length >= 20: return cols # Max length constraint
        idx = np.random.randint(length + 1)
        new_col = sample_dirichlet_pwm(1) # Shape (4, 1)
        # Convert to log space if sample_dirichlet_pwm returns log space (it does)
        # But wait, sample_dirichlet_pwm returns shape (4, L). For L=1 it is (4, 1).
        # CythonWeightMatrix expects (4, L).
        # We handle raw numpy arrays here.
        return np.insert(cols, idx, new_col[:, 0], axis=1)
    
    def _decorrelate_adaptive(self, motifs, weights, min_hood, max_steps=100, early_stop=20):
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
    
    def step(self):
        """Perform one nested sampling step."""
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
    
    def get_best_model(self):
        """Return the best model in ensemble."""
        best_idx = int(np.argmax(self.model_likelihoods))
        return self.models[best_idx], self.model_likelihoods[best_idx]
    
    def cleanup(self):
        """Shutdown thread pool."""
        if self.executor:
            self.executor.shutdown(wait=False)


import pickle
import os

class AsciiPlotter:
    """Simple ASCII graph plotter for monitoring."""
    def __init__(self, height=10, width=50):
        self.height = height
        self.width = width
        self.history = []
        
    def add(self, value):
        self.history.append(value)
        if len(self.history) > self.width:
            self.history.pop(0)
            
    def plot(self):
        if not self.history: return ""
        min_v = min(self.history)
        max_v = max(self.history)
        rnge = max_v - min_v if max_v != min_v else 1.0
        
        # Normalize to 0..height-1
        rows = [[' ' for _ in range(len(self.history))] for _ in range(self.height)]
        for x, val in enumerate(self.history):
            normalized = int((val - min_v) / rnge * (self.height - 1))
            rows[self.height - 1 - normalized][x] = '*'
            
        # Draw
        lines = []
        lines.append(f"    {max_v:.2f}")
        for r in rows:
            lines.append("    |" + "".join(r) + "|")
        lines.append(f"    {min_v:.2f}")
        return "\n".join(lines)

def save_checkpoint(trainer, filename, cycle):
    """Save trainer state to file."""
    # Extract raw data from Cython objects
    serialized_models = []
    for model in trainer.models:
        motifs_data = [m.get_columns().copy() for m in model['motifs']]
        serialized_models.append({'motifs': motifs_data, 'weights': model['weights'].copy()})
        
    state = {
        'cycle': cycle,
        'models': serialized_models,
        'model_likelihoods': trainer.model_likelihoods,
        'num_motifs': trainer.num_motifs,
        'motif_length': trainer.motif_length,
        'ensemble_size': trainer.ensemble_size
    }
    
    with open(filename, 'wb') as f:
        pickle.dump(state, f)
    print(f"  [Checkpoint saved to {filename}]")

def load_checkpoint(filename, sequences, n_jobs):
    """Load trainer from checkpoint."""
    print(f"Loading checkpoint from {filename}...")
    with open(filename, 'rb') as f:
        state = pickle.load(f)
        
    trainer = FastTrainer(sequences, state['num_motifs'], state['motif_length'], 
                         state['ensemble_size'], n_jobs=n_jobs)
    
    # Restore models
    trainer.models = []
    for s_model in state['models']:
        motifs = [CythonWeightMatrix(cols) for cols in s_model['motifs']]
        trainer.models.append({'motifs': motifs, 'weights': s_model['weights']})
        
    trainer.model_likelihoods = state['model_likelihoods']
    return trainer, state['cycle']

def main():
    parser = argparse.ArgumentParser(description='Fast Cython-optimized Mocca')
    parser.add_argument('-seqs', required=True, help='FASTA file')
    parser.add_argument('-numMotifs', type=int, default=1)
    parser.add_argument('-motifLength', type=int, default=10)
    parser.add_argument('-maxCycles', type=int, default=500)
    parser.add_argument('-ensembleSize', type=int, default=20)
    parser.add_argument('-threads', type=int, default=-1, help='Number of threads (-1=all)')
    parser.add_argument('-out', required=True, help='Output XMS file')
    parser.add_argument('-stopIqr', type=float, default=0.01)
    parser.add_argument('-checkpoint', help='File to save checkpoints to')
    parser.add_argument('-checkpointInterval', type=int, default=100)
    parser.add_argument('-restart', help='Checkpoint file to restart from')
    args = parser.parse_args()
    
    print(f"Loading sequences from {args.seqs}...")
    sequences = list(SeqIO.parse(args.seqs, "fasta"))
    print(f"Loaded {len(sequences)} sequences.")
    
    start_cycle = 0
    if args.restart:
        try:
            trainer, start_cycle = load_checkpoint(args.restart, sequences, args.threads)
            print(f"Resuming from Cycle {start_cycle}")
        except Exception as e:
            print(f"Failed to load checkpoint: {e}")
            return
    else:
        print("Initializing Cython Data-Oriented Threaded trainer...")
        trainer = FastTrainer(sequences, args.numMotifs, args.motifLength, 
                             args.ensembleSize, n_jobs=args.threads)
    
    print("Starting sampling...")
    import time
    start_time = time.time()
    
    plotter_l = AsciiPlotter(height=5, width=60)
    plotter_iqr = AsciiPlotter(height=5, width=60)
    
    for cycle in range(start_cycle, args.maxCycles):
        worst, hood = trainer.step()
        if cycle % 100 == 0:
            best_model, best_hood = trainer.get_best_model()

            # Calculate IQR
            likelihoods = np.array(trainer.model_likelihoods)
            q75, q25 = np.percentile(likelihoods, [75, 25])
            iqr = q75 - q25

            # Plotting
            plotter_l.add(best_hood)
            plotter_iqr.add(iqr)
            
            print(f"\nCycle {cycle}: L={best_hood:.2f} IQR={iqr:.4f}")
            print("Likelihood History:")
            print(plotter_l.plot())
            # print("IQR History:")
            # print(plotter_iqr.plot()) # Optional
            
            # Checkpoint
            if args.checkpoint and cycle % args.checkpointInterval == 0 and cycle > start_cycle:
                save_checkpoint(trainer, args.checkpoint, cycle)
            
            if iqr < args.stopIqr and cycle > start_cycle + 500:
                print(f"Converged at Cycle {cycle} (IQR {iqr:.4f} < {args.stopIqr})")
                break
    
    end_time = time.time()
    print(f"Finished in {end_time - start_time:.2f}s")
    
    # Save final checkpoint if requested
    if args.checkpoint:
         save_checkpoint(trainer, args.checkpoint, args.maxCycles)
    
    best_model, best_hood = trainer.get_best_model()
    print(f"Best Likelihood: {best_hood}")
    
    trainer.cleanup()
    
    with open(args.out, 'w') as f:
        f.write("<motifs>\n")
        f.write(f"<!-- Final Cycle: {args.maxCycles}, Best Likelihood: {best_hood} -->\n")
        for i, wm in enumerate(best_model['motifs']):
            f.write(f"  <motif id='{i}' weight='1.0'>\n")
            f.write(f"    <weightMatrix length='{wm.get_length()}'>\n")
            columns = wm.get_columns()
            for pos in range(wm.get_length()):
                probs = np.power(2.0, columns[:, pos])
                probs = probs / probs.sum()
                f.write(f"      {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}\n")
            f.write("    </weightMatrix>\n")
            f.write("  </motif>\n")
        f.write("</motifs>\n")
    
    print(f"Saved motifs to {args.out}")


if __name__ == "__main__":
    main()
