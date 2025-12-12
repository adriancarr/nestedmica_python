
import numpy as np
from nestedmica.model.core import Facette, LikelihoodCalculator
from nestedmica.maths import addLog2, log2, sumLog2

class LogisticSequenceFacette(Facette):
    def get_likelihood_calculator(self, data):
        # data is a Bio.SeqRecord or similar
        return LogisticSequenceLikelihoodCalculator(self, data)

class LogisticSequenceLikelihoodCalculator(LikelihoodCalculator):
    def __init__(self, facette, sequence_record):
        self.facette = facette
        self.sequence_record = sequence_record
        
        # Pre-process sequence to indices
        seq_str = str(sequence_record.seq).upper()
        self.indices = []
        for char in seq_str:
            if char == 'A': self.indices.append(0)
            elif char == 'C': self.indices.append(1)
            elif char == 'G': self.indices.append(2)
            elif char == 'T': self.indices.append(3)
            else: self.indices.append(-1)
        self.indices = np.array(self.indices, dtype=int)
            
        # Get label
        # BioPython annotations are a dictionary
        # In Mocca.java it checks "mocca.label" property
        self.label = 1 # Default to positive if not labeled? Or throw error.
        # Check annotations/features or description
        # Mocca.java parses description line "label=1"
        # We can simulate this or assume it's set in the record
        if "mocca.label" in sequence_record.annotations:
            self.label = int(sequence_record.annotations["mocca.label"])
        # Fallback parsing description if needed, implemented in seq.py hopefully

    def likelihood(self, contributions, weights):
        # contributions is a list of ContributionItems
        # weights is the mixture vector (numpy array)
        
        eta = 0.0
        
        # In Logistic model, we usually treat the mixture as 1-1 mapping or just use contributions
        # The Java code iterates contributions
        
        for i, ci in enumerate(contributions):
            wwm = ci.get_item() # WeightedWeightMatrix
            wm = wwm.get_weight_matrix()
            # wm.columns is (L, 4)
            
            # Simple scoring without "Views" complexity
            # score = sum of scores across all positions
            
            scores = self._scan_wm(wm)
            
            wml = wm.get_length()
            odds = wml * log2(0.25)
            
            # x = log_sum_exp(scores)
            x = sumLog2(scores) # Using our sumLog2 equivalent
            
            eta += wwm.get_weight() * (x - odds)
            
        pi = 1.0 / (1.0 + np.exp(-eta))
        
        if self.label >= 0:
            return log2(pi)
        else:
            return log2(1.0 - pi)

    def _scan_wm(self, wm):
        # wm.columns is (L, 4)
        # Sequence indices is (N,)
        
        # Vectorized scanning
        # Use numpy striding
        
        seq_len = len(self.indices)
        motif_len = wm.get_length()
        
        if seq_len < motif_len:
            return np.array([-np.inf]) # Or similar
            
        # Create a view of the sequence with windows
        # shape (N - L + 1, L)
        from numpy.lib.stride_tricks import as_strided
        
        # Be careful with stride tricks.
        # indices is 1D array of int
        
        itemsize = self.indices.itemsize
        num_windows = seq_len - motif_len + 1
        window_shape = (num_windows, motif_len)
        strides = (itemsize, itemsize)
        
        windows = as_strided(self.indices, shape=window_shape, strides=strides)
        
        # Check for -1 (invalid)
        # Any window containing -1 is invalid
        valid_mask = np.all(windows >= 0, axis=1)
        
        # wm.columns is (L, 4). Transpose to (4, L) or keep
        # We want to lookup scores.
        # wm.columns[pos, base]
        
        # windows[i, j] is the base at pos j in window i
        # We want sum(wm.columns[j, windows[i, j]]) for j in 0..L-1
        
        # Advanced indexing:
        # wm_cols = wm.columns # (L, 4)
        # values = wm_cols[np.arange(motif_len), windows] # Shape (N, L) ? No
        # values[i, j] should be wm_cols[j, windows[i, j]]
        
        # Correct indexing:
        # windows is (NumWindows, L)
        # col_indices is (1, L) -> (NumWindows, L)
        col_indices = np.arange(motif_len)
        
        # Replace -1 with 0 temporarily to avoid index error, but mask later
        safe_windows = windows.copy()
        safe_windows[~valid_mask] = 0
        
        # Lookup
        scores_matrix = wm.get_columns()[col_indices, safe_windows] # This might broadcast wrong?
        # wm.columns is (L, 4)
        # We want: for each window w, for each pos p, val = wm.columns[p, w[p]]
        
        # Let's verify indexing.
        # wm.columns[p, c] -> value
        # We use integer array indexing.
        # A[I, J]
        # I = col_indices (broadcast to NumWindows, L)
        # J = safe_windows (NumWindows, L)
        
        scores_matrix = wm.get_columns()[col_indices, safe_windows] 
        # Result shape (NumWindows, L)?
        # Actually indexing with two arrays of same shape returns array of that shape.
        # So yes, (NumWindows, L).
        
        window_scores = np.sum(scores_matrix, axis=1)
        
        # Apply mask
        window_scores[~valid_mask] = -np.inf
        
        return window_scores
