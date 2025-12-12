
import numpy as np
from nestedmica.model.core import Facette, LikelihoodCalculator
from nestedmica.model.motif import WeightMatrix

# Try to import Cython-optimized version, fallback to NumPy
try:
    from nestedmica.model.dp_likelihood import dp_likelihood_cython, scan_wm_cython
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False


class MotifFacette(Facette):
    def __init__(self, uncounted_expectation=1.0, rev_comp=False):
        self.uncounted_expectation = uncounted_expectation
        self.rev_comp = rev_comp
        self._calculator_cache = {}
        
    def get_likelihood_calculator(self, data):
        seq_id = id(data)
        if seq_id not in self._calculator_cache:
            self._calculator_cache[seq_id] = MotifUncountedLikelihood(self, data)
        return self._calculator_cache[seq_id]


class MotifUncountedLikelihood(LikelihoodCalculator):
    def __init__(self, facette, datum):
        self.facette = facette
        self.seq_str = str(datum.seq).upper()
        
        # Convert to indices
        char_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.indices = np.array([char_to_idx.get(c, -1) for c in self.seq_str], dtype=np.int64)
        
        # Background scores
        self.bg_score_per_base = -2.0
        self.bg_scores = np.full(len(self.indices), self.bg_score_per_base, dtype=np.float64)
        self.bg_hood = np.sum(self.bg_scores)
        
    def likelihood(self, contributions, weights):
        idx_len = len(self.indices)
        
        active_motifs = []
        active_transitions = []
        uncounted_exp = self.facette.uncounted_expectation
        
        for i, item in enumerate(contributions):
            weight = weights[i]
            if weight > 0.5:
                wwm = item.get_item()
                wm = wwm.get_weight_matrix()
                active_motifs.append(wm)
                trans = 1.0 * uncounted_exp / idx_len
                active_transitions.append(trans)
                
        num_motifs = len(active_motifs)
        if num_motifs == 0:
            return self.bg_hood
            
        # Pre-scan motifs
        motif_lengths = np.array([wm.get_length() for wm in active_motifs], dtype=np.int64)
        max_len = idx_len - min(motif_lengths) + 1 if len(motif_lengths) > 0 else 0
        
        motif_emission_scores = np.full((num_motifs, max_len), -1e100, dtype=np.float64)
        for m, wm in enumerate(active_motifs):
            if USE_CYTHON:
                columns = wm.get_columns().T.astype(np.float64, order='C')  # (motif_len, 4)
                scores = scan_wm_cython(self.indices, columns, wm.get_length())
            else:
                scores = self._scan_wm_numpy(wm)
            motif_emission_scores[m, :len(scores)] = scores
            
        # DP
        sum_trans = sum(active_transitions)
        if sum_trans >= 1.0:
            sum_trans = 0.9999
            
        base_penalty = np.log2(1.0 - sum_trans)
        motif_penalties = np.array([np.log2(t) for t in active_transitions], dtype=np.float64)
        
        if USE_CYTHON:
            result = dp_likelihood_cython(
                self.indices, self.bg_scores, base_penalty,
                motif_emission_scores, motif_lengths, motif_penalties, num_motifs
            )
        else:
            result = self._dp_numpy(idx_len, num_motifs, motif_emission_scores, 
                                    motif_lengths, motif_penalties, base_penalty)
        
        return result
    
    def _scan_wm_numpy(self, wm):
        """NumPy fallback for motif scanning."""
        seq_len = len(self.indices)
        motif_len = wm.get_length()
        if seq_len < motif_len:
            return np.array([], dtype=np.float64)
            
        from numpy.lib.stride_tricks import as_strided
        itemsize = self.indices.itemsize
        num_windows = seq_len - motif_len + 1
        windows = as_strided(self.indices, shape=(num_windows, motif_len), 
                            strides=(itemsize, itemsize))
        valid_mask = np.all(windows >= 0, axis=1)
        safe_windows = windows.copy()
        safe_windows[~valid_mask] = 0
        
        scores_matrix = wm.get_columns()[np.arange(motif_len), safe_windows]
        window_scores = np.sum(scores_matrix, axis=1)
        window_scores[~valid_mask] = -1e100
        
        return window_scores.astype(np.float64)
    
    def _dp_numpy(self, idx_len, num_motifs, motif_emission_scores, 
                  motif_lengths, motif_penalties, base_penalty):
        """NumPy fallback for DP."""
        matrix = np.zeros(idx_len + 1, dtype=np.float64)
        
        for i in range(1, idx_len + 1):
            score = matrix[i-1] + self.bg_scores[i-1] + base_penalty
            
            for m in range(num_motifs):
                wml = motif_lengths[m]
                if i >= wml:
                    emit = motif_emission_scores[m, i - wml]
                    if emit > -1e99:
                        from_score = matrix[i - wml]
                        path_score = from_score + emit + motif_penalties[m]
                        score = np.logaddexp2(score, path_score)
                        
            matrix[i] = score
        
        return matrix[idx_len]
