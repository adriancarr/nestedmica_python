# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
"""
Cython-optimized DP likelihood calculation for motif finding.
Compile with: cythonize -i dp_likelihood.pyx
"""

import numpy as np
cimport numpy as np
from libc.math cimport log2, pow

# Type definitions
ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t ITYPE_t

cdef double NEG_INF = -1e100


cdef inline double addlog2(double x, double y) nogil:
    """Fast log-space addition: log2(2^x + 2^y)"""
    cdef double diff
    if x > y:
        diff = y - x
        if diff > -50:
            return x + log2(1.0 + pow(2.0, diff))
        return x
    else:
        diff = x - y
        if diff > -50:
            return y + log2(1.0 + pow(2.0, diff))
        return y


def dp_likelihood_cython(
    np.ndarray[ITYPE_t, ndim=1] indices,
    np.ndarray[DTYPE_t, ndim=1] bg_scores,
    double base_penalty,
    np.ndarray[DTYPE_t, ndim=2] motif_emission_scores,
    np.ndarray[ITYPE_t, ndim=1] motif_lengths,
    np.ndarray[DTYPE_t, ndim=1] motif_penalties,
    int num_motifs
):
    """
    Cython-optimized DP likelihood calculation.
    
    Parameters:
    -----------
    indices : int64 array of sequence indices (A=0, C=1, G=2, T=3, N=-1)
    bg_scores : float64 array of background log2 scores per position
    base_penalty : float64 log2 penalty for background transition
    motif_emission_scores : 2D float64 array (num_motifs, max_positions)
    motif_lengths : int64 array of motif lengths
    motif_penalties : float64 array of log2 motif transition penalties
    num_motifs : int number of active motifs
    
    Returns:
    --------
    float64 : log2 likelihood of sequence
    """
    cdef int idx_len = len(indices)
    cdef np.ndarray[DTYPE_t, ndim=1] matrix = np.zeros(idx_len + 1, dtype=np.float64)
    
    cdef int i, m, wml
    cdef double score, emit, from_score, path_score
    
    matrix[0] = 0.0
    
    for i in range(1, idx_len + 1):
        # Background transition
        score = matrix[i-1] + bg_scores[i-1] + base_penalty
        
        # Motif transitions
        for m in range(num_motifs):
            wml = motif_lengths[m]
            if i >= wml:
                emit = motif_emission_scores[m, i - wml]
                if emit > NEG_INF:
                    from_score = matrix[i - wml]
                    path_score = from_score + emit + motif_penalties[m]
                    score = addlog2(score, path_score)
                    
        matrix[i] = score
    
    return matrix[idx_len]


def scan_wm_cython(
    np.ndarray[ITYPE_t, ndim=1] indices,
    np.ndarray[DTYPE_t, ndim=2] columns,
    int motif_len
):
    """
    Cython-optimized motif scanning.
    
    Parameters:
    -----------
    indices : int64 array of sequence indices
    columns : 2D float64 array of motif weights (motif_len, 4)
    motif_len : int motif length
    
    Returns:
    --------
    float64 array : log2 scores at each valid position
    """
    cdef int seq_len = len(indices)
    if seq_len < motif_len:
        return np.empty(0, dtype=np.float64)
    
    cdef int num_windows = seq_len - motif_len + 1
    cdef np.ndarray[DTYPE_t, ndim=1] window_scores = np.empty(num_windows, dtype=np.float64)
    
    cdef int w, p, idx
    cdef double score
    cdef bint valid
    
    for w in range(num_windows):
        score = 0.0
        valid = True
        for p in range(motif_len):
            idx = indices[w + p]
            if idx < 0:
                valid = False
                break
            score = score + columns[p, idx]
        
        if valid:
            window_scores[w] = score
        else:
            window_scores[w] = NEG_INF
    
    return window_scores
