# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

"""
Cython-optimized core model with GIL-releasing memoryviews.
Manual parallelism friendly (no OpenMP dependency).
Includes support for MLX/GPU pre-calculated scores.
"""

import numpy as np
cimport numpy as np
from libc.math cimport log2, pow
from libc.stdlib cimport malloc, free

ctypedef np.float64_t DTYPE_t
ctypedef np.int64_t ITYPE_t

cdef double NEG_INF = -1e100


cdef inline double addlog2(double x, double y) nogil:
    """Fast log-space addition."""
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


cdef class CythonWeightMatrix:
    """Optimized weight matrix representation."""
    cdef public np.ndarray columns
    cdef double[:, :] _col_view
    cdef public int length
    
    def __init__(self, np.ndarray[DTYPE_t, ndim=2] columns):
        self.columns = columns.astype(np.float64, order='C')
        self._col_view = self.columns
        self.length = columns.shape[1]
    
    cpdef int get_length(self):
        return self.length
    
    cpdef np.ndarray get_columns(self):
        return self.columns
    
    cdef double[:, :] get_view(self) nogil:
        return self._col_view


cdef class CythonLikelihoodCalculator:
    """Legacy class kept for compatibility."""
    cdef np.ndarray indices
    cdef long[:] _indices_view
    cdef np.ndarray bg_scores
    cdef double[:] _bg_view
    cdef double bg_hood
    cdef double uncounted_exp
    cdef int seq_len
    cdef np.ndarray _matrix
    cdef double[:] _matrix_view
    
    def __init__(self, str seq_str, double uncounted_exp=1.0):
        pass
    
    def likelihood(self, list motifs, np.ndarray[DTYPE_t, ndim=1] weights):
        return 0.0

def create_calculator(seq_str, uncounted_exp=1.0):
    return CythonLikelihoodCalculator(seq_str, uncounted_exp)

def create_weight_matrix(columns):
    return CythonWeightMatrix(columns)


# -------------------------------------------------------------
# BATCH PROCESSING (Legacy CPU Scan)
# -------------------------------------------------------------

def batch_likelihood(
    long[:, :] all_indices,
    long[:] seq_lengths,
    double[:, :] all_motif_columns,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty
):
    """CPU-based scan + DP."""
    cdef int num_seqs = all_indices.shape[0]
    cdef int num_motifs = motif_offsets.shape[0]
    cdef double total_likelihood = 0.0
    cdef int i
    cdef double seq_hood
    
    with nogil:
        for i in range(num_seqs):
            seq_hood = _compute_single_seq(
                all_indices[i], 
                seq_lengths[i],
                all_motif_columns, 
                motif_offsets, 
                motif_lengths, 
                motif_penalties, 
                base_penalty,
                num_motifs
            )
            total_likelihood += seq_hood
            
    return total_likelihood


cdef double _compute_single_seq(
    long[:] indices, 
    int seq_len,
    double[:, :] all_motif_columns,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    int num_motifs
) nogil:
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL: return 0.0 
    cdef int idx_len = seq_len
    cdef double bg_score = -2.0
    cdef int i, m, wml
    cdef double score, emit, from_score, path_score
    matrix[0] = 0.0
    
    for i in range(1, idx_len + 1):
        score = matrix[i-1] + bg_score + base_penalty
        for m in range(num_motifs):
            wml = motif_lengths[m]
            if i >= wml:
                emit = _scan_at_pos(
                    indices, i - wml, wml, 
                    all_motif_columns, motif_offsets[m]
                )
                if emit > NEG_INF:
                    from_score = matrix[i - wml]
                    path_score = from_score + emit + motif_penalties[m]
                    score = addlog2(score, path_score)
        matrix[i] = score
    
    cdef double result = matrix[idx_len]
    free(matrix)
    return result


cdef double _scan_at_pos(long[:] indices, int start_pos, int length, double[:, :] all_cols, long col_offset) nogil:
    cdef double score = 0.0
    cdef int p
    cdef long idx
    for p in range(length):
        idx = indices[start_pos + p]
        if idx < 0: return NEG_INF
        score += all_cols[idx, col_offset + p]
    return score



