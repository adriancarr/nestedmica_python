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


def reverse_complement_columns(np.ndarray[DTYPE_t, ndim=2] columns):
    """
    Compute reverse complement of a log2-space PWM.
    - Reverse the order of columns (right-to-left).
    - Swap A<->T (rows 0<->3) and C<->G (rows 1<->2).
    
    Args:
        columns: 4 x L array of log2 probabilities (A, C, G, T order).
        
    Returns:
        4 x L array representing the reverse complement motif.
    """
    # Reverse column order
    cdef np.ndarray[DTYPE_t, ndim=2] rc = columns[:, ::-1].copy()
    # Swap bases: A(0)<->T(3), C(1)<->G(2)
    cdef np.ndarray[DTYPE_t, ndim=2] swapped = np.empty_like(rc)
    swapped[0, :] = rc[3, :]  # A <- T
    swapped[1, :] = rc[2, :]  # C <- G
    swapped[2, :] = rc[1, :]  # G <- C
    swapped[3, :] = rc[0, :]  # T <- A
    return swapped


# -------------------------------------------------------------
# BATCH PROCESSING (Legacy CPU Scan)
# -------------------------------------------------------------

cdef inline long get_context_idx(long[:] indices, int pos, int order) nogil:
    """
    Compute context index from previous 'order' bases.
    Uses base-4 encoding: A=0, C=1, G=2, T=3
    Returns -1 if any base is invalid (N).
    """
    cdef long idx = 0
    cdef int i
    cdef long base
    for i in range(order):
        base = indices[pos - order + i]
        if base < 0:
            return -1
        idx = idx * 4 + base
    return idx


def batch_likelihood(
    long[:, :] all_indices,
    long[:] seq_lengths,
    double[:, :] all_motif_columns,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    double[:, :] bg_model,
    int bg_order
):
    """CPU-based scan + DP (forward strand only) with higher-order background."""
    cdef int num_seqs = all_indices.shape[0]
    cdef int num_motifs = motif_offsets.shape[0]
    cdef double total_likelihood = 0.0
    cdef int i
    cdef double seq_hood
    
    with nogil:
        for i in range(num_seqs):
            seq_hood = _compute_single_seq_bg(
                all_indices[i], 
                seq_lengths[i],
                all_motif_columns, 
                motif_offsets, 
                motif_lengths, 
                motif_penalties, 
                base_penalty,
                num_motifs,
                bg_model,
                bg_order
            )
            total_likelihood += seq_hood
            
    return total_likelihood


def batch_likelihood_dual(
    long[:, :] all_indices,
    long[:] seq_lengths,
    double[:, :] all_motif_columns_fwd,
    double[:, :] all_motif_columns_rc,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    double[:, :] bg_model,
    int bg_order
):
    """CPU-based scan + DP with BOTH strands and higher-order background."""
    cdef int num_seqs = all_indices.shape[0]
    cdef int num_motifs = motif_offsets.shape[0]
    cdef double total_likelihood = 0.0
    cdef int i
    cdef double seq_hood
    
    with nogil:
        for i in range(num_seqs):
            seq_hood = _compute_single_seq_dual_bg(
                all_indices[i], 
                seq_lengths[i],
                all_motif_columns_fwd,
                all_motif_columns_rc,
                motif_offsets, 
                motif_lengths, 
                motif_penalties, 
                base_penalty,
                num_motifs,
                bg_model,
                bg_order
            )
            total_likelihood += seq_hood
            
    return total_likelihood


cdef double _compute_single_seq_bg(
    long[:] indices, 
    int seq_len,
    double[:, :] all_motif_columns,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    int num_motifs,
    double[:, :] bg_model,
    int bg_order
) nogil:
    """DP with higher-order background model (forward strand only)."""
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL: return 0.0 
    cdef int idx_len = seq_len
    cdef int i, m, wml
    cdef double score, emit, from_score, path_score, bg_score
    cdef long ctx_idx, base_idx
    matrix[0] = 0.0
    
    for i in range(1, idx_len + 1):
        # Get background score for this position
        base_idx = indices[i-1]
        if base_idx < 0:
            bg_score = -2.0  # Fallback for N
        elif bg_order == 0:
            bg_score = bg_model[0, base_idx]
        elif i > bg_order:
            ctx_idx = get_context_idx(indices, i-1, bg_order)
            if ctx_idx < 0:
                bg_score = -2.0  # Fallback for N in context
            else:
                bg_score = bg_model[ctx_idx, base_idx]
        else:
            bg_score = -2.0  # Not enough context yet
        
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


# Keep old function for backward compatibility (uses fixed bg_score)
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


cdef double _compute_single_seq_dual(
    long[:] indices, 
    int seq_len,
    double[:, :] all_motif_columns_fwd,
    double[:, :] all_motif_columns_rc,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    int num_motifs
) nogil:
    """DP with both strands (legacy, fixed bg_score)."""
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL: return 0.0 
    cdef int idx_len = seq_len
    cdef double bg_score = -2.0
    cdef int i, m, wml
    cdef double score, emit_fwd, emit_rc, emit, from_score, path_score
    matrix[0] = 0.0
    
    for i in range(1, idx_len + 1):
        score = matrix[i-1] + bg_score + base_penalty
        for m in range(num_motifs):
            wml = motif_lengths[m]
            if i >= wml:
                emit_fwd = _scan_at_pos(
                    indices, i - wml, wml, 
                    all_motif_columns_fwd, motif_offsets[m]
                )
                emit_rc = _scan_at_pos(
                    indices, i - wml, wml, 
                    all_motif_columns_rc, motif_offsets[m]
                )
                if emit_fwd > emit_rc:
                    emit = emit_fwd
                else:
                    emit = emit_rc
                    
                if emit > NEG_INF:
                    from_score = matrix[i - wml]
                    path_score = from_score + emit + motif_penalties[m]
                    score = addlog2(score, path_score)
        matrix[i] = score
    
    cdef double result = matrix[idx_len]
    free(matrix)
    return result


cdef double _compute_single_seq_dual_bg(
    long[:] indices, 
    int seq_len,
    double[:, :] all_motif_columns_fwd,
    double[:, :] all_motif_columns_rc,
    long[:] motif_offsets,
    long[:] motif_lengths,
    double[:] motif_penalties,
    double base_penalty,
    int num_motifs,
    double[:, :] bg_model,
    int bg_order
) nogil:
    """DP with both strands and higher-order background model."""
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL: return 0.0 
    cdef int idx_len = seq_len
    cdef int i, m, wml
    cdef double score, emit_fwd, emit_rc, emit, from_score, path_score, bg_score
    cdef long ctx_idx, base_idx
    matrix[0] = 0.0
    
    for i in range(1, idx_len + 1):
        # Get background score for this position
        base_idx = indices[i-1]
        if base_idx < 0:
            bg_score = -2.0
        elif bg_order == 0:
            bg_score = bg_model[0, base_idx]
        elif i > bg_order:
            ctx_idx = get_context_idx(indices, i-1, bg_order)
            if ctx_idx < 0:
                bg_score = -2.0
            else:
                bg_score = bg_model[ctx_idx, base_idx]
        else:
            bg_score = -2.0
        
        score = matrix[i-1] + bg_score + base_penalty
        for m in range(num_motifs):
            wml = motif_lengths[m]
            if i >= wml:
                emit_fwd = _scan_at_pos(
                    indices, i - wml, wml, 
                    all_motif_columns_fwd, motif_offsets[m]
                )
                emit_rc = _scan_at_pos(
                    indices, i - wml, wml, 
                    all_motif_columns_rc, motif_offsets[m]
                )
                if emit_fwd > emit_rc:
                    emit = emit_fwd
                else:
                    emit = emit_rc
                    
                if emit > NEG_INF:
                    from_score = matrix[i - wml]
                    path_score = from_score + emit + motif_penalties[m]
                    score = addlog2(score, path_score)
        matrix[i] = score
    
    cdef double result = matrix[idx_len]
    free(matrix)
    return result


# -------------------------------------------------------------
# GAPPED MOTIF SUPPORT
# -------------------------------------------------------------

cdef double _get_bg_score_for_pos(
    long[:] indices,
    int pos,
    double[:, :] bg_model,
    int bg_order
) nogil:
    """Get background score at a single position."""
    cdef long base_idx = indices[pos]
    cdef long ctx_idx
    
    if base_idx < 0:
        return -2.0  # Fallback for N
    
    if bg_order == 0:
        return bg_model[0, base_idx]
    
    if pos >= bg_order:
        ctx_idx = get_context_idx(indices, pos, bg_order)
        if ctx_idx < 0:
            return -2.0
        return bg_model[ctx_idx, base_idx]
    
    return -2.0


cdef double _scan_gapped_motif_ptr(
    long[:] indices,
    int start_pos,
    double[:, :] all_block_cols,
    long[:] block_offsets,
    long[:] block_lengths,
    long *gap_lengths,
    int num_blocks,
    double[:, :] bg_model,
    int bg_order,
    int seq_total_len
) nogil:
    """
    Score a gapped motif at a fixed position with fixed gap lengths.
    Uses C pointer for gap_lengths to allow nogil operation.
    
    Args:
        indices: Sequence as base indices.
        start_pos: Starting position in sequence.
        all_block_cols: Concatenated block PWM columns (4, total_block_len).
        block_offsets: Offset into all_block_cols for each block.
        block_lengths: Length of each block.
        gap_lengths: C pointer to gap lengths (num_blocks - 1 elements).
        num_blocks: Number of blocks.
        bg_model: Background model.
        bg_order: Background Markov order.
        seq_total_len: Total sequence length for bounds checking.
        
    Returns:
        Log2 score of the gapped motif hit, or NEG_INF if invalid.
    """
    cdef double score = 0.0
    cdef int pos = start_pos
    cdef int b, g, p
    cdef long idx
    cdef double bg_score
    
    for b in range(num_blocks):
        # Score block b under PWM
        for p in range(block_lengths[b]):
            if pos >= seq_total_len:
                return NEG_INF
            idx = indices[pos]
            if idx < 0:
                return NEG_INF
            score += all_block_cols[idx, block_offsets[b] + p]
            pos += 1
        
        # Score gap after block b (hard gap = background score)
        if b < num_blocks - 1:
            for g in range(gap_lengths[b]):
                if pos >= seq_total_len:
                    return NEG_INF
                bg_score = _get_bg_score_for_pos(indices, pos, bg_model, bg_order)
                score += bg_score
                pos += 1
    
    return score


cdef int _gapped_motif_total_span(
    long[:] block_lengths,
    long[:] gap_lengths,
    int num_blocks
) nogil:
    """Calculate total span of a gapped motif."""
    cdef int span = 0
    cdef int b
    for b in range(num_blocks):
        span += block_lengths[b]
        if b < num_blocks - 1:
            span += gap_lengths[b]
    return span


def batch_likelihood_gapped(
    long[:, :] all_indices,
    long[:] seq_lengths,
    double[:, :] all_block_columns,
    long[:] block_offsets,
    long[:] block_lengths,
    long[:, :] gap_length_ranges,  # (num_gaps, 2) for min/max per gap
    double[:] gap_log_priors,      # Pre-computed log priors for each gap length
    long num_gap_lengths,          # Number of discrete gap lengths to marginalize
    double motif_penalty,
    double base_penalty,
    double[:, :] bg_model,
    int bg_order
):
    """
    Batch likelihood with gapped motifs (forward strand only).
    
    Marginalizes over all valid gap length combinations.
    
    For simplicity, this version assumes a single gapped motif.
    """
    cdef int num_seqs = all_indices.shape[0]
    cdef int num_blocks = block_offsets.shape[0]
    cdef double total_likelihood = 0.0
    cdef int i
    cdef double seq_hood
    
    with nogil:
        for i in range(num_seqs):
            seq_hood = _compute_single_seq_gapped(
                all_indices[i],
                seq_lengths[i],
                all_block_columns,
                block_offsets,
                block_lengths,
                gap_length_ranges,
                gap_log_priors,
                num_gap_lengths,
                motif_penalty,
                base_penalty,
                num_blocks,
                bg_model,
                bg_order
            )
            total_likelihood += seq_hood
    
    return total_likelihood


cdef double _compute_single_seq_gapped(
    long[:] indices,
    int seq_len,
    double[:, :] all_block_columns,
    long[:] block_offsets,
    long[:] block_lengths,
    long[:, :] gap_length_ranges,
    double[:] gap_log_priors,
    long num_gap_lengths,
    double motif_penalty,
    double base_penalty,
    int num_blocks,
    double[:, :] bg_model,
    int bg_order
) nogil:
    """
    DP for single gapped motif with marginalization over gap lengths.
    
    This handles a single gapped motif. The DP sums over all valid gap 
    configurations at each position.
    """
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL:
        return 0.0
    
    cdef int i, g, h
    cdef double score, emit, from_score, path_score, bg_score
    cdef long base_idx, ctx_idx
    cdef int min_span, span
    cdef long *gap_lens = NULL
    cdef int num_gaps = num_blocks - 1
    
    # Cache gap range values into local C variables
    cdef long gap0_min = 0, gap0_max = 0, gap1_min = 0, gap1_max = 0
    cdef long bl0 = 0, bl1 = 0, bl2 = 0
    
    if num_blocks >= 1:
        bl0 = block_lengths[0]
    if num_blocks >= 2:
        bl1 = block_lengths[1]
    if num_blocks >= 3:
        bl2 = block_lengths[2]
    
    if num_gaps >= 1:
        gap0_min = gap_length_ranges[0, 0]
        gap0_max = gap_length_ranges[0, 1]
    if num_gaps >= 2:
        gap1_min = gap_length_ranges[1, 0]
        gap1_max = gap_length_ranges[1, 1]
    
    # Allocate gap lengths array for iteration
    if num_gaps > 0:
        gap_lens = <long *> malloc(num_gaps * sizeof(long))
        if gap_lens == NULL:
            free(matrix)
            return 0.0
    
    matrix[0] = 0.0
    
    # Calculate minimum span (all gaps at minimum)
    min_span = 0
    for i in range(num_blocks):
        min_span += block_lengths[i]
    if num_gaps >= 1:
        min_span += gap0_min
    if num_gaps >= 2:
        min_span += gap1_min
    
    for i in range(1, seq_len + 1):
        # Background transition
        base_idx = indices[i - 1]
        if base_idx < 0:
            bg_score = -2.0
        elif bg_order == 0:
            bg_score = bg_model[0, base_idx]
        elif i > bg_order:
            ctx_idx = get_context_idx(indices, i - 1, bg_order)
            if ctx_idx < 0:
                bg_score = -2.0
            else:
                bg_score = bg_model[ctx_idx, base_idx]
        else:
            bg_score = -2.0
        
        score = matrix[i - 1] + bg_score + base_penalty
        
        # Gapped motif emission - marginalize over gap lengths
        if i >= min_span:
            if num_gaps == 0:
                # No gaps - simple contiguous motif
                span = min_span
                if i >= span:
                    emit = _scan_gapped_motif_ptr(
                        indices, i - span,
                        all_block_columns, block_offsets, block_lengths,
                        NULL, num_blocks, bg_model, bg_order, seq_len
                    )
                    if emit > NEG_INF:
                        from_score = matrix[i - span]
                        path_score = from_score + emit + motif_penalty
                        score = addlog2(score, path_score)
            
            elif num_gaps == 1:
                # Single gap - enumerate gap lengths
                for g in range(gap0_min, gap0_max + 1):
                    gap_lens[0] = g
                    span = bl0 + g + bl1
                    if i >= span:
                        emit = _scan_gapped_motif_ptr(
                            indices, i - span,
                            all_block_columns, block_offsets, block_lengths,
                            gap_lens, num_blocks, bg_model, bg_order, seq_len
                        )
                        if emit > NEG_INF:
                            # Add gap prior
                            emit += gap_log_priors[g - gap0_min]
                            from_score = matrix[i - span]
                            path_score = from_score + emit + motif_penalty
                            score = addlog2(score, path_score)
            
            elif num_gaps == 2:
                # Two gaps - enumerate both
                for g in range(gap0_min, gap0_max + 1):
                    gap_lens[0] = g
                    for h in range(gap1_min, gap1_max + 1):
                        gap_lens[1] = h
                        span = bl0 + g + bl1 + h + bl2
                        if i >= span:
                            emit = _scan_gapped_motif_ptr(
                                indices, i - span,
                                all_block_columns, block_offsets, block_lengths,
                                gap_lens, num_blocks, bg_model, bg_order, seq_len
                            )
                            if emit > NEG_INF:
                                # Add gap priors
                                emit += gap_log_priors[g - gap0_min]
                                emit += gap_log_priors[h - gap1_min]
                                from_score = matrix[i - span]
                                path_score = from_score + emit + motif_penalty
                                score = addlog2(score, path_score)
        
        matrix[i] = score
    
    cdef double result = matrix[seq_len]
    free(matrix)
    if gap_lens != NULL:
        free(gap_lens)
    return result




def batch_likelihood_gapped_dual(
    long[:, :] all_indices,
    long[:] seq_lengths,
    double[:, :] all_block_columns_fwd,
    double[:, :] all_block_columns_rc,
    long[:] block_offsets,
    long[:] block_lengths,
    long[:, :] gap_length_ranges,
    double[:] gap_log_priors,
    long num_gap_lengths,
    double motif_penalty,
    double base_penalty,
    double[:, :] bg_model,
    int bg_order
):
    """
    Batch likelihood with gapped motifs (both strands).
    
    Takes max of forward and reverse complement scores.
    """
    cdef int num_seqs = all_indices.shape[0]
    cdef int num_blocks = block_offsets.shape[0]
    cdef double total_likelihood = 0.0
    cdef int i
    cdef double seq_hood
    
    with nogil:
        for i in range(num_seqs):
            seq_hood = _compute_single_seq_gapped_dual(
                all_indices[i],
                seq_lengths[i],
                all_block_columns_fwd,
                all_block_columns_rc,
                block_offsets,
                block_lengths,
                gap_length_ranges,
                gap_log_priors,
                num_gap_lengths,
                motif_penalty,
                base_penalty,
                num_blocks,
                bg_model,
                bg_order
            )
            total_likelihood += seq_hood
    
    return total_likelihood


cdef double _compute_single_seq_gapped_dual(
    long[:] indices,
    int seq_len,
    double[:, :] all_block_columns_fwd,
    double[:, :] all_block_columns_rc,
    long[:] block_offsets,
    long[:] block_lengths,
    long[:, :] gap_length_ranges,
    double[:] gap_log_priors,
    long num_gap_lengths,
    double motif_penalty,
    double base_penalty,
    int num_blocks,
    double[:, :] bg_model,
    int bg_order
) nogil:
    """DP for gapped motif with both strands."""
    cdef double *matrix = <double *> malloc((seq_len + 1) * sizeof(double))
    if matrix == NULL:
        return 0.0
    
    cdef int i, g, h
    cdef double score, emit_fwd, emit_rc, emit, from_score, path_score, bg_score
    cdef long base_idx, ctx_idx
    cdef int min_span, span
    cdef long *gap_lens = NULL
    cdef int num_gaps = num_blocks - 1
    
    # Cache gap range values into local C variables
    cdef long gap0_min = 0, gap0_max = 0, gap1_min = 0, gap1_max = 0
    cdef long bl0 = 0, bl1 = 0, bl2 = 0
    
    if num_blocks >= 1:
        bl0 = block_lengths[0]
    if num_blocks >= 2:
        bl1 = block_lengths[1]
    if num_blocks >= 3:
        bl2 = block_lengths[2]
    
    if num_gaps >= 1:
        gap0_min = gap_length_ranges[0, 0]
        gap0_max = gap_length_ranges[0, 1]
    if num_gaps >= 2:
        gap1_min = gap_length_ranges[1, 0]
        gap1_max = gap_length_ranges[1, 1]
    
    if num_gaps > 0:
        gap_lens = <long *> malloc(num_gaps * sizeof(long))
        if gap_lens == NULL:
            free(matrix)
            return 0.0
    
    matrix[0] = 0.0
    
    # Calculate minimum span
    min_span = 0
    for i in range(num_blocks):
        min_span += block_lengths[i]
    if num_gaps >= 1:
        min_span += gap0_min
    if num_gaps >= 2:
        min_span += gap1_min
    
    for i in range(1, seq_len + 1):
        # Background transition
        base_idx = indices[i - 1]
        if base_idx < 0:
            bg_score = -2.0
        elif bg_order == 0:
            bg_score = bg_model[0, base_idx]
        elif i > bg_order:
            ctx_idx = get_context_idx(indices, i - 1, bg_order)
            if ctx_idx < 0:
                bg_score = -2.0
            else:
                bg_score = bg_model[ctx_idx, base_idx]
        else:
            bg_score = -2.0
        
        score = matrix[i - 1] + bg_score + base_penalty
        
        # Gapped motif emission - both strands
        if i >= min_span:
            if num_gaps == 0:
                span = min_span
                if i >= span:
                    emit_fwd = _scan_gapped_motif_ptr(
                        indices, i - span,
                        all_block_columns_fwd, block_offsets, block_lengths,
                        NULL, num_blocks, bg_model, bg_order, seq_len
                    )
                    emit_rc = _scan_gapped_motif_ptr(
                        indices, i - span,
                        all_block_columns_rc, block_offsets, block_lengths,
                        NULL, num_blocks, bg_model, bg_order, seq_len
                    )
                    emit = emit_fwd if emit_fwd > emit_rc else emit_rc
                    if emit > NEG_INF:
                        from_score = matrix[i - span]
                        path_score = from_score + emit + motif_penalty
                        score = addlog2(score, path_score)
            
            elif num_gaps == 1:
                for g in range(gap0_min, gap0_max + 1):
                    gap_lens[0] = g
                    span = bl0 + g + bl1
                    if i >= span:
                        emit_fwd = _scan_gapped_motif_ptr(
                            indices, i - span,
                            all_block_columns_fwd, block_offsets, block_lengths,
                            gap_lens, num_blocks, bg_model, bg_order, seq_len
                        )
                        emit_rc = _scan_gapped_motif_ptr(
                            indices, i - span,
                            all_block_columns_rc, block_offsets, block_lengths,
                            gap_lens, num_blocks, bg_model, bg_order, seq_len
                        )
                        emit = emit_fwd if emit_fwd > emit_rc else emit_rc
                        if emit > NEG_INF:
                            emit += gap_log_priors[g - gap0_min]
                            from_score = matrix[i - span]
                            path_score = from_score + emit + motif_penalty
                            score = addlog2(score, path_score)
            
            elif num_gaps == 2:
                for g in range(gap0_min, gap0_max + 1):
                    gap_lens[0] = g
                    for h in range(gap1_min, gap1_max + 1):
                        gap_lens[1] = h
                        span = bl0 + g + bl1 + h + bl2
                        if i >= span:
                            emit_fwd = _scan_gapped_motif_ptr(
                                indices, i - span,
                                all_block_columns_fwd, block_offsets, block_lengths,
                                gap_lens, num_blocks, bg_model, bg_order, seq_len
                            )
                            emit_rc = _scan_gapped_motif_ptr(
                                indices, i - span,
                                all_block_columns_rc, block_offsets, block_lengths,
                                gap_lens, num_blocks, bg_model, bg_order, seq_len
                            )
                            emit = emit_fwd if emit_fwd > emit_rc else emit_rc
                            if emit > NEG_INF:
                                emit += gap_log_priors[g - gap0_min]
                                emit += gap_log_priors[h - gap1_min]
                                from_score = matrix[i - span]
                                path_score = from_score + emit + motif_penalty
                                score = addlog2(score, path_score)
        
        matrix[i] = score
    
    cdef double result = matrix[seq_len]
    free(matrix)
    if gap_lens != NULL:
        free(gap_lens)
    return result


def generate_markov_sequence(long length, np.ndarray[DTYPE_t, ndim=2] transition, 
                             np.ndarray[DTYPE_t, ndim=1] init_probs, int order,
                             np.ndarray[DTYPE_t, ndim=1] rand_vals):
    """
    Fast Markov sequence generation using pre-generated random numbers.
    
    Args:
        length: Length of sequence to generate.
        transition: (4^order, 4) transition probabilities.
        init_probs: (4,) initial probabilities.
        order: Markov order.
        rand_vals: (length,) array of random floats [0, 1).
        
    Returns:
        Generated sequence string.
    """
    cdef int i, j, ctx_idx, next_base
    cdef double p, cumsum
    cdef char* seq_buf = <char*> malloc(length + 1)
    cdef double[:, :] trans_view = transition
    cdef double[:] init_view = init_probs
    cdef double[:] rand_view = rand_vals
    
    # 4 bases mapping: 0=A, 1=C, 2=G, 3=T
    cdef char[4] bases = [65, 67, 71, 84] # 'A', 'C', 'G', 'T'
    cdef int[4] base_map_rev = [0, 1, 2, 3] # simple mapping for internal logic
    
    # Generate initial context (first 'order' bases)
    # Actually, we use init_probs for the first 'order' bases?
    # Usually standard is: 
    # - 1st base: init_probs (unconditional)
    # - 2nd base: prob(b2|b1) if order >= 1
    # - ...
    # But BackgroundGenerator.MarkovBackground simplified implementation uses 
    # 'init_probs' (IID) for the first 'order' bases. We match that logic.
    
    for i in range(min(length, order)):
        p = rand_view[i]
        cumsum = 0.0
        next_base = 3 # Default T
        for j in range(4):
            cumsum += init_view[j]
            if p < cumsum:
                next_base = j
                break
        seq_buf[i] = bases[next_base]
        
    if length <= order:
        try:
            return seq_buf[:length].decode('ascii')
        finally:
            free(seq_buf)
            
    # Generate rest
    # We need to maintain context as integer index
    # Context index = b_{i-order} * 4^{order-1} + ... + b_{i-1} * 4^0
    
    # Initialize context index from the first 'order' bases
    # We need to reconstruct the index.
    # Note: bases array is 'A','C','G','T'. We need 0,1,2,3 from them.
    # ASCII: A=65, C=67, G=71, T=84
    # Simple switch or map
    
    # Pre-calculate powers of 4? No, simply shift.
    
    cdef int mask = 0
    if order > 0:
        # We need a mask to keep only last 'order' bases in the window
        # But actually context index is modulo 4^order? No.
        # Yes, shift left and mask out high bits usually.
        # Or just: ctx = ctx % (4**(order-1)) * 4 + new_base
        # 4^order might be large (order 5 => 1024).
        
        # Build initial context index
        ctx_idx = 0
        for i in range(order):
            # Map char back to 0-3
            # A->0, C->1, G->2, T->3
            char_code = seq_buf[i]
            val = 0
            if char_code == 67: val = 1
            elif char_code == 71: val = 2
            elif char_code == 84: val = 3
            
            ctx_idx = ctx_idx * 4 + val
            
    cdef int order_pow_minus_1 = 1
    for i in range(order - 1):
        order_pow_minus_1 *= 4
    
    # Main loop
    for i in range(order, length):
        p = rand_view[i]
        
        # Select base based on transition[ctx_idx]
        cumsum = 0.0
        next_base = 3
        for j in range(4):
            cumsum += trans_view[ctx_idx, j]
            if p < cumsum:
                next_base = j
                break
        
        seq_buf[i] = bases[next_base]
        
        # Update context
        if order > 0:
            # Remove oldest base: ctx_idx %= (4^(order-1))
            # But ctx_idx is for 'order' bases.
            # We need to remove the first one and add new one.
            # ctx_idx = (ctx_idx % order_pow_minus_1) * 4 + next_base
            # Wait, 4^order is the size.
            # context is bases [i-order : i].
            # next context is [i-order+1 : i+1].
            # So we drop the high digit.
            if order > 1:
                # Modulo 4^(order-1). But wait.
                # Example order=2. ctx 0..15. bases b1, b2.
                # idx = b1*4 + b2.
                # Next we want b2, b3.
                # idx_new = b2*4 + b3.
                # So we take idx % 4 and... no.
                # We take idx % (4^(order-1)) to drop b1.
                # b1 is the highest digit.
                # 4^(order-1) is the weight of b1? No, 4^(order-1) is weight of highest.
                # So idx_new = (idx % (4**(order-1))) * 4 + new
                pass
            
            if order == 1:
                ctx_idx = next_base
            else:
                # Calculate modulo
                # We need to remove the high part.
                # High part value was b_{oldest} * 4^(order-1)
                # We can just do modulo 4^(order-1) * 4? No.
                # We keep the lower part: modulo 4^(order-1).
                ctx_idx = (ctx_idx % order_pow_minus_1) * 4 + next_base

    try:
        return seq_buf[:length].decode('ascii')
    finally:
        free(seq_buf)
