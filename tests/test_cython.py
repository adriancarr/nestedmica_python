"""
Tests for Cython-optimized functions in cython_model.pyx.
"""

import unittest
import numpy as np
from nestedmica.model.cython_model import (
    CythonWeightMatrix,
    reverse_complement_columns,
    batch_likelihood,
    batch_likelihood_dual
)


class TestCythonWeightMatrix(unittest.TestCase):
    
    def test_creation(self):
        """CythonWeightMatrix should store columns correctly."""
        cols = np.random.random((4, 10))
        wm = CythonWeightMatrix(cols)
        
        self.assertEqual(wm.get_length(), 10)
        np.testing.assert_array_almost_equal(wm.get_columns(), cols)
    
    def test_immutable_copy(self):
        """Modifying original should not affect stored columns."""
        cols = np.random.random((4, 5))
        wm = CythonWeightMatrix(cols.copy())
        
        cols[0, 0] = 999.0  # Modify original
        
        self.assertNotEqual(wm.get_columns()[0, 0], 999.0)


class TestReverseComplement(unittest.TestCase):
    
    def test_rc_identity(self):
        """RC of RC should equal original."""
        cols = np.random.random((4, 8))
        
        rc1 = reverse_complement_columns(cols)
        rc2 = reverse_complement_columns(rc1)
        
        np.testing.assert_array_almost_equal(cols, rc2)
    
    def test_rc_column_order(self):
        """RC reverses column order."""
        cols = np.zeros((4, 3))
        cols[0, 0] = 1.0  # A at position 0
        cols[3, 2] = 1.0  # T at position 2
        
        rc = reverse_complement_columns(cols)
        
        # T at pos 0 becomes A at pos 2 (reversed and complemented)
        self.assertEqual(rc[0, 0], 1.0)  # A at pos 0 (was T at pos 2, complemented)
        self.assertEqual(rc[3, 2], 1.0)  # T at pos 2 (was A at pos 0, complemented)
    
    def test_rc_base_swap(self):
        """RC swaps A<->T and C<->G."""
        # Single position with only A
        cols = np.array([[1.0], [0.0], [0.0], [0.0]])  # 100% A
        
        rc = reverse_complement_columns(cols)
        
        # Should be 100% T now
        np.testing.assert_array_almost_equal(rc, [[0.0], [0.0], [0.0], [1.0]])


class TestBatchLikelihood(unittest.TestCase):
    
    def setUp(self):
        # Simple test setup
        self.indices = np.array([[0, 1, 2, 3, 0, 1, 2, 3]], dtype=np.int64)  # ACGTACGT
        self.lengths = np.array([8], dtype=np.int64)
        
        # Simple motif: uniform probability
        self.cols = np.full((4, 4), -2.0, dtype=np.float64)  # 4 bases x 4 positions
        self.offsets = np.array([0], dtype=np.int64)
        self.motif_lengths = np.array([4], dtype=np.int64)
        self.penalties = np.array([-2.0], dtype=np.float64)
        
        # Background model (order-0, uniform)
        self.bg_model = np.full((1, 4), -2.0, dtype=np.float64)
    
    def test_deterministic(self):
        """Same input should give same output."""
        result1 = batch_likelihood(
            self.indices, self.lengths, self.cols,
            self.offsets, self.motif_lengths, self.penalties,
            -0.01, self.bg_model, 0
        )
        result2 = batch_likelihood(
            self.indices, self.lengths, self.cols,
            self.offsets, self.motif_lengths, self.penalties,
            -0.01, self.bg_model, 0
        )
        
        self.assertEqual(result1, result2)
    
    def test_returns_float(self):
        """Likelihood should be a finite float."""
        result = batch_likelihood(
            self.indices, self.lengths, self.cols,
            self.offsets, self.motif_lengths, self.penalties,
            -0.01, self.bg_model, 0
        )
        
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))


class TestBatchLikelihoodDual(unittest.TestCase):
    
    def setUp(self):
        self.indices = np.array([[0, 1, 2, 3, 0, 1, 2, 3]], dtype=np.int64)
        self.lengths = np.array([8], dtype=np.int64)
        
        self.cols_fwd = np.full((4, 4), -2.0, dtype=np.float64)
        self.cols_rc = reverse_complement_columns(self.cols_fwd)
        self.offsets = np.array([0], dtype=np.int64)
        self.motif_lengths = np.array([4], dtype=np.int64)
        self.penalties = np.array([-2.0], dtype=np.float64)
        
        self.bg_model = np.full((1, 4), -2.0, dtype=np.float64)
    
    def test_dual_geq_single(self):
        """Dual-strand likelihood should be >= single-strand."""
        single = batch_likelihood(
            self.indices, self.lengths, self.cols_fwd,
            self.offsets, self.motif_lengths, self.penalties,
            -0.01, self.bg_model, 0
        )
        dual = batch_likelihood_dual(
            self.indices, self.lengths, self.cols_fwd, self.cols_rc,
            self.offsets, self.motif_lengths, self.penalties,
            -0.01, self.bg_model, 0
        )
        
        # In log space, >= means the value could be equal or greater
        self.assertGreaterEqual(dual, single - 0.01)  # Allow small tolerance


if __name__ == '__main__':
    unittest.main()
