"""
Tests for the higher-order Markov background model.
"""

import unittest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nestedmica.model.background import (
    learn_markov_background, 
    context_to_index, 
    print_background_stats
)


class TestBackgroundModel(unittest.TestCase):
    
    def test_order_0_uniform(self):
        """Order-0 on uniform sequence should give ~equal probabilities."""
        # Create uniform sequence
        seqs = [SeqRecord(Seq("ACGT" * 25), id="1")]  # 100bp, uniform
        
        bg = learn_markov_background(seqs, order=0)
        
        # Should be 1 context (order-0)
        self.assertEqual(bg.shape[0], 1)
        self.assertEqual(bg.shape[1], 4)
        
        # All probs should be similar (log2 of ~0.25 = -2)
        for i in range(4):
            self.assertAlmostEqual(bg[0, i], -2.0, delta=0.1)
    
    def test_order_3_shape(self):
        """Order-3 should produce 64 contexts."""
        seqs = [SeqRecord(Seq("ACGT" * 100), id="1")]
        
        bg = learn_markov_background(seqs, order=3)
        
        # 4^3 = 64 contexts
        self.assertEqual(bg.shape[0], 64)
        self.assertEqual(bg.shape[1], 4)
    
    def test_gc_bias_detection(self):
        """GC-rich sequence should have higher GC probabilities."""
        # Create GC-rich sequence
        gc_seq = "GCGCGCGC" * 50  # 400bp, 100% GC
        seqs = [SeqRecord(Seq(gc_seq), id="1")]
        
        bg = learn_markov_background(seqs, order=0)
        
        # G and C (indices 1, 2) should have higher log probs than A, T (0, 3)
        gc_prob = (np.power(2, bg[0, 1]) + np.power(2, bg[0, 2]))  # Linear scale
        at_prob = (np.power(2, bg[0, 0]) + np.power(2, bg[0, 3]))
        
        self.assertGreater(gc_prob, at_prob)
    
    def test_context_to_index(self):
        """Verify base-4 encoding of contexts."""
        base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        # AAA = 0*16 + 0*4 + 0 = 0
        self.assertEqual(context_to_index("AAA", base_map), 0)
        
        # TTT = 3*16 + 3*4 + 3 = 48 + 12 + 3 = 63
        self.assertEqual(context_to_index("TTT", base_map), 63)
        
        # ACG = 0*16 + 1*4 + 2 = 6
        self.assertEqual(context_to_index("ACG", base_map), 6)
    
    def test_handles_n_bases(self):
        """Background model should handle N bases gracefully."""
        seqs = [SeqRecord(Seq("ACGTNACGT"), id="1")]
        
        # Should not raise
        bg = learn_markov_background(seqs, order=0)
        
        # Should have valid shape
        self.assertEqual(bg.shape, (1, 4))


class TestPrintBackgroundStats(unittest.TestCase):
    
    def test_no_crash(self):
        """print_background_stats should run without error."""
        seqs = [SeqRecord(Seq("ACGT" * 25), id="1")]
        bg = learn_markov_background(seqs, order=3)
        
        # Should not raise
        print_background_stats(bg, order=3)


if __name__ == '__main__':
    unittest.main()
