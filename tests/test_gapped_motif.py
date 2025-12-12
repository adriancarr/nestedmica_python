"""
Unit tests for gapped motif data structures.
"""

import unittest
import numpy as np
from nestedmica.model.gapped_motif import (
    GapConfig, GapSpec, GappedMotif, GapPriorType
)


class TestGapConfig(unittest.TestCase):
    """Tests for GapConfig dataclass."""
    
    def test_default_config(self):
        """Test default configuration."""
        config = GapConfig()
        self.assertFalse(config.allow_gaps)
        self.assertEqual(config.min_gap_length, 0)
        self.assertEqual(config.max_gap_length, 20)
        self.assertEqual(config.max_num_gaps, 2)
        self.assertEqual(config.gap_prior, GapPriorType.GEOMETRIC)
    
    def test_validation(self):
        """Test validation constraints."""
        # min > max should fail
        with self.assertRaises(ValueError):
            GapConfig(min_gap_length=10, max_gap_length=5)
        
        # Negative min should fail
        with self.assertRaises(ValueError):
            GapConfig(min_gap_length=-1)
        
        # Invalid geometric param should fail
        with self.assertRaises(ValueError):
            GapConfig(gap_prior=GapPriorType.GEOMETRIC, gap_prior_param=1.5)
    
    def test_log_prior_uniform(self):
        """Test uniform prior calculation."""
        config = GapConfig(
            min_gap_length=5, 
            max_gap_length=10, 
            gap_prior=GapPriorType.UNIFORM
        )
        # 6 possible values: 5,6,7,8,9,10
        expected = -np.log(6)
        
        for g in range(5, 11):
            self.assertAlmostEqual(config.log_prior_gap_length(g), expected, places=5)
        
        # Out of range should be -inf
        self.assertEqual(config.log_prior_gap_length(4), -np.inf)
        self.assertEqual(config.log_prior_gap_length(11), -np.inf)
    
    def test_log_prior_geometric(self):
        """Test geometric prior prefers shorter gaps."""
        config = GapConfig(
            min_gap_length=1, 
            max_gap_length=10, 
            gap_prior=GapPriorType.GEOMETRIC,
            gap_prior_param=0.3
        )
        
        # Shorter gaps should have higher probability
        p1 = config.log_prior_gap_length(1)
        p5 = config.log_prior_gap_length(5)
        p10 = config.log_prior_gap_length(10)
        
        self.assertGreater(p1, p5)
        self.assertGreater(p5, p10)
    
    def test_sample_gap_length(self):
        """Test gap length sampling."""
        config = GapConfig(
            min_gap_length=3, 
            max_gap_length=7, 
            gap_prior=GapPriorType.UNIFORM
        )
        rng = np.random.default_rng(42)
        
        samples = [config.sample_gap_length(rng) for _ in range(100)]
        
        # All samples should be in valid range
        self.assertTrue(all(3 <= s <= 7 for s in samples))
        # Should see variety
        self.assertGreater(len(set(samples)), 1)


class TestGapSpec(unittest.TestCase):
    """Tests for GapSpec dataclass."""
    
    def test_valid_gap(self):
        """Test valid gap creation."""
        gap = GapSpec(length=5)
        self.assertEqual(gap.length, 5)
    
    def test_invalid_gap(self):
        """Test negative gap length fails."""
        with self.assertRaises(ValueError):
            GapSpec(length=-1)


class TestGappedMotif(unittest.TestCase):
    """Tests for GappedMotif class."""
    
    def setUp(self):
        """Create test PWM blocks."""
        # Simple uniform PWM blocks
        self.block1 = np.full((4, 5), -2.0)  # 5bp block
        self.block2 = np.full((4, 4), -2.0)  # 4bp block
        self.block3 = np.full((4, 3), -2.0)  # 3bp block
    
    def test_single_block(self):
        """Test ungapped motif (single block)."""
        motif = GappedMotif(blocks=[self.block1])
        
        self.assertEqual(motif.num_blocks, 1)
        self.assertEqual(motif.num_gaps, 0)
        self.assertEqual(motif.total_block_length, 5)
        self.assertEqual(motif.total_gap_length, 0)
        self.assertEqual(motif.total_span, 5)
        self.assertFalse(motif.is_gapped)
    
    def test_two_blocks_with_gap(self):
        """Test bipartite motif."""
        gaps = [GapSpec(length=6)]
        motif = GappedMotif(blocks=[self.block1, self.block2], gaps=gaps)
        
        self.assertEqual(motif.num_blocks, 2)
        self.assertEqual(motif.num_gaps, 1)
        self.assertEqual(motif.total_block_length, 9)  # 5 + 4
        self.assertEqual(motif.total_gap_length, 6)
        self.assertEqual(motif.total_span, 15)
        self.assertTrue(motif.is_gapped)
    
    def test_three_blocks_multiple_gaps(self):
        """Test tripartite motif."""
        gaps = [GapSpec(length=3), GapSpec(length=5)]
        motif = GappedMotif(
            blocks=[self.block1, self.block2, self.block3], 
            gaps=gaps
        )
        
        self.assertEqual(motif.num_blocks, 3)
        self.assertEqual(motif.num_gaps, 2)
        self.assertEqual(motif.total_block_length, 12)  # 5 + 4 + 3
        self.assertEqual(motif.total_gap_length, 8)  # 3 + 5
        self.assertEqual(motif.total_span, 20)
    
    def test_to_flat_columns(self):
        """Test flattening for Cython."""
        gaps = [GapSpec(length=5)]
        motif = GappedMotif(blocks=[self.block1, self.block2], gaps=gaps)
        
        columns, offsets, lengths = motif.to_flat_columns()
        
        self.assertEqual(columns.shape, (4, 9))  # 5 + 4 columns
        np.testing.assert_array_equal(offsets, [0, 5])
        np.testing.assert_array_equal(lengths, [5, 4])
    
    def test_get_gap_lengths(self):
        """Test gap length extraction."""
        gaps = [GapSpec(length=3), GapSpec(length=7)]
        motif = GappedMotif(
            blocks=[self.block1, self.block2, self.block3], 
            gaps=gaps
        )
        
        gap_lengths = motif.get_gap_lengths()
        np.testing.assert_array_equal(gap_lengths, [3, 7])
    
    def test_from_contiguous(self):
        """Test creating from contiguous PWM."""
        motif = GappedMotif.from_contiguous(self.block1)
        
        self.assertEqual(motif.num_blocks, 1)
        self.assertFalse(motif.is_gapped)
    
    def test_copy(self):
        """Test deep copy."""
        gaps = [GapSpec(length=5)]
        original = GappedMotif(blocks=[self.block1, self.block2], gaps=gaps)
        copied = original.copy()
        
        # Modify original
        original.blocks[0][0, 0] = 999.0
        original.gaps[0] = GapSpec(length=99)
        
        # Copy should be unchanged
        self.assertNotEqual(copied.blocks[0][0, 0], 999.0)
        self.assertEqual(copied.gaps[0].length, 5)
    
    def test_consensus(self):
        """Test consensus string generation."""
        # Create blocks with known consensus
        block1 = np.array([
            [0.0, -10, -10, -10],  # A at pos 0
            [-10, 0.0, -10, -10],  # C at pos 1
            [-10, -10, 0.0, -10],  # G at pos 2
            [-10, -10, -10, 0.0],  # T at pos 3
        ])
        block2 = np.array([
            [0.0, 0.0],  # A A
            [-10, -10],
            [-10, -10],
            [-10, -10],
        ])
        
        gaps = [GapSpec(length=5)]
        motif = GappedMotif(blocks=[block1, block2], gaps=gaps)
        
        consensus = motif.get_consensus()
        self.assertEqual(consensus, "ACGT[5]AA")
    
    def test_log_prior(self):
        """Test prior probability calculation."""
        config = GapConfig(
            allow_gaps=True,
            min_gap_length=3,
            max_gap_length=10,
            max_num_gaps=2,
            gap_prior=GapPriorType.UNIFORM
        )
        
        gaps = [GapSpec(length=5)]
        motif = GappedMotif(
            blocks=[self.block1, self.block2], 
            gaps=gaps,
            config=config
        )
        
        log_p = motif.log_prior()
        
        # Should include number of gaps prior + gap length prior
        # Number of gaps: uniform over 0,1,2 -> log(1/3)
        # Gap length: uniform over 3-10 (8 values) -> log(1/8)
        expected = -np.log(3) + (-np.log(8))
        self.assertAlmostEqual(log_p, expected, places=5)
    
    def test_validation_block_gap_mismatch(self):
        """Test block/gap count validation."""
        gaps = [GapSpec(length=5), GapSpec(length=3)]  # 2 gaps
        
        # 2 blocks need exactly 1 gap
        with self.assertRaises(ValueError):
            GappedMotif(blocks=[self.block1, self.block2], gaps=gaps)
    
    def test_empty_blocks_fails(self):
        """Test that empty blocks list fails."""
        with self.assertRaises(ValueError):
            GappedMotif(blocks=[])


if __name__ == '__main__':
    unittest.main()


class TestMCMCProposals(unittest.TestCase):
    """Tests for RJMCMC proposal functions."""
    
    def setUp(self):
        """Create test motifs and config."""
        from nestedmica.model.gapped_motif import (
            propose_insert_gap, propose_delete_gap, 
            propose_perturb_gap_length, propose_shift_block_boundary
        )
        self.propose_insert_gap = propose_insert_gap
        self.propose_delete_gap = propose_delete_gap
        self.propose_perturb_gap_length = propose_perturb_gap_length
        self.propose_shift_block_boundary = propose_shift_block_boundary
        
        self.config = GapConfig(
            allow_gaps=True,
            min_gap_length=3,
            max_gap_length=15,
            max_num_gaps=2,
            gap_prior=GapPriorType.UNIFORM
        )
        
        # 12-column single block
        self.single_block = np.full((4, 12), -2.0)
        self.rng = np.random.default_rng(42)
    
    def test_propose_insert_gap_success(self):
        """Test successful gap insertion."""
        motif = GappedMotif(blocks=[self.single_block], gaps=[], config=self.config)
        
        proposed, log_hastings = self.propose_insert_gap(motif, self.rng)
        
        self.assertIsNotNone(proposed)
        self.assertEqual(proposed.num_blocks, 2)
        self.assertEqual(proposed.num_gaps, 1)
        self.assertIsInstance(log_hastings, float)
    
    def test_propose_insert_gap_max_reached(self):
        """Test insertion fails when max gaps reached."""
        block = np.full((4, 6), -2.0)
        gaps = [GapSpec(length=5), GapSpec(length=5)]
        motif = GappedMotif(blocks=[block, block, block], gaps=gaps, config=self.config)
        
        proposed, _ = self.propose_insert_gap(motif, self.rng)
        self.assertIsNone(proposed)
    
    def test_propose_delete_gap_success(self):
        """Test successful gap deletion."""
        block1 = np.full((4, 6), -2.0)
        block2 = np.full((4, 6), -2.0)
        gaps = [GapSpec(length=5)]
        motif = GappedMotif(blocks=[block1, block2], gaps=gaps, config=self.config)
        
        proposed, log_hastings = self.propose_delete_gap(motif, self.rng)
        
        self.assertIsNotNone(proposed)
        self.assertEqual(proposed.num_blocks, 1)
        self.assertEqual(proposed.num_gaps, 0)
        self.assertEqual(proposed.total_block_length, 12)
    
    def test_propose_delete_gap_no_gaps(self):
        """Test deletion fails when no gaps exist."""
        motif = GappedMotif(blocks=[self.single_block], gaps=[], config=self.config)
        
        proposed, _ = self.propose_delete_gap(motif, self.rng)
        self.assertIsNone(proposed)
    
    def test_propose_perturb_gap_length(self):
        """Test gap length perturbation."""
        block1 = np.full((4, 6), -2.0)
        block2 = np.full((4, 6), -2.0)
        gaps = [GapSpec(length=8)]
        motif = GappedMotif(blocks=[block1, block2], gaps=gaps, config=self.config)
        
        # Run multiple times to ensure we get some successful proposals
        successes = 0
        changes = 0
        for _ in range(20):
            proposed, log_hastings = self.propose_perturb_gap_length(motif, self.rng)
            if proposed is not None:
                successes += 1
                self.assertEqual(proposed.num_gaps, 1)
                self.assertEqual(log_hastings, 0.0)  # Symmetric
                if proposed.gaps[0].length != 8:
                    changes += 1
        
        self.assertGreater(successes, 0)
        # At least some proposals should change the gap length
        self.assertGreater(changes, 0)
    
    def test_propose_shift_block_boundary(self):
        """Test block boundary shift."""
        block1 = np.full((4, 8), -2.0)
        block2 = np.full((4, 8), -2.0)
        gaps = [GapSpec(length=8)]
        motif = GappedMotif(blocks=[block1, block2], gaps=gaps, config=self.config)
        
        successes = 0
        for _ in range(20):
            proposed, log_hastings = self.propose_shift_block_boundary(motif, self.rng)
            if proposed is not None:
                successes += 1
                self.assertEqual(proposed.num_gaps, 1)
                # Either block size or gap length should have changed
                total_span = proposed.total_span
                # With shifts, total can change based on gap changes
                self.assertGreater(total_span, 0)
        
        self.assertGreater(successes, 0)
    
    def test_insert_delete_reversibility(self):
        """Test that insert followed by delete can restore original."""
        motif = GappedMotif(blocks=[self.single_block.copy()], gaps=[], config=self.config)
        
        # Insert gap
        inserted, _ = self.propose_insert_gap(motif, self.rng)
        self.assertIsNotNone(inserted)
        self.assertEqual(inserted.num_blocks, 2)
        
        # Delete gap
        restored, _ = self.propose_delete_gap(inserted, self.rng)
        self.assertIsNotNone(restored)
        self.assertEqual(restored.num_blocks, 1)
        self.assertEqual(restored.total_block_length, 12)


if __name__ == '__main__':
    unittest.main()

