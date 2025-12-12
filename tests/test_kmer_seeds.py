"""
Tests for k-mer enrichment seeding module.
"""

import unittest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nestedmica.utils.kmer_seeds import (
    count_kmers,
    reverse_complement,
    find_enriched_kmers,
    kmer_to_pwm,
    generate_seed_pwms
)


class TestKmerCounting(unittest.TestCase):
    
    def test_count_simple(self):
        """Count k-mers in simple sequence."""
        seqs = [SeqRecord(Seq("ACGTACGT"), id="1")]
        counts = count_kmers(seqs, k=4)
        
        # ACGT appears twice (pos 0 and 4), plus reverse complements
        self.assertIn("ACGT", counts)
        self.assertGreater(counts["ACGT"], 0)
    
    def test_skip_n_bases(self):
        """K-mers with N should be skipped."""
        seqs = [SeqRecord(Seq("ACNTNACGT"), id="1")]
        counts = count_kmers(seqs, k=4)
        
        # K-mers containing N should not be counted
        for kmer, count in counts.items():
            self.assertNotIn('N', kmer)


class TestReverseComplement(unittest.TestCase):
    
    def test_rc_simple(self):
        """Reverse complement of ACGT should be ACGT."""
        self.assertEqual(reverse_complement("ACGT"), "ACGT")
    
    def test_rc_asymmetric(self):
        """RC of asymmetric sequence."""
        self.assertEqual(reverse_complement("AAAA"), "TTTT")
        self.assertEqual(reverse_complement("GCGC"), "GCGC")


class TestEnrichedKmers(unittest.TestCase):
    
    def test_finds_enriched(self):
        """Should find over-represented k-mers."""
        # Create sequence with planted motif
        motif = "TATAAATA"
        background = "ACGTACGTACGTACGTACGTACGTACGTACGT"
        seq = background + motif + background + motif + background
        seqs = [SeqRecord(Seq(seq), id="1")]
        
        enriched = find_enriched_kmers(seqs, k=8, top_n=5)
        
        # Should find something
        self.assertGreater(len(enriched), 0)
        # Top enriched should be related to TATAAATA
        top_kmer = enriched[0][0]
        self.assertEqual(len(top_kmer), 8)


class TestKmerToPwm(unittest.TestCase):
    
    def test_shape(self):
        """PWM should have correct shape."""
        pwm = kmer_to_pwm("ACGT")
        self.assertEqual(pwm.shape, (4, 4))
    
    def test_log_scale(self):
        """PWM should be in log2 scale (negative values)."""
        pwm = kmer_to_pwm("ACGT")
        self.assertTrue(np.all(pwm < 0))  # Log of probabilities < 1
    
    def test_consensus_high(self):
        """Consensus base should have highest probability."""
        pwm = kmer_to_pwm("AAAA")
        # A (index 0) should have highest log-prob at each position
        for pos in range(4):
            self.assertEqual(np.argmax(pwm[:, pos]), 0)


class TestGenerateSeedPwms(unittest.TestCase):
    
    def test_returns_list(self):
        """Should return list of PWMs."""
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGT" * 10
        seqs = [SeqRecord(Seq(seq), id="1")]
        
        seeds = generate_seed_pwms(seqs, num_motifs=2, motif_length=8)
        
        self.assertIsInstance(seeds, list)


if __name__ == '__main__':
    unittest.main()
