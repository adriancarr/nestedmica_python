"""
Unit tests for nestedmica.synthetic package.
"""

import unittest
import tempfile
import os
import json
import numpy as np
from Bio import SeqIO

from nestedmica.synthetic.background import (
    UniformBackground, MarkovBackground, LearnedBackground,
    FixedLength, NormalLength, LogNormalLength, LearnedLength
)
from nestedmica.synthetic.motifs import (
    ConsensusMotif, PWMMotif, GappedSyntheticMotif
)
from nestedmica.synthetic.planting import UniformPlanting, PositionalPlanting, PlantedMotif
from nestedmica.synthetic.datasets import (
    SyntheticDataset, SimpleBenchmark, GappedMotifBenchmark
)


class TestLengthDistributions(unittest.TestCase):
    """Tests for sequence length distributions."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
    
    def test_fixed_length(self):
        """Test fixed length distribution."""
        dist = FixedLength(200)
        self.assertEqual(dist.sample(self.rng), 200)
        lengths = dist.sample_n(10, self.rng)
        self.assertTrue(np.all(lengths == 200))
    
    def test_normal_length(self):
        """Test normal length distribution."""
        dist = NormalLength(mean=200, std=30, min_len=100, max_len=300)
        lengths = dist.sample_n(100, self.rng)
        self.assertTrue(np.all(lengths >= 100))
        self.assertTrue(np.all(lengths <= 300))
        # Mean should be roughly 200
        self.assertAlmostEqual(np.mean(lengths), 200, delta=20)
    
    def test_lognormal_length(self):
        """Test log-normal length distribution."""
        dist = LogNormalLength(mean=5.3, sigma=0.3, min_len=50, max_len=500)
        lengths = dist.sample_n(100, self.rng)
        self.assertTrue(np.all(lengths >= 50))
        self.assertTrue(np.all(lengths <= 500))


class TestBackgroundGenerators(unittest.TestCase):
    """Tests for background sequence generators."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
    
    def test_uniform_background(self):
        """Test uniform background with default GC."""
        bg = UniformBackground(gc_content=0.5)
        seq = bg.generate(1000, self.rng)
        self.assertEqual(len(seq), 1000)
        gc = (seq.count('G') + seq.count('C')) / len(seq)
        self.assertAlmostEqual(gc, 0.5, delta=0.1)
    
    def test_uniform_gc_bias(self):
        """Test uniform background with GC bias."""
        bg = UniformBackground(gc_content=0.7)
        seq = bg.generate(1000, self.rng)
        gc = (seq.count('G') + seq.count('C')) / len(seq)
        self.assertAlmostEqual(gc, 0.7, delta=0.1)
    
    def test_markov_background(self):
        """Test Markov background generation."""
        bg = MarkovBackground(order=2, gc_content=0.5)
        seq = bg.generate(500, self.rng)
        self.assertEqual(len(seq), 500)
        self.assertTrue(all(b in 'ACGT' for b in seq))
    
    def test_generate_n(self):
        """Test batch generation."""
        bg = UniformBackground(0.5)
        seqs = bg.generate_n(100, 10, self.rng)
        self.assertEqual(len(seqs), 10)
        self.assertTrue(all(len(s) == 100 for s in seqs))


class TestLearnedBackground(unittest.TestCase):
    """Tests for learned background from FASTA."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
        
        # Create temporary FASTA with balanced composition
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        for i in range(20):
            # Balanced sequence
            seq = 'ACGT' * 37 + 'AC'  # 150 bp, 50% GC
            self.temp_fasta.write(f">seq_{i}\n{seq}\n")
        self.temp_fasta.close()
    
    def tearDown(self):
        os.unlink(self.temp_fasta.name)
    
    def test_learned_from_fasta(self):
        """Test learning background from FASTA."""
        bg = LearnedBackground.from_fasta(self.temp_fasta.name, order=1)
        # Should have roughly 50% GC from balanced input
        self.assertAlmostEqual(bg.gc_content, 0.5, delta=0.15)
    
    def test_generate_sequence(self):
        """Test that generated sequences are valid DNA."""
        bg = LearnedBackground.from_fasta(self.temp_fasta.name, order=2)
        seq = bg.generate(500, self.rng)
        self.assertEqual(len(seq), 500)
        # All bases should be valid
        self.assertTrue(all(b in 'ACGT' for b in seq))


class TestMotifs(unittest.TestCase):
    """Tests for motif classes."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
    
    def test_consensus_motif(self):
        """Test consensus motif sampling."""
        motif = ConsensusMotif("ACGTACGT")
        instance = motif.sample(self.rng)
        self.assertEqual(len(instance), 8)
        # Should be mostly ACGTACGT due to high consensus probability
        self.assertTrue(all(b in 'ACGT' for b in instance))
    
    def test_consensus_with_iupac(self):
        """Test consensus with IUPAC codes."""
        motif = ConsensusMotif("ACGTNNN")
        self.assertEqual(motif.length, 7)
        instance = motif.sample(self.rng)
        self.assertEqual(len(instance), 7)
    
    def test_pwm_motif(self):
        """Test PWM motif."""
        pwm = np.array([
            [0.9, 0.1, 0.1, 0.1, 0.25],
            [0.03, 0.8, 0.1, 0.1, 0.25],
            [0.03, 0.05, 0.7, 0.1, 0.25],
            [0.04, 0.05, 0.1, 0.7, 0.25]
        ])
        motif = PWMMotif(pwm, name="test")
        self.assertEqual(motif.length, 5)
        instance = motif.sample(self.rng)
        self.assertEqual(len(instance), 5)
    
    def test_gapped_motif_from_pattern(self):
        """Test gapped motif pattern parsing."""
        motif = GappedSyntheticMotif.from_pattern("ACGT[3-8]TGCA")
        self.assertEqual(len(motif.blocks), 2)
        self.assertEqual(motif.gap_ranges, [(3, 8)])
        
        instance, gaps = motif.sample(self.rng)
        self.assertTrue(3 <= gaps[0] <= 8)
        # Length = 4 + gap + 4
        self.assertEqual(len(instance), 4 + gaps[0] + 4)


class TestPlanting(unittest.TestCase):
    """Tests for motif planting strategies."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
    
    def test_uniform_planting(self):
        """Test uniform planting."""
        planter = UniformPlanting(density=1.0, min_edge_distance=5)
        seq = "A" * 100
        _, planted = planter.plant(seq, 10, 0, self.rng)
        
        self.assertIsNotNone(planted)
        self.assertGreaterEqual(planted.start, 5)
        self.assertLessEqual(planted.end, 95)
    
    def test_density_zero(self):
        """Test that density=0 plants nothing."""
        planter = UniformPlanting(density=0.0)
        seq = "A" * 100
        _, planted = planter.plant(seq, 10, 0, self.rng)
        self.assertIsNone(planted)
    
    def test_strand_bias(self):
        """Test strand bias."""
        planter = UniformPlanting(density=1.0, strand_bias=1.0)
        seq = "A" * 100
        for _ in range(10):
            _, planted = planter.plant(seq, 10, 0, self.rng)
            if planted:
                self.assertEqual(planted.strand, '+')

    def test_positional_planting(self):
        """Test positional planting bias."""
        # Biased towards center (50) of 100bp seq
        planter = PositionalPlanting(center_fraction=0.5, std_fraction=0.1, density=1.0)
        seq = "A" * 100
        
        positions = []
        for _ in range(50):
            _, planted = planter.plant(seq, 10, 0, self.rng)
            if planted:
                positions.append(planted.start)
        
        # Mean position should be roughly 45 (center 50 - 5 for motif/2? No, start)
        # Center of motif (start+5) should be near 50. So start near 45.
        mean_pos = np.mean(positions)
        self.assertTrue(35 < mean_pos < 55)


class TestShuffling(unittest.TestCase):
    """Tests for sequence shuffling."""
    
    def setUp(self):
        self.rng = np.random.default_rng(42)
        
    def count_dinucleotides(self, seq):
        counts = {}
        for i in range(len(seq) - 1):
            pair = seq[i:i+2]
            counts[pair] = counts.get(pair, 0) + 1
        return counts
    
    def test_dinucl_shuffle_counts(self):
        """Test that dinucleotide shuffle preserves counts."""
        from nestedmica.synthetic.background import dinucl_shuffle
        
        seq = "ACGT" * 20 + "AAAA" + "CCCC"
        shuffled = dinucl_shuffle(seq, self.rng)
        
        counts_orig = self.count_dinucleotides(seq)
        counts_shuff = self.count_dinucleotides(shuffled)
        
        self.assertEqual(counts_orig, counts_shuff)
        self.assertNotEqual(seq, shuffled) # Unlikely to be same
        
    def test_shuffled_background(self):
        """Test ShuffledBackground generator."""
        from nestedmica.synthetic.background import ShuffledBackground
        
        # Create dummy sequences
        # Need a sequence with multiple possible Eulerian paths
        # "AATTAATTAA" allows varying order of AA/TT loops
        seqs = ["AATTAATTAA" * 10] 
        bg = ShuffledBackground(seqs)
        
        # Generate same length
        gen = bg.generate(len(seqs[0]), self.rng)
        self.assertEqual(len(gen), len(seqs[0]))
        # With high probability (almost certainty for length 100), shuffled != original
        self.assertNotEqual(gen, seqs[0])
        
        # Check dinucleotides match source
        counts_orig = self.count_dinucleotides(seqs[0])
        counts_gen = self.count_dinucleotides(gen)
        self.assertEqual(counts_orig, counts_gen)
        
    def test_shuffled_padding(self):
        """Test generation with padding."""
        from nestedmica.synthetic.background import ShuffledBackground
        seqs = ["ACGT"] # 4 bp
        bg = ShuffledBackground(seqs)
        
        # Ask for 10 bp
        gen = bg.generate(10, self.rng)
        self.assertEqual(len(gen), 10)
        # First 4 should match counts
        self.assertEqual(self.count_dinucleotides(gen[:4]), self.count_dinucleotides(seqs[0]))


class TestDatasets(unittest.TestCase):
    """Tests for dataset generators."""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.rng = np.random.default_rng(42)
    
    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_simple_benchmark(self):
        """Test simple benchmark generation."""
        output = os.path.join(self.temp_dir, "simple.fa")
        
        dataset = SimpleBenchmark(
            motif="ACGTACGT",
            n_sequences=20,
            seq_length=100,
            density=1.0,
            seed=42
        )
        dataset.generate(output)
        
        # Check FASTA was created
        self.assertTrue(os.path.exists(output))
        records = list(SeqIO.parse(output, "fasta"))
        self.assertEqual(len(records), 20)
        self.assertEqual(len(records[0].seq), 100)
    
    def test_gapped_benchmark(self):
        """Test gapped motif benchmark."""
        output = os.path.join(self.temp_dir, "gapped.fa")
        truth = os.path.join(self.temp_dir, "gapped_truth.json")
        
        dataset = GappedMotifBenchmark(
            pattern="ACGT[5-10]TGCA",
            n_sequences=10,
            seq_length=150,
            density=1.0,
            seed=42
        )
        dataset.generate(output)
        dataset.save_truth(truth)
        
        # Check output files
        self.assertTrue(os.path.exists(output))
        self.assertTrue(os.path.exists(truth))
        
        # Verify truth structure
        with open(truth) as f:
            data = json.load(f)
        self.assertIn('motifs', data)
        self.assertIn('sequences', data)
    
    def test_variable_length(self):
        """Test variable-length sequence generation."""
        output = os.path.join(self.temp_dir, "variable.fa")
        
        dataset = SyntheticDataset(
            n_sequences=50,
            seq_length=NormalLength(mean=200, std=50),
            seed=42
        )
        dataset.add_motif("GATAAGA")
        dataset.generate(output)
        
        records = list(SeqIO.parse(output, "fasta"))
        lengths = [len(r.seq) for r in records]
        
        # Should have variable lengths
        # Should have variable lengths
        self.assertGreater(max(lengths) - min(lengths), 20)


class TestBenchmarkGenerator(unittest.TestCase):
    """Tests for BenchmarkGenerator."""
    
    def setUp(self):
        self.test_dir = "test_benchmark_output"
        os.makedirs(self.test_dir, exist_ok=True)
        
    def tearDown(self):
        import shutil
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    def test_standard_suite(self):
        from nestedmica.synthetic.datasets import BenchmarkGenerator, SimpleBenchmark
        
        gen = BenchmarkGenerator(self.test_dir)
        
        config = {'motif': 'ACGT'}
        gen.create_standard_suite(config, train_size=10, test_size=5)
        
        # Check structure
        train_fa = os.path.join(self.test_dir, 'train', 'sequences.fa')
        test_fa = os.path.join(self.test_dir, 'test', 'sequences.fa')
        
        self.assertTrue(os.path.exists(train_fa))
        self.assertTrue(os.path.exists(test_fa))
        
        # Verify count using SeqIO
        from Bio import SeqIO
        records = list(SeqIO.parse(train_fa, "fasta"))
        self.assertEqual(len(records), 10)
        
        records = list(SeqIO.parse(test_fa, "fasta"))
        self.assertEqual(len(records), 5)


class TestCollisionDetection(unittest.TestCase):
    """Test multi-motif collision detection."""
    
    def test_avoid_overlap(self):
        """Test that multiple planted motifs do not overlap."""
        from nestedmica.synthetic.datasets import SyntheticDataset
        from nestedmica.synthetic.background import UniformBackground, FixedLength
        from nestedmica.synthetic.planting import UniformPlanting
        
        # Sequence length 20. Motifs are 8bp each.
        # If we plant 2 motifs, they consume 16bp. Space is tight.
        dataset = SyntheticDataset(
            n_sequences=50,
            seq_length=FixedLength(20),
            background=UniformBackground(0.5),
            seed=42
        )
        
        dataset.add_motif("AAAAAAAA") # 8bp
        dataset.add_motif("CCCCCCCC") # 8bp
        
        dataset.set_planting(UniformPlanting(density=1.0, avoid_overlap=True))
        
        # We need to capture the ground truth which is usually saved to file 
        # or accessible via dataset.ground_truth after generation.
        # generate() writes to file but also populates self.ground_truth.
        
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.fa') as tmp:
            dataset.generate(tmp.name)
            
            # Check collisions in ground truth
            for item in dataset.ground_truth:
                planted = item['planted']
                # Should have up to 2 motifs. Check overlap.
                intervals = []
                for p in planted:
                    intervals.append((p['start'], p['end']))
                
                # Check all pairs
                for i in range(len(intervals)):
                    for j in range(i + 1, len(intervals)):
                        s1, e1 = intervals[i]
                        s2, e2 = intervals[j]
                        # Overlap condition: start1 < end2 AND start2 < end1
                        overlap = (s1 < e2) and (s2 < e1)
                        if overlap:
                            self.fail(f"Found overlap: {s1}-{e1} vs {s2}-{e2} in {item['id']}")


if __name__ == '__main__':
    unittest.main()
