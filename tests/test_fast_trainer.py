import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from nestedmica.trainer.fast import FastTrainer

class TestFastTrainer(unittest.TestCase):
    def setUp(self):
        self.seqs = [
            SeqRecord(Seq("ACGTACGT"), id="1"),
            SeqRecord(Seq("TGCATGCA"), id="2")
        ]
        
    def test_initialization(self):
        trainer = FastTrainer(self.seqs, num_motifs=2, motif_length=5, ensemble_size=10, n_jobs=1)
        self.assertEqual(len(trainer.models), 10)
        self.assertEqual(trainer.num_sequences, 2)
        # Check sequence encoding ('A':0, 'C':1, 'G':2, 'T':3)
        # Seq 1: ACGT... -> 0, 1, 2, 3...
        self.assertEqual(trainer.all_indices[0, 0], 0)
        self.assertEqual(trainer.all_indices[0, 1], 1)
        
    def test_step(self):
        trainer = FastTrainer(self.seqs, num_motifs=1, motif_length=4, ensemble_size=5, n_jobs=1)
        initial_models = len(trainer.models)
        
        worst, hood = trainer.step()
        
        # Ensemble size should remain constant
        self.assertEqual(len(trainer.models), initial_models)
        # Returned worst model should be valid
        self.assertIn('motifs', worst)
        self.assertIn('weights', worst)
        self.assertIsInstance(hood, float)
        
    def test_variable_length_methods(self):
        # Test Delete (Zap)
        # 4 cols -> 3 cols
        cols = np.zeros((4, 4))
        # Zap has min length constraint 5, so let's make it 6
        cols = np.zeros((4, 6))
        
        # We need access to private methods or expose them. 
        # Since they are logically private but practically accessible in Python:
        trainer = FastTrainer(self.seqs, 1, 6, 1, n_jobs=1)
        zapped = trainer._delete_column(cols)
        self.assertEqual(zapped.shape[1], 5)
        
        # Test Insert (Indel)
        # 6 cols -> 7 cols
        indeling = trainer._insert_column(cols)
        self.assertEqual(indeling.shape[1], 7)
    def test_evidence_calculation(self):
        # Create trainer with ensemble size 1 for simple math
        # X_0 = 1, X_1 = e^-1, X_2 = e^-2
        # w_1 = 0.5 * (X_0 - X_2) = 0.5 * (1 - e^-2)
        # log_w1 approx log(0.5) + log(0.86) = -1 + (-0.21) = -1.21 (natural)
        # We work in log2.
        
        trainer = FastTrainer(self.seqs, 1, 5, 10, n_jobs=1)
        self.assertEqual(trainer.log_evidence, -np.float64('inf'))
        self.assertEqual(trainer.step_count, 0)
        
        # Run one step
        worst, hood = trainer.step()
        
        # Evidence should now be finite (log2 sum)
        self.assertGreater(trainer.log_evidence, -1e9) 
        self.assertEqual(trainer.step_count, 1)
        
        # Run another step
        prev_evidence = trainer.log_evidence
        trainer.step()
        self.assertNotEqual(trainer.log_evidence, prev_evidence)

if __name__ == '__main__':
    unittest.main()
