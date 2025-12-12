import unittest
import os
import shutil
from nestedmica.utils.checkpoint import save_checkpoint, load_checkpoint
from nestedmica.trainer.fast import FastTrainer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class MockTrainer:
    def __init__(self):
        self.models = [{'motifs': [], 'weights': []}]
        self.model_likelihoods = [0.0]
        self.num_motifs = 1
        self.motif_length = 10
        self.ensemble_size = 1
        
class TestCheckpoint(unittest.TestCase):
    def setUp(self):
        self.test_file = "test_cp.pkl"
        
    def tearDown(self):
        if os.path.exists(self.test_file):
            os.remove(self.test_file)
            
    def test_save_load(self):
        # Create a real trainer to test deep structure
        seqs = [SeqRecord(Seq("ACGT"), id="1")]
        # We need to mock n_jobs=1 to avoid thread pool in tests
        trainer = FastTrainer(seqs, num_motifs=1, motif_length=5, ensemble_size=2, n_jobs=1)
        
        # Modify state to verify persistence
        trainer.model_likelihoods[0] = -123.45
        
        save_checkpoint(trainer, self.test_file, cycle=99)
        
        # Load back
        loaded_trainer, cycle = load_checkpoint(self.test_file, seqs, n_jobs=1, trainer_cls=FastTrainer)
        
        self.assertEqual(cycle, 99)
        self.assertEqual(loaded_trainer.num_motifs, 1)
        self.assertEqual(loaded_trainer.motif_length, 5)
        self.assertEqual(loaded_trainer.model_likelihoods, trainer.model_likelihoods)
        self.assertEqual(len(loaded_trainer.models), 2)

if __name__ == '__main__':
    unittest.main()
