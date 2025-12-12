"""
Tests for motif export functions.
"""

import unittest
import tempfile
import os
import numpy as np
from nestedmica.model.cython_model import CythonWeightMatrix
from nestedmica.utils.export import (
    export_xms, export_meme, export_pfm, export_transfac
)


class MockMotif:
    """Mock motif for testing export functions."""
    def __init__(self, length=8):
        # Create uniform PWM in log2 scale
        self.columns = np.full((4, length), -2.0, dtype=np.float64)
        
    def get_length(self):
        return self.columns.shape[1]
    
    def get_columns(self):
        return self.columns


class TestExportXms(unittest.TestCase):
    
    def test_creates_file(self):
        """XMS export should create file with correct structure."""
        motifs = [MockMotif(8)]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xms', delete=False) as f:
            filepath = f.name
        
        try:
            export_xms(motifs, filepath, log_evidence=-1000.0)
            
            with open(filepath, 'r') as f:
                content = f.read()
            
            self.assertIn('<motifs>', content)
            self.assertIn('</motifs>', content)
            self.assertIn('GlobalLogEvidence', content)
            self.assertIn('weightMatrix', content)
        finally:
            os.unlink(filepath)


class TestExportMeme(unittest.TestCase):
    
    def test_creates_meme_file(self):
        """MEME export should create valid MEME format."""
        motifs = [MockMotif(8)]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.meme', delete=False) as f:
            filepath = f.name
        
        try:
            export_meme(motifs, filepath)
            
            with open(filepath, 'r') as f:
                content = f.read()
            
            self.assertIn('MEME version', content)
            self.assertIn('ALPHABET=', content)
            self.assertIn('MOTIF', content)
            self.assertIn('letter-probability matrix', content)
        finally:
            os.unlink(filepath)


class TestExportPfm(unittest.TestCase):
    
    def test_creates_pfm_file(self):
        """PFM export should create JASPAR format."""
        motifs = [MockMotif(8)]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pfm', delete=False) as f:
            filepath = f.name
        
        try:
            export_pfm(motifs, filepath)
            
            with open(filepath, 'r') as f:
                content = f.read()
            
            self.assertIn('>', content)  # Header line
            self.assertIn('A [', content)
            self.assertIn('C [', content)
            self.assertIn('G [', content)
            self.assertIn('T [', content)
        finally:
            os.unlink(filepath)


class TestExportTransfac(unittest.TestCase):
    
    def test_creates_transfac_file(self):
        """TRANSFAC export should create valid format."""
        motifs = [MockMotif(8)]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            filepath = f.name
        
        try:
            export_transfac(motifs, filepath)
            
            with open(filepath, 'r') as f:
                content = f.read()
            
            self.assertIn('AC', content)
            self.assertIn('ID', content)
            self.assertIn('P0', content)
            self.assertIn('//', content)
        finally:
            os.unlink(filepath)


if __name__ == '__main__':
    unittest.main()
