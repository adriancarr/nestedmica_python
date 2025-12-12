"""
Tests for CLI applications (mocca_fast, discover, synthetic).
"""

import unittest
import subprocess
import sys
import os
import shutil
import contextlib
import io
from unittest.mock import patch

# Import synthetic logic for in-process testing
from nestedmica.synthetic.cli import main as synthetic_main


class TestMoccaFastCLI(unittest.TestCase):
    
    def test_help_exits_zero(self):
        """mocca_fast --help should exit 0."""
        result = subprocess.run(
            [sys.executable, '-m', 'nestedmica.apps.mocca_fast', '-h'],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn('FASTA', result.stdout)
    
    def test_missing_required_args(self):
        """mocca_fast with no args should exit non-zero."""
        result = subprocess.run(
            [sys.executable, '-m', 'nestedmica.apps.mocca_fast'],
            capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)


class TestDiscoverCLI(unittest.TestCase):
    
    def test_help_exits_zero(self):
        """discover --help should exit 0."""
        result = subprocess.run(
            [sys.executable, '-m', 'nestedmica.apps.discover', '-h'],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn('Auto-Discovery', result.stdout)
    
    def test_missing_required_args(self):
        """discover with no args should exit non-zero."""
        result = subprocess.run(
            [sys.executable, '-m', 'nestedmica.apps.discover'],
            capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)


class TestSyntheticCLI(unittest.TestCase):
    """Integration tests for the synthetic data CLI."""
    
    def setUp(self):
        self.test_dir = "test_cli_output"
        os.makedirs(self.test_dir, exist_ok=True)
        self.held_stdout = io.StringIO()
        
    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    def run_cli(self, args):
        """Run synthetic CLI with arguments in-process."""
        # Patch sys.argv and stdout
        # Note: sys.argv[0] is usually script name, argparse ignores it or uses it as prog
        argv = ['nestedmica.synthetic'] + args
        
        with patch('sys.argv', argv):
            with patch('sys.stdout', self.held_stdout):
                try:
                    synthetic_main()
                    return 0
                except SystemExit as e:
                    return e.code if isinstance(e.code, int) else 1
                except Exception:
                    return 1

    def test_generate_simple(self):
        """Test simple generation."""
        output = os.path.join(self.test_dir, "simple.fa")
        code = self.run_cli([
            'generate', 
            '-o', output, 
            '-n', '10', 
            '--motif', 'ACGT'
        ])
        
        self.assertEqual(code, 0) # run_cli returns 0 on success
        self.assertTrue(os.path.exists(output))
        truth = output.replace('.fa', '_truth.json')
        self.assertTrue(os.path.exists(truth))

    def test_generate_positional(self):
        """Test positional planting option."""
        output = os.path.join(self.test_dir, "pos.fa")
        code = self.run_cli([
            'generate', 
            '-o', output, 
            '-n', '10', 
            '--motif', 'ACGT',
            '--planting', 'center'
        ])
        self.assertTrue(os.path.exists(output))
        
    def test_benchmark_command(self):
        """Test benchmark subcommand."""
        output_dir = os.path.join(self.test_dir, "bench")
        code = self.run_cli([
            'benchmark',
            '-o', output_dir,
            '--motif', 'ACGT',
            '--train', '10',
            '--test', '5'
        ])
        
        self.assertTrue(os.path.exists(os.path.join(output_dir, 'train', 'sequences.fa')))
        self.assertTrue(os.path.exists(os.path.join(output_dir, 'test', 'sequences.fa')))

    def test_invalid_bg_method(self):
        """Test usage of invalid choices."""
        output = os.path.join(self.test_dir, "bad.fa")
        # Capture stderr to suppress argparse errors
        with contextlib.redirect_stderr(io.StringIO()):
            try:
                code = self.run_cli([
                    'generate', 
                    '-o', output, 
                    '--bg-method', 'invalid_method'
                ])
                # Should have raised SystemExit(2)
                self.assertNotEqual(code, 0)
            except SystemExit:
                pass


if __name__ == '__main__':
    unittest.main()
