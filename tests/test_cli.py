"""
Tests for CLI applications.
"""

import unittest
import subprocess
import sys


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


if __name__ == '__main__':
    unittest.main()
