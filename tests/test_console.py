import unittest
from nestedmica.utils.console import AsciiPlotter

class TestAsciiPlotter(unittest.TestCase):
    def test_empty_plot(self):
        p = AsciiPlotter(height=5, width=10)
        self.assertEqual(p.plot(), "")
        
    def test_simple_plot(self):
        p = AsciiPlotter(height=5, width=5)
        for i in range(5):
            p.add(i)
        
        output = p.plot()
        # Verify min/max labels
        self.assertIn("4.00", output)
        self.assertIn("0.00", output)
        # Verify structure
        lines = output.split('\n')
        self.assertEqual(len(lines), 7) # 5 + 2 labels

    def test_overflow(self):
        p = AsciiPlotter(height=5, width=3)
        p.add(1); p.add(2); p.add(3); p.add(4)
        # Should only have 2, 3, 4
        self.assertEqual(len(p.history), 3)
        self.assertEqual(p.history[0], 2)

if __name__ == '__main__':
    unittest.main()
