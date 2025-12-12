"""
Console utilities for Nested MICA.
"""

from typing import List, Optional

class AsciiPlotter:
    """
    Simple ASCII graph plotter for monitoring numerical trends in the console.
    
    Attributes:
        height (int): Height of the graph in lines.
        width (int): Width of the graph (history buffer size).
        history (List[float]): Buffer of recent values.
    """
    
    def __init__(self, height: int = 10, width: int = 50):
        """
        Initialize the plotter.

        Args:
            height (int): Graph height in text lines.
            width (int): Number of data points to keep in history.
        """
        self.height = height
        self.width = width
        self.history: List[float] = []
        
    def add(self, value: float) -> None:
        """
        Add a value to the plot history.

        Args:
            value (float): The numerical value to plot.
        """
        self.history.append(value)
        if len(self.history) > self.width:
            self.history.pop(0)
            
    def plot(self) -> str:
        """
        Generate the ASCII representation of the current history.

        Returns:
            str: Multi-line string containing the plot.
        """
        if not self.history: 
            return ""
            
        min_v = min(self.history)
        max_v = max(self.history)
        rnge = max_v - min_v if max_v != min_v else 1.0
        
        # Normalize to 0..height-1
        rows = [[' ' for _ in range(len(self.history))] for _ in range(self.height)]
        for x, val in enumerate(self.history):
            normalized = int((val - min_v) / rnge * (self.height - 1))
            # Clamp in case of floating point weirdness
            normalized = max(0, min(self.height - 1, normalized))
            rows[self.height - 1 - normalized][x] = '*'
            
        # Draw axes labels and box
        lines = []
        lines.append(f"    {max_v:.2f}")
        for r in rows:
            lines.append("    |" + "".join(r) + "|")
        lines.append(f"    {min_v:.2f}")
        return "\n".join(lines)
