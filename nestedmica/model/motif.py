
import numpy as np

class WeightMatrix:
    def __init__(self, columns, offset=0):
        # columns is a numpy array of shape (L, AlphabetSize)
        # e.g. (10, 4) for DNA
        self.columns = columns
        self.offset = offset
        
    def get_columns(self):
        return self.columns
        
    def get_length(self):
        return self.columns.shape[0]

    def get_score(self, sequence_window):
        # sequence_window: integer array of indices
        # simple score: sum of logs
        # This assumes columns are log-probabilities or weights
        # If columns are probabilities, we need to log them.
        # NestedMICA usually works with weights directly or log-odds.
        rows = np.arange(self.get_length())
        return np.sum(self.columns[rows, sequence_window])

class NMWeightMatrix(WeightMatrix):
    # Same as WeightMatrix but explicit about offset logic if needed
    pass

class WeightedWeightMatrix:
    def __init__(self, wm, weight):
        self.wm = wm
        self.weight = weight
        
    def get_weight_matrix(self):
        return self.wm
        
    def get_weight(self):
        return self.weight
