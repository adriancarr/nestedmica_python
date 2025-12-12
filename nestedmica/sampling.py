
import numpy as np
from nestedmica.maths import log2 #, addLog2

class ContributionPrior:
    def probability(self, item):
        raise NotImplementedError
    
    def variate(self):
        raise NotImplementedError

class ContributionSampler:
    def sample(self, item, uncle_vector=None):
        raise NotImplementedError

class PenalizedVariate:
    def __init__(self, variate, balance_penalty=0.0, silent=False):
        self.variate = variate
        self.balance_penalty = balance_penalty
        self.silent = silent

class SimplePrior(ContributionPrior):
    # Valid for what? Let's assume WeightedWeightMatrix for now.
    def probability(self, item):
        return 0.0 # Log prob
        
    def variate(self):
        # Return a random item
        # Placeholder
        return None

class MixPolicy:
    def sample(self, mixture):
        # mixture is numpy array
        # Mutate it
        idx = np.random.randint(len(mixture))
        mixture[idx] = np.random.random() # Logic needed here
        
    def sample_component(self, mixture, component):
        pass

    def prior(self, mixture):
        return 0.0

    def variate(self, mixture):
        mixture[:] = np.random.random(mixture.shape)
        
class BinaryMixPolicy(MixPolicy):
    def __init__(self, expected_usage_fraction=0.5):
        self.expected_usage_fraction = expected_usage_fraction
        
    def variate(self, mixture):
        # mixture is numpy array (N,)
        # Set values to 0 or 1 based on expected usage
        r = np.random.random(len(mixture))
        mixture[:] = (r < self.expected_usage_fraction).astype(float)
        
    def sample(self, mixture):
        # Sample one component
        idx = np.random.randint(len(mixture))
        # Flip it? Or sample from prior?
        # Java BinaryMixPolicy.sample: flips a random bit.
        current = mixture[idx]
        mixture[idx] = 1.0 - current


# Specific for PWMs
class WeightMatrixPrior(ContributionPrior):
    def __init__(self, length, alphabet_size=4):
        self.length = length
        self.alphabet_size = alphabet_size
        
    def probability(self, wwm):
        # wwm is WeightedWeightMatrix
        # Return log prior
        # Uniform prior: const
        return 0.0
        
    def variate(self):
        # Generate random PWM
        from nestedmica.model.motif import WeightMatrix, WeightedWeightMatrix
        # cols = np.random.dirichlet([1.0]*self.alphabet_size, size=self.length).T
        # My WeightMatrix impl used (L, 4).
        cols = np.random.random((self.length, self.alphabet_size))
        cols /= cols.sum(axis=1, keepdims=True)
        
        # Convert to log-probabilities for internal use
        log_cols = log2(cols + 1e-10) # 1e-10 for numerical stability
        
        return WeightedWeightMatrix(WeightMatrix(log_cols), 1.0) # Weight 1.0

class SimplePWMSampler(ContributionSampler):
    def sample(self, item, uncle_vector=None):
        # item is WeightedWeightMatrix
        # Perturb it
        from nestedmica.model.motif import WeightMatrix, WeightedWeightMatrix
        
        orig_wm = item.get_weight_matrix()
        cols = orig_wm.get_columns().copy()
        
        # Pick a column and re-sample from Dirichlet
        col_idx = np.random.randint(len(cols))
        
        # New column (prob space)
        new_col_probs = np.random.dirichlet([1.0]*4)
        new_col_log = log2(new_col_probs + 1e-10)
        
        cols[col_idx] = new_col_log
        
        new_wm = WeightMatrix(cols, orig_wm.offset)
        new_wwm = WeightedWeightMatrix(new_wm, item.get_weight())
        
        return PenalizedVariate(new_wwm)
