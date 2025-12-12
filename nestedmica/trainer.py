
import numpy as np
import copy
from nestedmica.maths import log2, addLog2, exp2, addLog
from nestedmica.model.core import SimpleMultiICAModel, ContributionItem, SimpleContributionItem

class TrainableState:
    def __init__(self, model):
        self.model = model # SimpleMultiICAModel
        self.checkpoints = []
        
    def checkpoint(self):
        # Shallow checkpoint: only save references to mixing_matrix and contributions
        # These are the only things that change during decorrelation
        checkpoint_data = {
            'mixing_matrix': self.model.mixing_matrix.copy(),  # Shallow copy of array
            'contributions': {}
        }
        # Deep copy only the contributions (small objects)
        for cg, items in self.model.contributions.items():
            checkpoint_data['contributions'][cg] = items.copy()  # Shallow list copy
        self.checkpoints.append(checkpoint_data)
        
    def rollback(self):
        if self.checkpoints:
            data = self.checkpoints.pop()
            self.model.mixing_matrix = data['mixing_matrix']
            self.model.contributions = data['contributions']
        else:
            raise Exception("No checkpoint to rollback to")
            
    def commit(self):
        if self.checkpoints:
            self.checkpoints.pop() # Discard backup
        
    # Delegate methods
    def likelihood(self):
        return self.model.likelihood()
        
    def get_mixture(self, d):
        return self.model.get_mixture(d)
        
    def get_contributions(self, cg):
        return self.model.get_contributions(cg)
    
    def set_contribution_maybe_clean(self, cg, c, item, silent):
        self.checkpoint()
        self.model.set_contribution(cg, c, item)
        # Note: Java calls checkpoint implicitly? No, caller calls commit/rollback.
        # Actually my logic above:
        # Caller loop:
        #   change model
        #   if bad: rollback
        #   else: commit
        # So I shouldn't checkpoint inside setter unless implicit.
        # Java RandomDecopTrainer:
        #   newModel.setContributionMaybeClean(...)
        #   if (accept) commit else rollback.
        # This implies setContribution DOES NOT checkpoint itself, OR it does and returns.
        # Wait, Java `TrainableState` handles transactions.
        # If I look at `RandomDecopTrainer.java`:
        #   newModel.permuteContributions(...)
        #   if (...) commit else rollback
        # This implies `permuteContributions` starts a transaction.
        
        # Let's assume explicit transaction management in `decorrelateState` for now:
        # self.checkpoint() -> make change -> likelihood() -> commit/rollback
        pass

class Trainer:
    def __init__(self, facette_map, data, components, priors, samplers, mix_policy, ensemble_size):
        self.facette_map = facette_map
        self.data = data
        self.components = components
        self.priors = priors
        self.samplers = samplers
        self.mix_policy = mix_policy
        self.ensemble_size = ensemble_size
        
        self.models = [] # List of TrainableState
        self.model_likelihoods = []
        
        self.step = 0
        self.facettes = facette_map.get_facettes()
        self.cgs = facette_map.get_contribution_groups()
        
    def init(self):
        # Initialize ensemble
        for _ in range(self.ensemble_size):
            model = self.direct_sample_model()
            self.models.append(model)
            self.model_likelihoods.append(model.likelihood())
            
    def direct_sample_model(self):
        # Create a fresh model from prior
        model_impl = SimpleMultiICAModel(self.facette_map, self.data, self.components)
        
        # Sample mixture
        for d in range(len(self.data)):
            self.mix_policy.variate(model_impl.get_mixture(d))
            
        # Sample contributions
        for i, cg in enumerate(self.cgs):
            contribs = model_impl.get_contributions(cg)
            for c in range(self.components):
                # Sample from prior
                item = self.priors[i].variate()
                model_impl.set_contribution(cg, c, SimpleContributionItem(item))
                
        return TrainableState(model_impl)
        
    def next_step(self):
        self.step += 1
        
        # 1. Find worst model
        min_idx = np.argmin(self.model_likelihoods)
        min_hood = self.model_likelihoods[min_idx]
        worst_model = self.models.pop(min_idx)
        _ = self.model_likelihoods.pop(min_idx)
        
        # 2. Re-sample (clone random survivor and evolve)
        survivor_idx = np.random.randint(len(self.models))
        survivor = self.models[survivor_idx]
        
        # Clone it
        new_model_impl = copy.deepcopy(survivor.model)
        new_state = TrainableState(new_model_impl)
        
        # Decorrelate
        # We need to ensure new likelihod > min_hood
        # RandomDecopTrainer logic essentially does MCMC to move it
        # But we need to start with at least min_hood. 
        # Since we cloned a survivor, its likelihood is > min_hood (by definition of sorted removal, except if all equal).
        
        current_hood = survivor.likelihood() # Should be same as cached
        
        # Evolve
        new_hood = self.decorrelate_state(new_state, min_hood)
        
        self.models.append(new_state)
        self.model_likelihoods.append(new_hood)
        
        # Return "WeightedModel" (just the discarded one and its weight)
        log_width = np.log(self.ensemble_size) - np.log(self.ensemble_size + 1)
        log_weight = min_hood + self.step * log_width # roughly
        
        return worst_model, log_weight

    def decorrelate_state(self, state, min_likelihood):
        raise NotImplementedError

class RandomDecopTrainer(Trainer):
    def decorrelate_state(self, state, min_likelihood):
        # Simple Metropolis walk restricted to L > min_likelihood
        
        # Steps - increased from 20 to 100 to match Java's adaptive ~200 steps
        steps = 100
        
        current_hood = state.likelihood()
        
        for _ in range(steps):
            state.checkpoint()
            
            # Propose move
            # 50% mix, 50% contribution
            if np.random.random() < 0.5:
                # Mixture move
                idx = np.random.randint(len(self.data))
                self.mix_policy.sample(state.get_mixture(idx))
            else:
                # Contribution move
                cg_idx = np.random.randint(len(self.cgs))
                cg = self.cgs[cg_idx]
                c = np.random.randint(self.components)
                
                # Get current item
                # Uncle vector? omitted for simplicity
                
                # Sample
                item = state.get_contributions(cg)[c]
                # item is ContributionItem, get_item() gives content
                
                pv = self.samplers[cg_idx].sample(item.get_item())
                # pv.variate is new item
                
                state.model.set_contribution(cg, c, SimpleContributionItem(pv.variate))
                
            new_hood = state.likelihood()
            
            if new_hood > min_likelihood:
                state.commit()
                current_hood = new_hood
            else:
                state.rollback()
                
        return current_hood
