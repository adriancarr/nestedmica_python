
import numpy as np

# matrix.py could just be this:
# from numpy import ndarray as Matrix1D
# from numpy import ndarray as Matrix2D
# But for now we just use numpy directly.

class ContributionGroup:
    def __init__(self, name):
        self.name = name
    
    def __hash__(self):
        return hash(self.name)
        
    def __eq__(self, other):
        return self.name == other.name

class Facette:
    def get_likelihood_calculator(self, data):
        raise NotImplementedError

class LikelihoodCalculator:
    def likelihood(self, contributions, mixture):
        raise NotImplementedError

class FacetteMap:
    def __init__(self, contribution_groups, facettes):
        self.contribution_groups = contribution_groups
        self.facettes = facettes
        # Simple mapping for now, Java had more complex mapping
        
    def get_contribution_groups(self):
        return self.contribution_groups
        
    def get_facettes(self):
        return self.facettes

    def get_contribution_for_facette(self, facette):
        # Simplified: assume 1-to-1 or just return the first group
        # In Java: SimpleFacetteMap had a specific mapping.
        # For MOCCA, there is usually 1 group ("motifs") and 1 facette per seqDB.
        return self.contribution_groups[0] 

class Datum:
    def __init__(self, id, facetted_data):
        self.id = id
        self.facetted_data = facetted_data
        
    def get_facetted_data(self):
        return self.facetted_data

class MultiICAModel:
    def get_components(self):
        raise NotImplementedError
        
    def get_mixture(self, datum_index):
        raise NotImplementedError
        
    def get_contributions(self, group):
        raise NotImplementedError
    
    def likelihood(self):
        raise NotImplementedError

class SimpleMultiICAModel(MultiICAModel):
    def __init__(self, facette_map, data, components, mix_matrix=None, contributions=None):
        self.facette_map = facette_map
        self.data = data
        self.components = components
        
        if mix_matrix is None:
            self.mixing_matrix = np.zeros((len(data), components))
        else:
            self.mixing_matrix = mix_matrix
            
        if contributions is None:
            self.contributions = {}
            for cg in facette_map.get_contribution_groups():
                # Store as list of items or numpy array?
                # Java uses ObjectMatrix1D containing ContributionItems
                # Let's use a list of ContributionItems for now.
                self.contributions[cg] = [None] * components
        else:
            self.contributions = contributions

    def get_components(self):
        return self.components
        
    def get_mixture(self, datum_index):
        return self.mixing_matrix[datum_index]
        
    def get_contributions(self, group):
        return self.contributions[group]
        
    def set_contribution(self, group, component, item):
        self.contributions[group][component] = item
        
    def likelihood(self):
        facettes = self.facette_map.get_facettes()
        L = 0.0
        for d, datum in enumerate(self.data):
            mixture = self.mixing_matrix[d]
            for f, facette in enumerate(facettes):
                facet_data = datum.get_facetted_data()[f]
                calc = facette.get_likelihood_calculator(facet_data)
                cg = self.facette_map.get_contribution_for_facette(facette)
                # In python implementation, we assume we can pass the contributions and mixture
                L += calc.likelihood(self.get_contributions(cg), mixture)
        return L

class ContributionItem:
    def get_item(self):
        raise NotImplementedError

class SimpleContributionItem(ContributionItem):
    def __init__(self, item):
        self.item = item
        
    def get_item(self):
        return self.item
