"""
Checkpointing utilities for Nested MICA.
"""

import pickle
import os
from typing import Any, Dict, List, Tuple
from nestedmica.model.cython_model import CythonWeightMatrix
# Circular dependency avoidance: We will import FastTrainer inside load_checkpoint if needed,
# or we assume the caller handles the class instantiation if we just return state.
# Ideally, load_checkpoint should re-instantiate the trainer.
# To avoid circular imports, we can pass the class as an argument or import inside the function.

def save_checkpoint(trainer: Any, filename: str, cycle: int) -> None:
    """
    Save trainer state to file.

    Args:
        trainer (FastTrainer): The trainer instance to save.
        filename (str): Path to save the checkpoint pickle.
        cycle (int): Current cycle number.
    """
    # Extract raw data from Cython objects
    serialized_models = []
    for model in trainer.models:
        motifs_data = [m.get_columns().copy() for m in model['motifs']]
        serialized_models.append({'motifs': motifs_data, 'weights': model['weights'].copy()})
        
    state = {
        'cycle': cycle,
        'models': serialized_models,
        'model_likelihoods': trainer.model_likelihoods,
        'num_motifs': trainer.num_motifs,
        'motif_length': trainer.motif_length,
        'ensemble_size': trainer.ensemble_size,
        'log_evidence': trainer.log_evidence,
        'step_count': trainer.step_count
    }
    
    with open(filename, 'wb') as f:
        pickle.dump(state, f)
    print(f"  [Checkpoint saved to {filename}]")

def load_checkpoint(filename: str, sequences: List[Any], n_jobs: int, trainer_cls: Any) -> Tuple[Any, int]:
    """
    Load trainer from checkpoint.

    Args:
        filename (str): Path to the checkpoint pickle.
        sequences (List): List of BioPython sequences (or whatever the trainer expects).
        n_jobs (int): Number of threads to use.
        trainer_cls (class): The trainer class to instantiate (e.g., FastTrainer).

    Returns:
        Tuple[FastTrainer, int]: Reconstructed trainer and the cycle number.
    """
    print(f"Loading checkpoint from {filename}...")
    with open(filename, 'rb') as f:
        state = pickle.load(f)
        
    trainer = trainer_cls(sequences, state['num_motifs'], state['motif_length'], 
                         state['ensemble_size'], n_jobs=n_jobs)
    
    # Restore models
    trainer.models = []
    for s_model in state['models']:
        motifs = [CythonWeightMatrix(cols) for cols in s_model['motifs']]
        trainer.models.append({'motifs': motifs, 'weights': s_model['weights']})
        
    trainer.model_likelihoods = state['model_likelihoods']
    
    # Restore Evidence state if present (backward compatibility)
    if 'log_evidence' in state:
        trainer.log_evidence = state['log_evidence']
        trainer.step_count = state['step_count']
        
    return trainer, state['cycle']
