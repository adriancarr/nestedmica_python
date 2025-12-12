#!/usr/bin/env python3
"""
inspect_checkpoint.py - View state of a running NestedMICA job.
Usage: python3 inspect_checkpoint.py complex_deep.pkl [--out snapshot.xms]
"""

import pickle
import argparse
import numpy as np
import sys

def main():
    parser = argparse.ArgumentParser(description='Inspect NestedMICA Checkpoint')
    parser.add_argument('checkpoint', help='Pickle file to inspect')
    parser.add_argument('--out', help='Optional: Write best model to XMS file')
    args = parser.parse_args()
    
    print(f"Loading {args.checkpoint}...")
    try:
        with open(args.checkpoint, 'rb') as f:
            state = pickle.load(f)
    except FileNotFoundError:
        print("Error: Checkpoint file not found.")
        sys.exit(1)
        
    cycle = state['cycle']
    likelihoods = np.array(state['model_likelihoods'])
    best_idx = np.argmax(likelihoods)
    best_hood = likelihoods[best_idx]
    
    # Calculate IQR
    q75, q25 = np.percentile(likelihoods, [75, 25])
    iqr = q75 - q25
    
    print("-" * 40)
    print(f"Status at Cycle {cycle}")
    print("-" * 40)
    print(f"Ensemble Size:   {state['ensemble_size']}")
    if 'log_evidence' in state:
        print(f"Log Evidence:    {state['log_evidence']:.2f}")
    print(f"Best Likelihood: {best_hood:.2f}")
    print(f"Likelihood IQR:  {iqr:.4f}")    print(f"Mean Likelihood: {np.mean(likelihoods):.2f}")
    print("-" * 40)
    
    if args.out:
        best_model = state['models'][best_idx]
        with open(args.out, 'w') as f:
            f.write("<motifs>\n")
            f.write(f"<!-- Snapshot Cycle: {cycle}, Best Likelihood: {best_hood} -->\n")
            
            # The 'motifs' in pickle are raw numpy arrays (cols), not objects if saved by our refactored saver
            # Actually, let's check nestedmica.utils.checkpoint.save_checkpoint:
            # motifs_data = [m.get_columns().copy() for m in model['motifs']]
            # So they are List[np.ndarray]
            
            for i, cols in enumerate(best_model['motifs']):
                length = cols.shape[1]
                f.write(f"  <motif id='{i}' weight='1.0'>\n")
                f.write(f"    <weightMatrix length='{length}'>\n")
                
                for pos in range(length):
                    probs = np.power(2.0, cols[:, pos])
                    probs = probs / probs.sum()
                    f.write(f"      {probs[0]:.4f} {probs[1]:.4f} {probs[2]:.4f} {probs[3]:.4f}\n")
                    
                f.write("    </weightMatrix>\n")
                f.write("  </motif>\n")
            f.write("</motifs>\n")
        print(f"Saved best model snapshot to {args.out}")

if __name__ == "__main__":
    main()
