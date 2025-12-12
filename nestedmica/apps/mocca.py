
import sys
import argparse
import numpy as np

from nestedmica.seq import load_fasta
from nestedmica.model.core import FacetteMap, ContributionGroup, Datum
from nestedmica.model.motif_generative import MotifFacette
from nestedmica.model.motif import WeightMatrix, WeightedWeightMatrix
from nestedmica.sampling import WeightMatrixPrior, SimplePWMSampler, BinaryMixPolicy
from nestedmica.trainer import RandomDecopTrainer

def main():
    parser = argparse.ArgumentParser(description="NestedMICA Python Port - Mocca")
    parser.add_argument("-seqs", required=True, help="FASTA file of sequences")
    parser.add_argument("-numMotifs", type=int, default=1, help="Number of motifs to find")
    parser.add_argument("-motifLength", type=int, default=10, help="Length of motifs")
    parser.add_argument("-ensembleSize", type=int, default=50, help="Size of nested sampling ensemble")
    parser.add_argument("-out", default="motifs.xms", help="Output file")
    parser.add_argument("-maxCycles", type=int, default=1000, help="Max cycles")
    
    args = parser.parse_args()
    
    print(f"Loading sequences from {args.seqs}...")
    seq_db = load_fasta(args.seqs)
    print(f"Loaded {len(seq_db)} sequences.")
    
    # Setup Data
    # Convert seq_db to Datum[]
    # In Mocca.java, Datum holds objects for each facette.
    # Here we have 1 facette (Sequence).
    data_list = []
    for seq_id in seq_db.ids():
        seq_rec = seq_db.get_sequence(seq_id)
        # Facetted data is list of objects.
        # We have 1 facette.
        data_list.append(Datum(seq_id, [seq_rec]))
    
    data = data_list
    
    # Setup Model Structure
    facette = MotifFacette(uncounted_expectation=1.0)
    cg = ContributionGroup("motifs")
    facette_map = FacetteMap([cg], [facette])
    # Mapping logic in core.py: get_contribution_for_facette returns first group. Correct.
    
    # Setup Priors and Samplers
    # We have 1 contribution group.
    # Prior for PWMs
    prior = WeightMatrixPrior(length=args.motifLength)
    # Sampler
    sampler = SimplePWMSampler()
    
    mix_policy = BinaryMixPolicy(expected_usage_fraction=0.5)
    
    print("Initializing trainer...")
    trainer = RandomDecopTrainer(
        facette_map,
        data,
        args.numMotifs,
        [prior],
        [sampler],
        mix_policy,
        args.ensembleSize
    )
    
    trainer.init()
    
    print("Starting sampling...")
    for i in range(args.maxCycles):
        # Run one step
        discarded, weight = trainer.next_step()
        
        if i % 100 == 0:
            best_model = trainer.models[-1] # Approximation
            hood = best_model.likelihood()
            print(f"Cycle {i}: L={hood:.2f} Evidence={weight:.2f}")
            
    print("Finished.")
    
    # Save results
    best_trainable = trainer.models[-1]
    # Find true best
    best_hood = -np.inf
    best_model_state = None
    for ts in trainer.models:
        if ts.likelihood() > best_hood:
            best_hood = ts.likelihood()
            best_model_state = ts
            
    print(f"Best Likelihood: {best_hood}")
    
    # Dump motifs
    # best_model_state.get_contributions(cg) -> list of ContributionItems
    contribs = best_model_state.get_contributions(cg)
    with open(args.out, "w") as f:
        f.write("<motifs>\n")
        for i, item in enumerate(contribs):
            wwm = item.get_item()
            wm = wwm.get_weight_matrix()
            cols = wm.get_columns() # (L, 4) or (4, L)
            f.write(f"  <motif id='{i}' weight='{wwm.get_weight()}'>\n")
            f.write(f"    <weightMatrix length='{wm.get_length()}'>\n")
            # Write columns
            for row in cols:
                # wm is now in log-space, convert back to probabilities for XMS
                f.write("      " + " ".join(f"{2**x:.4f}" for x in row) + "\n")
            f.write("    </weightMatrix>\n")
            f.write("  </motif>\n")
        f.write("</motifs>\n")
    print(f"Saved motifs to {args.out}")

if __name__ == "__main__":
    main()
