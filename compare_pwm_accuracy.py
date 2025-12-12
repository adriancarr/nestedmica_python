import re
import numpy as np

def parse_xms(filename):
    """
    Rudimentary XMS parser to extract probability matrices.
    Returns list of dicts: {'id': int, 'matrix': np.array (4 x L)}
    """
    motifs = []
    current_matrix = []
    in_matrix = False
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if "<weightMatrix" in line:
                in_matrix = True
                current_matrix = []
            elif "</weightMatrix>" in line:
                in_matrix = False
                if current_matrix:
                    # Shape is L x 4 in file usually, need to check
                    # The file format shows rows for columns? No, looking at previous cat:
                    # It has multiple lines, each with 4 floats.
                    # So it's L lines of 4 floats.
                    mat = np.array(current_matrix).T # Turn into 4 x L
                    motifs.append(mat)
            elif in_matrix:
                parts = line.split()
                if len(parts) == 4:
                    vals = [float(x) for x in parts]
                    current_matrix.append(vals)
    return motifs

def get_consensus_and_avg_prob(matrix):
    """
    Returns consensus string and average probability of the dominant bases.
    Map: 0=A, 1=C, 2=G, 3=T
    """
    bases = ['A', 'C', 'G', 'T']
    consensus = []
    probs = []
    
    for i in range(matrix.shape[1]):
        col = matrix[:, i]
        idx = np.argmax(col)
        consensus.append(bases[idx])
        probs.append(col[idx])
        
    return "".join(consensus), np.mean(probs)

def main():
    print("Parsing complex_varying_result.xms...")
    motifs = parse_xms("complex_varying_result.xms")
    
    # Ground Truth Data
    # Name: (Sequence, Expected Probability of Dominant Base)
    # Expected Prob = Conservation + (1-Conservation)/4
    targets = {
        'Strong': ('ACGTACGT', 0.95 + 0.0125),   # 0.9625
        'Medium': ('TTTGGCCAAA', 0.80 + 0.05),   # 0.8500
        'Weak':   ('TATATATATATATA', 0.60 + 0.10) # 0.7000
    }
    
    print("\n" + "="*60)
    print(f"{'Motif':<10} | {'Expected':<12} | {'Found':<12} | {'Exp. Prob':<10} | {'Found Prob':<10} | {'Status'}")
    print("-" * 60)
    
    found_map = {}
    
    # Simple matching logic
    for m in motifs:
        seq, avg_prob = get_consensus_and_avg_prob(m)
        
        # Identify
        label = "Unknown"
        match_stats = None
        
        # Check against targets (allowing for substring matching for the Weak one)
        for name, (t_seq, t_prob) in targets.items():
            if t_seq in seq or seq in t_seq:
                if name == 'Weak':
                    # Clean up the weak one 'tATAT...' -> 'ATATA...' for checking
                    pass 
                label = name
                match_stats = (t_seq, t_prob)
                break
        
        if label != "Unknown":
            t_seq, t_prob = match_stats
            print(f"{label:<10} | {len(t_seq):<12} | {len(seq):<12} | {t_prob:.4f}     | {avg_prob:.4f}     | {'✅' if abs(t_prob - avg_prob) < 0.05 else '⚠️'}")
        else:
            print(f"Unknown    | {'?':<12} | {len(seq):<12} | {'?':<10} | {avg_prob:.4f}     | ❓")

if __name__ == "__main__":
    main()
