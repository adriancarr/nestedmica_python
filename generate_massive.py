import numpy as np
import random

def generate_massive(filename, n_seqs=500000, seq_len=50, n_motifs=10, motif_len=10):
    print(f"Generating {n_seqs} sequences of length {seq_len}...")
    
    # Random background
    bases = np.array(['A', 'C', 'G', 'T'])
    # 500k * 50 * 1 bytes = 25MB ints. Fast.
    ints = np.random.randint(0, 4, size=(n_seqs, seq_len), dtype=np.int8)
    
    # Generate motifs
    motifs = []
    for _ in range(n_motifs):
        # motif ints
        motifs.append(np.random.randint(0, 4, size=motif_len, dtype=np.int8))
        
    print(f"Embedding motifs...")
    # Vectorized embedding? 
    # Iterate seqs is slow in Python.
    # 500k iterations is 0.5s. Fine.
    
    # Pre-select motifs and positions
    motif_indices = np.random.randint(0, n_motifs, size=n_seqs)
    positions = np.random.randint(0, seq_len - motif_len, size=n_seqs)
    
    # Loop to embed (vectorizing this is tricky without fancy indexing)
    # Just loop. 500k is fine.
    for i in range(n_seqs):
        m = motifs[motif_indices[i]]
        p = positions[i]
        ints[i, p:p+motif_len] = m
    
    print(f"Writing to {filename}...")
    # Fast I/O
    with open(filename, 'w') as f:
        for i in range(n_seqs):
            # Convert to string
            row = ints[i]
            # Fast mapping
            chars = "".join([bases[x] for x in row])
            f.write(f">seq_{i}\n{chars}\n")
            
    print("Done.")

if __name__ == "__main__":
    generate_massive("massive_test.fa")
