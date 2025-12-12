import numpy as np
import random

def generate_huge(filename, n_seqs=50000, seq_len=50, n_motifs=10, motif_len=10):
    print(f"Generating {n_seqs} sequences of length {seq_len}...")
    
    # Random background letters
    bases = np.array(['A', 'C', 'G', 'T'])
    # Fast generation using numpy integers then mapping
    # 0-3
    ints = np.random.randint(0, 4, size=(n_seqs, seq_len))
    
    # Generate random consensus motifs
    motifs = []
    for _ in range(n_motifs):
        m_ints = np.random.randint(0, 4, size=motif_len)
        motifs.append(m_ints)
        
    # Embed motifs
    # For each sequence, embed 1 random motif? Or all?
    # To mimic real biological data:
    # Each sequence has ~1 motif occurrence. Some 0, some >1.
    # We'll embed 1 random motif per sequence to ensure signal.
    
    embed_count = 0
    for i in range(n_seqs):
        # Pick random motif
        m_idx = np.random.randint(0, n_motifs)
        m = motifs[m_idx]
        
        # Pick random pos
        pos = np.random.randint(0, seq_len - motif_len)
        
        # Embed
        ints[i, pos:pos+motif_len] = m
        embed_count += 1
        
    print(f"Embedded {embed_count} motif instances.")
    
    # Write to FASTA
    # Optimize writing
    print(f"Writing to {filename}...")
    with open(filename, 'w') as f:
        for i in range(n_seqs):
            # Convert ints to string
            # Fast way:
            seq_chars = bases[ints[i]]
            seq_str = "".join(seq_chars)
            f.write(f">seq_{i}\n{seq_str}\n")
            
    print("Done.")

if __name__ == "__main__":
    generate_huge("huge_test.fa")
