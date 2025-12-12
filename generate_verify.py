import numpy as np
import random

def generate_verification(filename, n_seqs=100, seq_len=50):
    motif = "ACGTACGTAC" # Length 10
    print(f"Generating {n_seqs} sequences with embedded motif '{motif}'...")
    
    with open(filename, 'w') as f:
        for i in range(n_seqs):
            # Random background
            seq = "".join(random.choices("ACGT", k=seq_len))
            # Embed
            pos = random.randint(0, seq_len - len(motif))
            seq_list = list(seq)
            for j, char in enumerate(motif):
                seq_list[pos+j] = char
            seq_final = "".join(seq_list)
            
            f.write(f">seq_{i}\n{seq_final}\n")
    print("Done.")

if __name__ == "__main__":
    generate_verification("verify.fa")
