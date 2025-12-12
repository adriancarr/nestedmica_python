"""
Generate a test dataset where the motif is planted ONLY on the reverse strand.
If RC handling works, the algorithm should still find the motif.
If not, it will fail to find anything significant.
"""
import numpy as np
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def reverse_complement(seq):
    """Simple reverse complement."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(comp[b] for b in reversed(seq))

def generate():
    random.seed(123)
    np.random.seed(123)
    
    # The "true" motif (as it would appear on forward strand)
    # We will plant its REVERSE COMPLEMENT in the sequences
    TRUE_MOTIF = "ACGTACGT"  # 8bp, very conserved
    RC_MOTIF = reverse_complement(TRUE_MOTIF)  # "ACGTACGT" -> "ACGTACGT" (palindrome!)
    
    # Actually let's use a non-palindromic motif for a real test
    TRUE_MOTIF = "GGCCAATTCC"  # 10bp
    RC_MOTIF = reverse_complement(TRUE_MOTIF)  # "GGAATTGGCC"
    
    print(f"True motif (forward): {TRUE_MOTIF}")
    print(f"RC motif (planted):   {RC_MOTIF}")
    
    seqs = []
    for i in range(100):
        # Random background
        L = 150
        bases = ['A', 'C', 'G', 'T']
        s_list = [random.choice(bases) for _ in range(L)]
        
        # Plant the RC version (90% of sequences)
        if random.random() < 0.9:
            pos = random.randint(20, 100)
            for j, b in enumerate(RC_MOTIF):
                s_list[pos + j] = b
        
        s_str = "".join(s_list)
        seqs.append(SeqRecord(Seq(s_str), id=f"seq_{i}", description="rc_test"))
    
    out_file = "rc_test.fa"
    SeqIO.write(seqs, out_file, "fasta")
    print(f"Generated {out_file} with motif planted on REVERSE strand only.")
    print(f"Expected discovery: Should find '{TRUE_MOTIF}' (or its RC) if both-strand works.")

if __name__ == "__main__":
    generate()
