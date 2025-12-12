import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def generate():
    # 2 Motifs: 
    # M1 (Len 8): "AAAAAATT"
    # M2 (Len 12): "GGGGGGCCCCCC"
    
    seqs = []
    for i in range(50):
        # random background length 100
        bases = ['A','C','G','T']
        s = "".join(random.choice(bases) for _ in range(100))
        
        # Plant M1
        pos1 = random.randint(0, 40)
        s = s[:pos1] + "AAAAAATT" + s[pos1+8:]
        
        # Plant M2
        pos2 = random.randint(50, 90)
        s = s[:pos2] + "GGGGGGCCCCCC" + s[pos2+12:]
        
        seqs.append(SeqRecord(Seq(s), id=f"seq_{i}", description="synthetic"))
        
    SeqIO.write(seqs, "synthetic_discover_test.fa", "fasta")
    print("Generated synthetic_discover_test.fa with 2 planted motifs.")

if __name__ == "__main__":
    generate()
