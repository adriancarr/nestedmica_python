import numpy as np
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def create_pwm(consensus, conservation):
    """Create a PWM from a consensus string with given conservation (0-1)."""
    length = len(consensus)
    pwm = np.zeros((4, length))
    bases = ['A', 'C', 'G', 'T']
    base_map = {b: i for i, b in enumerate(bases)}
    
    for i, char in enumerate(consensus):
        target = base_map[char]
        # Distribute probability
        # Conserved base gets: conservation + (1-conservation)/4
        # Others get: (1-conservation)/4
        noise = (1.0 - conservation) / 4.0
        prob = conservation + noise
        
        for b_idx in range(4):
            if b_idx == target:
                pwm[b_idx, i] = prob
            else:
                pwm[b_idx, i] = noise
                
    # Normalize columns to sum to 1
    pwm = pwm / pwm.sum(axis=0, keepdims=True)
    return pwm

def sample_motif(pwm):
    """Sample a sequence string from a PWM."""
    bases = ['A', 'C', 'G', 'T']
    length = pwm.shape[1]
    res = []
    for i in range(length):
        probs = pwm[:, i]
        choice = np.random.choice(bases, p=probs)
        res.append(choice)
    return "".join(res)

def generate():
    random.seed(42)
    np.random.seed(42)
    
    # Ground Truth Definitions
    # M1: Strong, Length 8 (Easy)
    pwm1 = create_pwm("ACGTACGT", 0.95)
    # M2: Medium, Length 10 (Moderate)
    pwm2 = create_pwm("TTTGGCCAAA", 0.8)
    # M3: Weak, Length 14 (Hard - close to background noise)
    pwm3 = create_pwm("TATATATATATATA", 0.6)
    
    seqs = []
    # 500 Sequences
    for i in range(500):
        # Background: Random uniform
        L = 200
        bases = ['A', 'C', 'G', 'T']
        s_list = [random.choice(bases) for _ in range(L)]
        
        # Plant M1 (100% chance for strong signal)
        if random.random() < 1.0:
            m = sample_motif(pwm1)
            pos = random.randint(0, 50)
            s_list[pos:pos+len(m)] = list(m)
            
        # Plant M2 (80% chance)
        if random.random() < 0.8:
            m = sample_motif(pwm2)
            pos = random.randint(60, 120)
            s_list[pos:pos+len(m)] = list(m)

        # Plant M3 (70% chance)
        if random.random() < 0.7:
            m = sample_motif(pwm3)
            pos = random.randint(130, 180)
            s_list[pos:pos+len(m)] = list(m)
            
        s_str = "".join(s_list)
        seqs.append(SeqRecord(Seq(s_str), id=f"seq_{i}", description="complex_synthetic"))
        
    out_file = "complex_varying_test.fa"
    SeqIO.write(seqs, out_file, "fasta")
    print(f"Generated {out_file} with 3 motifs (Strong, Medium, Weak).")

if __name__ == "__main__":
    generate()
