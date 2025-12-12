import numpy as np
import random

def create_pwm(length, bias_char=None):
    """Create a sharp PWM of given length."""
    pwm = np.zeros((length, 4))
    chars = "ACGT"
    
    for i in range(length):
        # Pick a dominant base
        if bias_char and random.random() < 0.3:
             dom = "ACGT".index(bias_char)
        else:
             dom = random.randint(0, 3)
             
        # Assign high probability (0.8-0.95) to dominant
        prob = random.uniform(0.8, 0.95)
        pwm[i, dom] = prob
        
        # Distribute rest
        remain = 1.0 - prob
        others = [j for j in range(4) if j != dom]
        for j in others:
            pwm[i, j] = remain / 3.0
            
    return pwm

def sample_from_pwm(pwm):
    """Sample a string from a PWM."""
    seq = []
    chars = "ACGT"
    for row in pwm:
        # row sums to 1.0 (approx)
        row = row / row.sum()
        char = np.random.choice(list(chars), p=row)
        seq.append(char)
    return "".join(seq)

def get_consensus(pwm):
    chars = "ACGT"
    return "".join([chars[np.argmax(row)] for row in pwm])

def generate_complex_dataset(filename, truth_filename, n_seqs=5000, seq_len=100):
    print(f"Generating {n_seqs} sequences (100bp) with complex motif model...")
    
    # Define 5 ground truth motifs with different lengths
    # sizes: 6, 8, 10, 12, 14
    motif_specs = [
        (6, None),
        (8, 'A'), # A-rich
        (10, 'G'), # G-rich
        (12, None), 
        (14, 'T') # T-rich
    ]
    
    motifs = []
    for l, bias in motif_specs:
        motifs.append(create_pwm(l, bias))
        
    # Write truth to file
    with open(truth_filename, 'w') as f:
        for i, m in enumerate(motifs):
            cons = get_consensus(m)
            f.write(f"Motif {i} (L={len(m)}): {cons}\n")
            print(f"  Motif {i}: {cons}")
            
    with open(filename, 'w') as f:
        for i in range(n_seqs):
            # Background
            seq_bg = "".join(random.choices("ACGT", k=seq_len))
            seq_list = list(seq_bg)
            
            # Number of motifs to insert: 0-4
            n_inserts = random.randint(0, 4)
            
            occupied = [] # Keep track of positions to avoid overlap
            
            for _ in range(n_inserts):
                # Pick a random motif type
                m_idx = random.randint(0, len(motifs) - 1)
                pwm = motifs[m_idx]
                motif_len = len(pwm)
                
                # Sample instance
                instance = sample_from_pwm(pwm)
                
                # Find position
                # Try 10 times to find non-overlapping spot
                for _try in range(10):
                    pos = random.randint(0, seq_len - motif_len)
                    # Check overlap
                    overlap = False
                    for (occ_start, occ_end) in occupied:
                        if not (pos + motif_len <= occ_start or pos >= occ_end):
                            overlap = True
                            break
                    
                    if not overlap:
                        # Insert
                        for k, char in enumerate(instance):
                            seq_list[pos+k] = char
                        occupied.append((pos, pos + motif_len))
                        break
            
            final_seq = "".join(seq_list)
            f.write(f">seq_{i} label=1\n{final_seq}\n")
            
    print(f"Done. Saved to {filename} and {truth_filename}.")

if __name__ == "__main__":
    generate_complex_dataset("complex_test.fa", "complex_truth.txt", n_seqs=10000)
