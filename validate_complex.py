import re
import sys
import numpy as np
import xml.etree.ElementTree as ET

def parse_truth(filename):
    motifs = {}
    with open(filename, 'r') as f:
        for line in f:
            # Motif 0 (L=6): CTCATG
            m = re.search(r'Motif (\d+) \(L=(\d+)\): ([ACGT]+)', line)
            if m:
                idx = int(m.group(1))
                seq = m.group(3)
                motifs[idx] = seq
    return motifs

def parse_xms(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    motifs = []
    
    for wm in root.iter():
        if 'weightmatrix' in wm.tag.lower():
            text = wm.text
            if text and text.strip():
                # Simplified format (space separated values)
                values = [float(x) for x in text.split()]
                # Handle attributes with potential namespaces or just get length
                length_attr = None
                for k, v in wm.attrib.items():
                    if 'length' in k:
                        length_attr = int(v)
                        break
                
                if length_attr is None:
                    length_attr = len(values)//4
                    
                length = length_attr
                if len(values) == length * 4:
                    matrix = np.array(values).reshape((length, 4))
                    motifs.append(matrix)
            else:
                # Verbose XML format (Java output)
                columns = []
                # Handle namespaces in child tags too
                for col in wm:
                    if 'column' in col.tag:
                        weights = {'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0}
                        for w in col:
                            if 'weight' in w.tag:
                                sym = None
                                for k, v in w.attrib.items():
                                    if 'symbol' in k:
                                        sym = v
                                        break
                                if sym:
                                    val = float(w.text)
                                    if 'adenine' in sym: weights['A'] = val
                                    elif 'cytosine' in sym: weights['C'] = val
                                    elif 'guanine' in sym: weights['G'] = val
                                    elif 'thymine' in sym: weights['T'] = val
                        columns.append([weights['A'], weights['C'], weights['G'], weights['T']])
                
                if columns:
                    matrix = np.array(columns)
                    motifs.append(matrix)
    return motifs

def get_consensus(pwm):
    bases = "ACGT"
    seq = ""
    for row in pwm:
        idx = np.argmax(row)
        seq += bases[idx]
    return seq

def score_match(pwm, target_seq):
    """Align target sequence to PWM and return score (log prob)."""
    # Simple best-hit alignment
    target_indices = [{'A':0,'C':1,'G':2,'T':3}[c] for c in target_seq]
    
    best_score = -np.inf
    
    # PWM length
    PL = len(pwm)
    # Target length
    TL = len(target_seq)
    
    # Try all overlaps?
    # Actually, we want to know if the PWM represents the target.
    # If standard alignment.
    # Or just check consensus similarity.
    # Let's use Consensus Sequence Similarity (Levenshtein) or Substring check.
    
    pwm_cons = get_consensus(pwm)
    
    # Check if target is substring of pwm_cons or vice versa
    if target_seq in pwm_cons: return 100 + len(target_seq)
    if pwm_cons in target_seq: return 100 + len(pwm_cons)
    
    # Basic edit distance?
    import difflib
    ratio = difflib.SequenceMatcher(None, target_seq, pwm_cons).ratio()
    return ratio

def main():
    if len(sys.argv) < 3:
        print("Usage: validate_complex.py <xms_file> <truth_file>")
        sys.exit(1)
        
    xms_motifs = parse_xms(sys.argv[1])
    truth_motifs = parse_truth(sys.argv[2])
    
    print(f"Found {len(xms_motifs)} motifs in XMS.")
    print(f"Loaded {len(truth_motifs)} ground truth motifs.")
    
    matches = {}
    used_xms = set()
    
    # For each truth, find best XMS match
    correct_count = 0
    
    print("\nMatching Results:")
    print(f"{'Truth ID':<10} {'Truth Seq':<20} {'Length':<8} | {'Best Match':<20} {'Length':<8} {'Score':<8}")
    print("-" * 80)
    
    for idx in sorted(truth_motifs.keys()):
        truth_seq = truth_motifs[idx]
        best_match_idx = -1
        best_score = -1.0
        best_cons = ""
        
        for i, m in enumerate(xms_motifs):
            cons = get_consensus(m)
            score = score_match(m, truth_seq)
            if score > best_score:
                best_score = score
                best_match_idx = i
                best_cons = cons
        
        len_match = len(best_cons)
        print(f"{idx:<10} {truth_seq:<20} {len(truth_seq):<8} | {best_cons:<20} {len_match:<8} {best_score:.2f}")
        
        # Heuristic for success: Exact match or close substring
        if best_score > 0.8: # high ratio
             correct_count += 1
             
    print(f"\nResult: {correct_count}/{len(truth_motifs)} motifs recovered correctly.")

if __name__ == "__main__":
    main()
