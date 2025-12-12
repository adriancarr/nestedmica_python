#!/usr/bin/env python3
"""
Compare Python and Java recovered motifs against the 10 planted ground truth motifs.
"""

import xml.etree.ElementTree as ET
import numpy as np

# Ground truth motifs from generate_hard_synthetic.py (seed=42)
GROUND_TRUTH_CONSENSUS = [
    "GTAGGTAAGC",  # Motif 0 (strength=0.95)
    "GGGGTATTTG",  # Motif 1 (strength=0.89)
    "CACTTCCCTT",  # Motif 2 (strength=0.83)
    "AATCCATAAG",  # Motif 3 (strength=0.77)
    "GGCTTTTGCC",  # Motif 4 (strength=0.71)
    "GCGTGTTAGA",  # Motif 5 (strength=0.64)
    "GGAAGCTATC",  # Motif 6 (strength=0.58)
    "CCACACTTGT",  # Motif 7 (strength=0.52)
    "GTATGGCATC",  # Motif 8 (strength=0.46)
    "TTCCCCCTCA",  # Motif 9 (strength=0.40)
]

STRENGTHS = [0.95, 0.89, 0.83, 0.77, 0.71, 0.64, 0.58, 0.52, 0.46, 0.40]

BASES = ['A', 'C', 'G', 'T']

def consensus_to_pwm(consensus, strength=0.9):
    """Convert consensus to approximate PWM."""
    pwm = np.zeros((len(consensus), 4))
    for i, base in enumerate(consensus):
        base_idx = BASES.index(base)
        for j in range(4):
            if j == base_idx:
                pwm[i, j] = 0.25 + 0.75 * strength
            else:
                pwm[i, j] = (1.0 - (0.25 + 0.75 * strength)) / 3.0
    return pwm

def parse_python_xms(filename):
    """Parse Python output format."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    motifs = []
    for motif in root.findall('motif'):
        wm = motif.find('weightMatrix')
        rows = []
        for line in wm.text.strip().split('\n'):
            rows.append([float(x) for x in line.split()])
        motifs.append(np.array(rows))
    return motifs

def parse_java_xms(filename):
    """Parse Java output format with namespace."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    ns = {'ns': 'http://biotiffin.org/XMS/'}
    
    motifs = []
    for motif in root.findall('ns:motif', ns):
        wm = motif.find('ns:weightmatrix', ns)
        cols = []
        for col in sorted(wm.findall('ns:column', ns), key=lambda x: int(x.attrib['pos'])):
            weights = {}
            for w in col.findall('ns:weight', ns):
                sym = w.attrib['symbol']
                val = float(w.text)
                weights[sym] = val
            cols.append([
                weights.get('adenine', 0),
                weights.get('cytosine', 0),
                weights.get('guanine', 0),
                weights.get('thymine', 0)
            ])
        motifs.append(np.array(cols))
    return motifs

def consensus(pwm):
    """Get consensus sequence from PWM."""
    return ''.join(BASES[np.argmax(row)] for row in pwm)

def pwm_correlation(pwm1, pwm2):
    """Compute mean Pearson correlation between two PWMs."""
    if len(pwm1) != len(pwm2):
        return -1
    correlations = []
    for i in range(len(pwm1)):
        corr = np.corrcoef(pwm1[i], pwm2[i])[0, 1]
        if not np.isnan(corr):
            correlations.append(corr)
    return np.mean(correlations) if correlations else 0

def information_content(pwm):
    """Calculate mean IC per position."""
    ic = []
    for row in pwm:
        row = np.array(row)
        row = row / row.sum()
        row = np.clip(row, 1e-10, 1)
        ic.append(2 + np.sum(row * np.log2(row)))
    return np.mean(ic)

def find_best_match(recovered_motifs, ground_truth_pwm):
    """Find which recovered motif best matches the ground truth."""
    best_idx = -1
    best_corr = -1
    for i, rec in enumerate(recovered_motifs):
        if len(rec) == len(ground_truth_pwm):
            corr = pwm_correlation(rec, ground_truth_pwm)
            if corr > best_corr:
                best_corr = corr
                best_idx = i
    return best_idx, best_corr

def main():
    print("=" * 80)
    print("CHALLENGING SYNTHETIC TEST - MOTIF RECOVERY COMPARISON")
    print("=" * 80)
    
    # Load recovered motifs
    try:
        py_motifs = parse_python_xms("hard_python_50.xms")
    except Exception as e:
        print(f"Error loading Python output: {e}")
        py_motifs = []
    
    try:
        java_motifs = parse_java_xms("hard_java_50.xms")
    except Exception as e:
        print(f"Error loading Java output: {e}")
        java_motifs = []
    
    print(f"\nLoaded {len(py_motifs)} Python motifs and {len(java_motifs)} Java motifs")
    
    print("\n" + "=" * 80)
    print("GROUND TRUTH vs RECOVERED MOTIFS")
    print("=" * 80)
    
    print(f"\n{'Ground Truth':<12} {'Strength':<10} {'Python Match':<15} {'Py Corr':<10} {'Java Match':<15} {'Ja Corr':<10}")
    print("-" * 80)
    
    py_matches = []
    ja_matches = []
    
    for gt_idx, (gt_cons, strength) in enumerate(zip(GROUND_TRUTH_CONSENSUS, STRENGTHS)):
        gt_pwm = consensus_to_pwm(gt_cons, strength)
        
        # Find best Python match
        if py_motifs:
            py_best_idx, py_corr = find_best_match(py_motifs, gt_pwm)
            if py_best_idx >= 0:
                py_cons = consensus(py_motifs[py_best_idx])
            else:
                py_cons = "---"
                py_corr = 0
        else:
            py_cons = "N/A"
            py_corr = 0
        
        # Find best Java match
        if java_motifs:
            ja_best_idx, ja_corr = find_best_match(java_motifs, gt_pwm)
            if ja_best_idx >= 0:
                ja_cons = consensus(java_motifs[ja_best_idx])
            else:
                ja_cons = "---"
                ja_corr = 0
        else:
            ja_cons = "N/A"
            ja_corr = 0
        
        py_matches.append(py_corr)
        ja_matches.append(ja_corr)
        
        print(f"Motif {gt_idx} ({gt_cons[:6]}...) {strength:<10.2f} {py_cons:<15} {py_corr:<10.3f} {ja_cons:<15} {ja_corr:<10.3f}")
    
    print("-" * 80)
    print(f"{'MEAN':<23} {'':<10} {'':<15} {np.mean(py_matches):<10.3f} {'':<15} {np.mean(ja_matches):<10.3f}")
    
    print("\n" + "=" * 80)
    print("RECOVERED MOTIF DETAILS")
    print("=" * 80)
    
    print("\n### Python Recovered Motifs ###")
    for i, pwm in enumerate(py_motifs):
        ic = information_content(pwm)
        print(f"  Motif {i}: {consensus(pwm)} (IC={ic:.2f} bits)")
    
    print("\n### Java Recovered Motifs ###")
    for i, pwm in enumerate(java_motifs):
        ic = information_content(pwm)
        print(f"  Motif {i}: {consensus(pwm)} (IC={ic:.2f} bits)")

if __name__ == "__main__":
    main()
