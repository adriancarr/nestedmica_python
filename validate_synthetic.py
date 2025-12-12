#!/usr/bin/env python3
"""
Compare recovered motifs against ground truth planted PWM.
"""

import xml.etree.ElementTree as ET
import numpy as np

# Ground truth PWM (A, C, G, T order)
GROUND_TRUTH = np.array([
    [0.90, 0.03, 0.03, 0.04],  # pos 0: A-rich
    [0.03, 0.90, 0.04, 0.03],  # pos 1: C-rich
    [0.03, 0.04, 0.90, 0.03],  # pos 2: G-rich
    [0.04, 0.03, 0.03, 0.90],  # pos 3: T-rich
    [0.90, 0.03, 0.03, 0.04],  # pos 4: A-rich
    [0.03, 0.90, 0.04, 0.03],  # pos 5: C-rich
    [0.03, 0.04, 0.90, 0.03],  # pos 6: G-rich
    [0.04, 0.03, 0.03, 0.90],  # pos 7: T-rich
    [0.50, 0.50, 0.00, 0.00],  # pos 8: A or C
    [0.00, 0.00, 0.50, 0.50],  # pos 9: G or T
])

def parse_python_xms(filename):
    """Parse Python output format."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    for motif in root.findall('motif'):
        wm = motif.find('weightMatrix')
        rows = []
        for line in wm.text.strip().split('\n'):
            rows.append([float(x) for x in line.split()])
        return np.array(rows)
    return None

def parse_java_xms(filename):
    """Parse Java output format with namespace."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    ns = {'ns': 'http://biotiffin.org/XMS/'}
    
    for motif in root.findall('ns:motif', ns):
        wm = motif.find('ns:weightmatrix', ns)
        cols = []
        for col in sorted(wm.findall('ns:column', ns), key=lambda x: int(x.attrib['pos'])):
            weights = {}
            for w in col.findall('ns:weight', ns):
                sym = w.attrib['symbol']
                val = float(w.text)
                weights[sym] = val
            # Order: A, C, G, T
            cols.append([
                weights.get('adenine', 0),
                weights.get('cytosine', 0),
                weights.get('guanine', 0),
                weights.get('thymine', 0)
            ])
        return np.array(cols)
    return None

def consensus(pwm):
    """Get consensus sequence from PWM."""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(bases[np.argmax(row)] for row in pwm)

def position_correlation(recovered, truth):
    """Compute per-position Pearson correlation."""
    correlations = []
    for i in range(min(len(recovered), len(truth))):
        corr = np.corrcoef(recovered[i], truth[i])[0, 1]
        correlations.append(corr)
    return correlations

def information_content(pwm):
    """Calculate per-position information content."""
    ic = []
    for row in pwm:
        row = np.array(row)
        row = row / row.sum()  # normalize
        row = np.clip(row, 1e-10, 1)  # avoid log(0)
        ic.append(2 + np.sum(row * np.log2(row)))
    return ic

def try_reverse_complement(recovered, truth):
    """Check if reverse complement matches better."""
    # Reverse the positions and swap A<->T, C<->G
    rc = recovered[::-1, [3, 2, 1, 0]]  # Reverse rows and swap columns
    
    fwd_corr = np.mean(position_correlation(recovered, truth))
    rc_corr = np.mean(position_correlation(rc, truth))
    
    return fwd_corr, rc_corr, rc_corr > fwd_corr

def main():
    print("=" * 60)
    print("SYNTHETIC DATA VALIDATION")
    print("=" * 60)
    
    print("\n### Ground Truth PWM ###")
    print(f"Consensus: {consensus(GROUND_TRUTH)}")
    print(f"Length: {len(GROUND_TRUTH)}")
    
    # Parse outputs - using v3 for Python
    python_pwm = parse_python_xms("synthetic_python_v3.xms")
    java_pwm = parse_java_xms("synthetic_java.xms")
    
    print("\n### Python Recovered ###")
    print(f"Consensus: {consensus(python_pwm)}")
    py_fwd, py_rc, py_is_rc = try_reverse_complement(python_pwm, GROUND_TRUTH)
    print(f"Mean Correlation (forward): {py_fwd:.3f}")
    print(f"Mean Correlation (rev-comp): {py_rc:.3f}")
    py_ic = np.mean(information_content(python_pwm))
    print(f"Mean IC: {py_ic:.2f} bits")
    
    print("\n### Java Recovered ###")
    print(f"Consensus: {consensus(java_pwm)}")
    ja_fwd, ja_rc, ja_is_rc = try_reverse_complement(java_pwm, GROUND_TRUTH)
    print(f"Mean Correlation (forward): {ja_fwd:.3f}")
    print(f"Mean Correlation (rev-comp): {ja_rc:.3f}")
    ja_ic = np.mean(information_content(java_pwm))
    print(f"Mean IC: {ja_ic:.2f} bits")
    
    print("\n### Position-by-Position Analysis ###")
    print(f"{'Pos':<4} {'Truth':<6} {'Python':<8} {'Java':<8} {'Py Corr':<10} {'Ja Corr':<10}")
    print("-" * 56)
    
    py_corrs = position_correlation(python_pwm, GROUND_TRUTH)
    ja_corrs = position_correlation(java_pwm, GROUND_TRUTH)
    
    for i in range(10):
        truth_cons = ['A', 'C', 'G', 'T'][np.argmax(GROUND_TRUTH[i])]
        py_cons = ['A', 'C', 'G', 'T'][np.argmax(python_pwm[i])]
        ja_cons = ['A', 'C', 'G', 'T'][np.argmax(java_pwm[i])]
        
        print(f"{i:<4} {truth_cons:<6} {py_cons:<8} {ja_cons:<8} {py_corrs[i]:<10.3f} {ja_corrs[i]:<10.3f}")
    
    print("\n### Summary ###")
    print(f"Ground Truth: {consensus(GROUND_TRUTH)}")
    print(f"Python:       {consensus(python_pwm)} (Corr: {np.mean(py_corrs):.3f})")
    print(f"Java:         {consensus(java_pwm)} (Corr: {np.mean(ja_corrs):.3f})")

if __name__ == "__main__":
    main()
