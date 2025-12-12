
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools

def load_motifs(filename):
    """Parses XMS file and returns a list of PWMs (numpy arrays)."""
    tree = ET.parse(filename)
    root = tree.getroot()
    # Handle optional namespace in root tag
    ns = {}
    if '}' in root.tag:
        ns = {'ns': root.tag.split('}')[0].strip('{')}
    
    motifs = []
    
    # XMS format structure varies slightly between versions (simple vs full)
    # Simple (Python port): <motif><weightMatrix>...
    # Full (Java): <motif><weightmatrix>... (properties and columns)
    
    # Try finding 'motif' tags
    findall_path = 'ns:motif' if ns else 'motif'
    
    for motif_node in root.findall(findall_path, namespaces=ns):
        wm_node = motif_node.find('ns:weightMatrix' if ns else 'weightMatrix', namespaces=ns)
        if wm_node is None:
             wm_node = motif_node.find('ns:weightmatrix' if ns else 'weightmatrix', namespaces=ns) # Java lower case 'weightmatrix'
        
        if wm_node is None:
            continue
            
        # Parse matrix
        # Python simple format: text content
        if wm_node.text and wm_node.text.strip():
            rows = []
            for line in wm_node.text.strip().split('\n'):
                rows.append([float(x) for x in line.split()])
            motifs.append(np.array(rows))
        else:
            # Java full format: <column><weight>...
            cols = []
            col_nodes = wm_node.findall('ns:column' if ns else 'column', namespaces=ns)
            # Sort by pos just in case
            col_nodes.sort(key=lambda x: int(x.attrib['pos']))
            
            for col in col_nodes:
                weights = {}
                for w in col.findall('ns:weight' if ns else 'weight', namespaces=ns):
                    sym = w.attrib['symbol']
                    val = float(w.text)
                    weights[sym] = val
                # Order: A, C, G, T? Or match Python alphabet?
                # Java alphabet is defined in file usually, or standard DNA
                # Let's assume standard DNA: adenine, cytosine, guanine, thymine
                # Python port uses Biopython dict or similar.
                # Standard order usually A, C, G, T.
                ordered_weights = [
                    weights.get('adenine', 0.0),
                    weights.get('cytosine', 0.0),
                    weights.get('guanine', 0.0),
                    weights.get('thymine', 0.0)
                ]
                cols.append(ordered_weights)
            motifs.append(np.array(cols)) # Java makes it cols x 4, Python might be rows x 4?
            
    return motifs

def calculate_ic(pwm):
    """Calculates Information Content (bits) for a PWM."""
    # PWM is assumed to be Position x 4 (normalized probabilities)
    # IC = 2 + sum(p * log2(p))
    pwm = np.array(pwm)
    # Normalize if needed
    row_sums = pwm.sum(axis=1, keepdims=True)
    pwm_norm = pwm / row_sums
    
    ic_per_col = 2 + np.sum(pwm_norm * np.nan_to_num(np.log2(pwm_norm)), axis=1)
    return ic_per_col

def kullback_leibler(p, q):
    """KL Divergence between two distributions P and Q."""
    p = np.asarray(p)
    q = np.asarray(q)
    return np.sum(np.where(p != 0, p * np.log2(p / q), 0))

def compare_motif_pair(m1, m2):
    """Compares two motifs. Returns avg euclidean distance and avg KL divergence."""
    # Ensure same length
    l = min(len(m1), len(m2))
    m1 = m1[:l]
    m2 = m2[:l]
    
    # Euclidean
    dist = np.linalg.norm(m1 - m2)
    
    # Pearson Correlation (per column average)
    corrs = []
    for c in range(l):
        corrs.append(np.corrcoef(m1[c], m2[c])[0,1])
    
    avg_corr = np.mean(corrs)
    
    return dist, avg_corr

def plot_logos(python_motifs, java_motifs, filename="comparison_plot.png"):
    n_py = len(python_motifs)
    n_ja = len(java_motifs)
    
    fig, axes = plt.subplots(max(n_py, n_ja), 2, figsize=(12, 3 * max(n_py, n_ja)))
    
    # Helper to plot simple bar chart "logo"
    def plot_bar_logo(ax, pwm, title):
        pwm = np.array(pwm)
        x = range(len(pwm))
        bottom = np.zeros(len(pwm))
        colors = ['green', 'blue', 'orange', 'red'] # A C G T
        bases = ['A', 'C', 'G', 'T']
        
        # Plot stacked bars for probabilities
        for i in range(4):
            ax.bar(x, pwm[:, i], bottom=bottom, color=colors[i], label=bases[i])
            bottom += pwm[:, i]
            
        ax.set_ylim(0, 1)
        ax.set_title(title)
        ax.set_xticks(x)
        if hasattr(ax, 'legend'):
           pass # ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1))

    # Match motifs? A simple heuristic matching
    # Compute similarity matrix
    sim_matrix = np.zeros((n_py, n_ja))
    for i in range(n_py):
        for j in range(n_ja):
            _, corr = compare_motif_pair(python_motifs[i], java_motifs[j])
            sim_matrix[i, j] = corr
            
    # Simple greedy matching
    matches = [] # (py_idx, java_idx)
    used_java = set()
    for i in range(n_py):
        best_j = -1
        best_corr = -1
        for j in range(n_ja):
            if j not in used_java and sim_matrix[i, j] > best_corr:
                best_corr = sim_matrix[i, j]
                best_j = j
        if best_j != -1:
            matches.append((i, best_j))
            used_java.add(best_j)
            
    # Plot matches
    for idx, (py_idx, ja_idx) in enumerate(matches):
        plot_bar_logo(axes[idx, 0], python_motifs[py_idx], f"Python Motif {py_idx}")
        plot_bar_logo(axes[idx, 1], java_motifs[ja_idx], f"Java Motif {ja_idx} (Corr: {sim_matrix[py_idx, ja_idx]:.2f})")

    # Add legend to the last plot
    axes[0,0].legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Plot saved to {filename}")

if __name__ == "__main__":
    py_motifs = load_motifs("output_v3.xms")
    java_motifs = load_motifs("java_output.xms")
    
    print(f"Loaded {len(py_motifs)} Python motifs and {len(java_motifs)} Java motifs.")
    
    # Calculate and print metrics
    print("\n--- Metrics ---")
    
    # Check alphabet mapping for Python output.
    # Python output seems to be just numbers. 
    # Mocca.java (Python port) likely outputs in A, C, G, T order if it follows BioJava defaults,
    # but let's verify if 'nestedmica.model.motif.NMWeightMatrix' enforces this.
    # Looking at the code text from previous steps, Mocca.py output format:
    # 0.2942 0.2717 0.3460 0.0882  <- this is likely A C G T
    
    # Let's perform matching and print quantitative stats
    sim_matrix = np.zeros((len(py_motifs), len(java_motifs)))
    dist_matrix = np.zeros((len(py_motifs), len(java_motifs)))
    
    for i, pm in enumerate(py_motifs):
        for j, jm in enumerate(java_motifs):
            d, c = compare_motif_pair(pm, jm)
            sim_matrix[i, j] = c
            dist_matrix[i, j] = d
            
    # Print best matches
    for i in range(len(py_motifs)):
        best_j = np.argmax(sim_matrix[i])
        corr = sim_matrix[i, best_j]
        
        pm = py_motifs[i]
        jm = java_motifs[best_j]
        
        ic_py = np.mean(calculate_ic(pm))
        ic_java = np.mean(calculate_ic(jm))
        
        print(f"Python Motif {i} matches Java Motif {best_j}")
        print(f"  Pearson Correlation: {corr:.3f}")
        print(f"  Euclidean Distance : {dist_matrix[i, best_j]:.3f}")
        print(f"  Mean IC (Bits)     : Py={ic_py:.2f}, Java={ic_java:.2f}") 
        print(f"  Delta IC           : {ic_java - ic_py:.2f}")
        pass
        
    plot_logos(py_motifs, java_motifs)
