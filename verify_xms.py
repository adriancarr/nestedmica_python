import sys
import xml.etree.ElementTree as ET
import numpy as np

def parse_xms(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    
    motifs = []
    
    for wm in root.iter('weightMatrix'):
        text = wm.text.strip()
        # Split by whitespace
        values = [float(x) for x in text.split()]
        # Reshape (L, 4)
        # Length is in attrib
        length = int(wm.attrib.get('length', len(values)//4))
        
        matrix = np.array(values).reshape((length, 4))
        motifs.append(matrix)
        
    return motifs

def get_consensus(pwm):
    bases = "ACGT"
    seq = ""
    for row in pwm: # shape (L, 4)
        idx = np.argmax(row)
        seq += bases[idx]
    return seq

def main():
    if len(sys.argv) < 3:
        print("Usage: verify_xms.py <xms_file> <expected_seq>")
        sys.exit(1)
        
    xms_file = sys.argv[1]
    expected = sys.argv[2]
    
    try:
        motifs = parse_xms(xms_file)
    except Exception as e:
        print(f"Error parsing XMS: {e}")
        sys.exit(1)
        
    found_match = False
    print(f"Found {len(motifs)} motifs.")
    for i, m in enumerate(motifs):
        cons = get_consensus(m)
        print(f"Motif {i}: {cons}")
        if cons == expected:
            found_match = True
            
    if found_match:
        print("SUCCESS: Found expected motif.")
    else:
        print("FAILURE: Did not find expected motif.")
        sys.exit(1)

if __name__ == "__main__":
    main()
