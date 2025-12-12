# Release v1.2.0: Higher-Order Backgrounds & Reverse Complement

This release adds two biologically important features for improved motif specificity.

## New Features 🚀

### 1. Higher-Order Markov Background (default order-3)
- **Learns** dinucleotide/trinucleotide patterns from input sequences
- **Reduces** false positives in biased sequences (GC-rich, AT-rich)
- **CLI**: `-bgOrder 3` (configurable 0-5)

### 2. Reverse Complement Scanning (default ON)
- **Scans** both DNA strands (forward + reverse complement)
- **Essential** for TF binding sites which are strand-independent
- **CLI**: Enabled by default via `both_strands=True`

## Usage
```bash
# Default: order-3 background, both strands
python3 -m nestedmica.apps.mocca_fast -seqs data.fa -numMotifs 2 -out result.xms

# Custom background order
python3 -m nestedmica.apps.discover -seqs data.fa -out discovered.xms -bgOrder 2
```

---

# Release v1.1.0: Auto-Discovery & Bayesian Evidence

This major release introduces automated motif discovery and rigorous statistical model selection.

## New Features 🚀

### 1. Automated Motif Discovery
- **Tool**: `python3 -m nestedmica.apps.discover`
- **Functionality**: Automatically determines the optimal number of motifs ($N$) and their lengths.
- **Method**: Sequential hypothesis testing using Bayes Factors derived from Global Log Evidence.

### 2. Bayesian Model Selection
- **Global Log Evidence ($\ln Z$)**: Full calculation of the model evidence integral.
- **Interpretation**: Allows direct comparison of models with different complexities (e.g., 3 motifs vs 4 motifs).
- **Stopping Criterion**: The algorithm stops adding motifs when $\Delta \ln Z$ drops or becomes minimal.

### 3. Smart Ensemble Sizing
- **Heuristic**: `EnsembleSize = 50 * log10(TotalBases)`
- **Benefit**: Automatically scales computational effort to the dataset size, preventing under-sampling on large genomes.

## Improvements 🛠️
- **Checkpointing**: Full support for saving/restoring Evidence state.
- **Performance**: Threaded batch processing remains ~3x faster than legacy Java.
- **Validation**: Verified on complex synthetic datasets with variable signal strength.

