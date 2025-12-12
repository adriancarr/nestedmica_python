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

## Usage
```bash
# Auto-discover optimal motifs
python3 -m nestedmica.apps.discover -seqs data.fa -out discovered.xms

# Inspect Evidence
# Output at end of run:
# "GlobalLogEvidence: -195937.69"
```
