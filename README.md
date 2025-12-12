# NestedMICA Python

A high-performance Python implementation of the Nested MICA motif discovery algorithm, featuring Cython optimization, parallel processing, and rigorous Bayesian model selection.

**Version**: 1.2.0

## Features

| Feature | Description |
|---------|-------------|
| **Auto-Discovery** | Automatically determines optimal number of motifs using Bayesian Evidence |
| **Both-Strand Scanning** | Scans forward and reverse complement for TF binding sites |
| **Higher-Order Background** | 3rd-order Markov model reduces false positives in biased sequences |
| **Adaptive MCMC** *(experimental)* | Dynamically adjusts proposal distribution for faster convergence |
| **Variable Length Motifs** | Insert/Delete moves for dynamic motif sizing |
| **Parallel Processing** | Multi-threaded Cython for ~3x speedup over Java |
| **Checkpointing** | Save/restore training state for long-running jobs |

---

## Installation

Requires Python 3.9+ and a C compiler.

```bash
pip install numpy biopython cython
python3 setup.py build_ext --inplace
```

---

## Quick Start

### Auto-Discovery (Recommended)

Let the algorithm find the optimal number and length of motifs:

```bash
python3 -m nestedmica.apps.discover -seqs data.fa -out discovered.xms -threads 8
```

### Manual Run (Advanced)

Find exactly N motifs of a specific length:

```bash
python3 -m nestedmica.apps.mocca_fast \
  -seqs data.fa \
  -numMotifs 3 \
  -motifLength 12 \
  -threads 8 \
  -out result.xms
```

---

## CLI Options

### `discover` (Auto-Discovery)

| Option | Default | Description |
|--------|---------|-------------|
| `-seqs` | *required* | Input FASTA file |
| `-out` | *required* | Output XMS file |
| `-minMotifs` | 1 | Minimum N to test |
| `-maxMotifs` | 10 | Maximum N to test |
| `-bayesThreshold` | 5.0 | Log Bayes Factor threshold to accept new motif |
| `-ensembleSize` | *auto* | Particle count (auto: 50×log₁₀(bases)) |
| `-threads` | -1 | Number of threads (-1 = all cores) |
| `-bgOrder` | 3 | Background Markov order (0-5) |
| `-adaptiveMCMC` | True | Use adaptive proposal distribution |
| `-noAdaptiveMCMC` | — | Disable adaptive proposals |

### `mocca_fast` (Manual Run)

| Option | Default | Description |
|--------|---------|-------------|
| `-seqs` | *required* | Input FASTA file |
| `-out` | *required* | Output XMS file |
| `-numMotifs` | 1 | Number of motifs to find |
| `-motifLength` | 10 | Initial motif length |
| `-maxCycles` | 500 | Maximum sampling cycles |
| `-ensembleSize` | 20 | Particle count |
| `-stopIqr` | 0.01 | Convergence threshold (IQR of ensemble) |
| `-threads` | -1 | Number of threads (-1 = all cores) |
| `-bgOrder` | 3 | Background Markov order (0-5) |
| `-adaptiveMCMC` | True | Use adaptive proposal distribution |
| `-checkpoint` | — | File to save checkpoints |
| `-restart` | — | Checkpoint to resume from |

---

## Interpreting Results

### Log Evidence (Model Selection)

The output includes **Global Log Evidence (LogZ)**:
- **Definition**: Log probability of the model given data
- **Usage**: Compare models by computing Bayes Factor = Δ LogZ

| Δ LogZ | Interpretation |
|--------|----------------|
| > 5 | Strong evidence for more complex model |
| 1-5 | Moderate evidence |
| < 1 | Insufficient evidence (prefer simpler model) |

### Output Format (XMS)

```xml
<motifs>
  <motif id='0' weight='1.0'>
    <weightMatrix length='10'>
      0.95 0.02 0.02 0.01  <!-- Position 1: A=95%, C=2%, G=2%, T=1% -->
      ...
    </weightMatrix>
  </motif>
</motifs>
```

---

## Advanced Topics

### Background Model

The `-bgOrder` option controls sequence composition modeling:

| Order | Parameters | Use Case |
|-------|------------|----------|
| 0 | 4 | Uniform (legacy mode) |
| 1 | 16 | Dinucleotide bias |
| 2 | 64 | Trinucleotide patterns |
| **3** | **256** | **Default (recommended)** |
| 4+ | 1024+ | Large datasets only |

### Checkpointing

For long runs or preemptible instances (AWS Spot):

```bash
# Save every 100 cycles
python3 -m nestedmica.apps.mocca_fast \
  -seqs data.fa -out result.xms \
  -checkpoint state.pkl -checkpointInterval 100

# Resume from checkpoint
python3 -m nestedmica.apps.mocca_fast \
  -seqs data.fa -out result.xms \
  -restart state.pkl
```

---

## Directory Structure

```
nestedmica/
├── apps/
│   ├── discover.py      # Auto-discovery tool
│   └── mocca_fast.py    # Core trainer wrapper
├── model/
│   ├── background.py    # Markov background model
│   └── cython_model.pyx # Optimized DP likelihood
├── trainer/
│   └── fast.py          # Nested Sampling engine
└── utils/
    ├── checkpoint.py    # Save/restore state
    └── console.py       # ASCII visualization
```

---

## License and Credits

Python port of **NestedMICA** (v0.8.0), originally developed by **Thomas Down** and **Genome Research Ltd**.

- **Original Author**: Thomas Down (@derkholm)
- **Copyright**: (c) 2004-2007 Genome Research Ltd.
- **License**: LGPL v2.1+
- **Source**: [Sanger FTP](https://ftp.sanger.ac.uk/pub/resources/software/nmica/)

### Citation

If you use this software, please cite:

> Down TA, Bergman CM. *NestedMICA: sensitive inference of over-represented motifs in nucleic acid sequence.*  
> Nucleic Acids Res. 2007;35(20):e137. doi: 10.1093/nar/gkm787.

