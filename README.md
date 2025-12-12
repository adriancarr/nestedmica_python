# NestedMICA Python

A high-performance Python implementation of the Nested MICA motif discovery algorithm, featuring Cython optimization, parallel processing, and rigorous Bayesian model selection.

## Features
- **Auto-Discovery**: Automatically determines the optimal number of motifs ($N$) and their lengths using Bayesian Evidence.
- **Speed**: Optimized with Cython and multithreading (~22s for 10k sequences vs 7s Java).
- **Variable Length**: Supports `Insert` and `Delete` moves for dynamic motif sizing.
- **Rigorous Significance**: Calculates Global Log Evidence ($\ln Z$) for model comparison (Bayes Factors).
- **Smart Sizing**: Auto-configures ensemble size based on dataset magnitude.
- **Checkpointing**: Save/Restore training state for long-running jobs (AWS Spot).

## Installation

Requires Python 3.9+ and a C compiler.

```bash
pip install numpy biopython cython
python3 setup.py build_ext --inplace
```

## Usage

### 1. Auto-Discovery (Recommended)
Let the algorithm decide the number of motifs and their size.

```bash
python3 -m nestedmica.apps.discover \
  -seqs data.fa \
  -out discovered.xms \
  -threads 8
```

**Optional Arguments:**
- `-minMotifs 3`: Start searching from N=3.
- `-maxMotifs 10`: Stop searching at N=10.
- `-ensembleSize 100`: Manually set precision (Default: Auto-calculated like $50 \times \log_{10}(\text{bases})$).

### 2. Manual Run (Advanced)
If you know exactly what you want (e.g. "Find exactly 5 motifs of length 10").

```bash
python3 -m nestedmica.apps.mocca_fast \
  -seqs data.fa \
  -numMotifs 5 \
  -motifLength 10 \
  -thread 8 \
  -out output.xms
```

### 3. Interpreting Results (Log Evidence)
The output now includes the **Global Log Evidence (`LogZ`)**.
- **Definition**: The probability of the *entire model* given the data.
- **Usage**: To compare two runs (e.g. 4 motifs vs 5 motifs), subtract their LogZ values to get the **Bayes Factor**.
    - $\Delta \ln Z > 5$: Strong evidence for the more complex model.
    - $\Delta \ln Z < 0$: The extra motif is just noise.

## Directory Structure
- `nestedmica/`: Python source package.
    - `apps/discover.py`: Auto-discovery tool.
    - `apps/mocca_fast.py`: Core trainer wrapper.
    - `trainer/fast.py`: Cython-optimized Nested Sampling engine.
- `source/`: Legacy Java source code (reference).

## License and Credits

This software is a Python port of **NestedMICA** (v0.8.0), originally developed by **Thomas Down** and **Genome Research Ltd**.

### Original Software
- **Author**: Thomas Down (@derkholm)
- **Copyright**: (c) 2004-2007 Genome Research Ltd.
- **License**: LGPL v2.1+
- **Source**: [Sanger FTP](https://ftp.sanger.ac.uk/pub/resources/software/nmica/)

### Citation
If you use this software, please cite the original NestedMICA publication:
> Down TA, Bergman CM. *NestedMICA: sensitive inference of over-represented motifs in nucleic acid sequence.*  
> Nucleic Acids Res. 2007;35(20):e137. doi: 10.1093/nar/gkm787.
