# NestedMICA Python


A high-performance Python implementation of the Nested MICA motif discovery algorithm, featuring Cython optimization, parallel processing, and variable-length motif support.

## Features
- **Speed**: Optimized with Cython and multithreading (~22s for 10k sequences vs 7s Java).
- **Variable Length**: Supports `Indel` and `Zap` moves to find optimal motif lengths (unlike the legacy Java release).
- **Checkpointing**: Save/Restore training state for long-running jobs (AWS Spot).
- **Parallelism**: Efficient threaded batch likelihood calculation (GIL-released).

## Installation

Requires Python 3.9+ and a C compiler.

```bash
pip install numpy biopython cython
python3 setup.py build_ext --inplace
```

## Usage

### Fast Parallel Run
```bash
python3 -m nestedmica.apps.mocca_fast \
  -seqs complex_test.fa \
  -numMotifs 5 \
  -motifLength 10 \
  -maxCycles 5000 \
  -ensembleSize 50 \
  -threads 8 \
  -out output.xms
```

### Deep Convergence (with Checkpointing)
```bash
python3 -m nestedmica.apps.mocca_fast \
  -seqs large.fa \
  -stopIqr 0.001 \
  -checkpoint run.pkl \
  -out deep_run.xms
```

## Directory Structure
- `nestedmica/`: Python source package.
- `source/`: Original Java source code (patched for baseline comparison).
- `validate_complex.py`: Verification scripts.

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

This Python implementation preserves the core algorithms (Nested Sampling, Mosaic Background) while introducing new optimization techniques (Cython, Threading) and variable-length sampling improvements.

