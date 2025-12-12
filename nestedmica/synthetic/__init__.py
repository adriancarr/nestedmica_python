"""
NestedMICA Synthetic Data Package.

This package provides tools for generating synthetic DNA sequences with
planted motifs for benchmarking and validation.

Key features:
- Configurable background models (uniform, Markov, learned from FASTA)
- Variable-length sequence generation
- Motif import from JASPAR, HOCOMOCO, MEME formats
- Support for gapped/bipartite motifs
- Ground truth export in JSON format
"""

from nestedmica.synthetic.background import (
    BackgroundGenerator, UniformBackground, MarkovBackground,
    LearnedBackground, ShuffledBackground,
    LengthDistribution, FixedLength, NormalLength, LogNormalLength, LearnedLength
)

from nestedmica.synthetic.motifs import (
    ConsensusMotif,
    PWMMotif,
    # GappedSyntheticMotif,
)

from nestedmica.synthetic.planting import (
    UniformPlanting, PlantedMotif, ClusteredPlanting, PositionalPlanting
)
from nestedmica.synthetic.datasets import (
    SyntheticDataset, SimpleBenchmark, GappedMotifBenchmark,
    MultiMotifBenchmark, RealisticChIPBenchmark, BenchmarkGenerator
)

__all__ = [
    # Background
    'UniformBackground',
    'MarkovBackground', 
    'LearnedBackground',
    'NormalLength',
    'LogNormalLength',
    'LearnedLength',
    # Motifs
    'ConsensusMotif',
    'PWMMotif',
    # Planting
    'UniformPlanting',
    # Datasets
    'SimpleBenchmark',
    'GappedMotifBenchmark',
]
