"""
High-level dataset generators for synthetic benchmarks.

This module provides preset dataset generators for common benchmarking
scenarios, combining background, motifs, and planting strategies.
"""

import json
import numpy as np
from typing import List, Dict, Any, Optional, Union
from dataclasses import dataclass, asdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from nestedmica.synthetic.background import (
    BackgroundGenerator, UniformBackground, LearnedBackground,
    LengthDistribution, FixedLength, NormalLength
)
from nestedmica.synthetic.motifs import (
    SyntheticMotif, ConsensusMotif, PWMMotif, GappedSyntheticMotif
)
from nestedmica.synthetic.planting import (
    UniformPlanting, PlantedMotif
)


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


@dataclass
class DatasetConfig:
    """Configuration for a synthetic dataset."""
    n_sequences: int
    seq_length: Union[int, str]  # Fixed or "variable"
    background_type: str
    motifs: List[Dict[str, Any]]
    density: float
    seed: int


class SyntheticDataset:
    """
    Base class for synthetic dataset generation.
    
    Combines background generation, motif definition, and planting
    into a complete benchmark dataset with ground truth.
    """
    
    def __init__(self, n_sequences: int = 100,
                 seq_length: Union[int, LengthDistribution] = 200,
                 background: Optional[BackgroundGenerator] = None,
                 seed: int = 42):
        """
        Initialize dataset generator.
        
        Args:
            n_sequences: Number of sequences to generate.
            seq_length: Fixed length or LengthDistribution.
            background: Background generator (default: UniformBackground).
            seed: Random seed.
        """
        self.n_sequences = n_sequences
        self.seq_length = seq_length if isinstance(seq_length, LengthDistribution) else FixedLength(seq_length)
        self.background = background if background is not None else UniformBackground()
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        
        self.motifs: List[SyntheticMotif] = []
        self.planting = UniformPlanting()
        self.ground_truth: List[Dict[str, Any]] = []
    
    def add_motif(self, motif: Union[SyntheticMotif, str]) -> 'SyntheticDataset':
        """
        Add a motif to plant.
        
        Args:
            motif: SyntheticMotif or consensus string.
            
        Returns:
            Self for chaining.
        """
        if isinstance(motif, str):
            # Check for gapped pattern
            if '[' in motif and ']' in motif:
                motif = GappedSyntheticMotif.from_pattern(motif)
            else:
                motif = ConsensusMotif(motif)
        self.motifs.append(motif)
        return self
    
    def set_planting(self, planting: UniformPlanting) -> 'SyntheticDataset':
        """Set planting strategy."""
        self.planting = planting
        return self
    
    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq))
    
    def generate(self, output_path: str) -> str:
        """
        Generate dataset and write to FASTA.
        
        Args:
            output_path: Path to output FASTA file.
            
        Returns:
            Path to generated file.
        """
        records = []
        self.ground_truth = []
        
        lengths = self.seq_length.sample_n(self.n_sequences, self.rng)
        
        for i in range(self.n_sequences):
            seq_len = int(lengths[i])
            
            # Generate background
            seq = self.background.generate(seq_len, self.rng)
            seq_list = list(seq)
            
            planted_info = {'id': f'seq_{i}', 'planted': []}
            occupied = []
            
            # Plant motifs
            for motif_id, motif in enumerate(self.motifs):
                # Get planting position
                _, planted = self.planting.plant(seq, motif.length, motif_id, self.rng, occupied)
                
                if planted is not None:
                    # Sample motif instance
                    if isinstance(motif, GappedSyntheticMotif):
                        instance, gap_lengths = motif.sample(self.rng)
                        planted.gap_lengths = gap_lengths
                    else:
                        instance = motif.sample(self.rng)
                    
                    # Reverse complement if needed
                    if planted.strand == '-':
                        instance = self._reverse_complement(instance)
                    
                    # Plant into sequence
                    for j, base in enumerate(instance):
                        if planted.start + j < len(seq_list):
                            seq_list[planted.start + j] = base
                    
                    planted.sequence = instance
                    planted_info['planted'].append(asdict(planted))
                    occupied.append((planted.start, planted.end))
            
            final_seq = ''.join(seq_list)
            records.append(SeqRecord(Seq(final_seq), id=f'seq_{i}', description=''))
            self.ground_truth.append(planted_info)
        
        # Write FASTA
        SeqIO.write(records, output_path, 'fasta')
        return output_path
    
    def save_truth(self, output_path: str) -> str:
        """
        Save ground truth to JSON.
        
        Args:
            output_path: Path to output JSON file.
            
        Returns:
            Path to generated file.
        """
        truth = {
            'config': {
                'n_sequences': self.n_sequences,
                'seed': self.seed,
                'n_motifs': len(self.motifs)
            },
            'motifs': [
                {'id': i, 'consensus': m.consensus, 'length': m.length}
                for i, m in enumerate(self.motifs)
            ],
            'sequences': self.ground_truth
        }
        
        with open(output_path, 'w') as f:
            json.dump(truth, f, indent=2, cls=NumpyEncoder)
        
        return output_path


class SimpleBenchmark(SyntheticDataset):
    """
    Simple benchmark with a single strong motif.
    
    Good for basic algorithm validation.
    """
    
    def __init__(self, motif: Union[str, SyntheticMotif] = "ACGTACGT",
                 n_sequences: int = 100, seq_length: int = 200,
                 density: float = 0.8, gc_content: float = 0.5,
                 seed: int = 42, **kwargs):
        super().__init__(
            n_sequences=n_sequences,
            seq_length=kwargs.get('length_dist', seq_length),
            background=UniformBackground(gc_content),
            seed=seed
        )
        self.add_motif(motif)
        self.set_planting(UniformPlanting(density=density))


class GappedMotifBenchmark(SyntheticDataset):
    """
    Benchmark for gapped/bipartite motif discovery.
    
    Uses pattern notation: "ACGT[5-10]TGCA"
    """
    
    def __init__(self, pattern: str = "ACGTACGT[5-10]TGCATGCA",
                 n_sequences: int = 100, seq_length: int = 250,
                 density: float = 0.8, gc_content: float = 0.5,
                 seed: int = 42, **kwargs):
        super().__init__(
            n_sequences=n_sequences,
            seq_length=kwargs.get('length_dist', seq_length),
            background=UniformBackground(gc_content),
            seed=seed
        )
        self.add_motif(GappedSyntheticMotif.from_pattern(pattern))
        self.set_planting(UniformPlanting(density=density))


class MultiMotifBenchmark(SyntheticDataset):
    """
    Benchmark with multiple distinct motifs.
    
    Tests ability to discover multiple motifs simultaneously.
    """
    
    def __init__(self, motifs: List[str] = None,
                 n_sequences: int = 500, seq_length: int = 200,
                 density: float = 0.6, gc_content: float = 0.5,
                 seed: int = 42, **kwargs):
        super().__init__(
            n_sequences=n_sequences,
            seq_length=kwargs.get('length_dist', seq_length),
            background=UniformBackground(gc_content),
            seed=seed
        )
        
        if motifs is None:
            motifs = ["ACGTACGT", "TATAAAT", "GATAAGA"]
        
        for m in motifs:
            self.add_motif(m)
        
        self.set_planting(UniformPlanting(density=density))


class RealisticChIPBenchmark(SyntheticDataset):
    """
    Realistic ChIP-seq-like benchmark.
    
    Features:
    - Variable-length sequences (like real peaks)
    - GC-biased background
    - Optional learned background from real data
    """
    
    def __init__(self, motif: str = "GATAA",
                 n_sequences: int = 500,
                 mean_length: int = 200, std_length: int = 50,
                 density: float = 0.7, gc_content: float = 0.55,
                 background_fasta: Optional[str] = None,
                 seed: int = 42):
        
        length_dist = NormalLength(mean_length, std_length)
        
        if background_fasta:
            background = LearnedBackground.from_fasta(background_fasta)
        else:
            background = UniformBackground(gc_content)
        
        super().__init__(
            n_sequences=n_sequences,
            seq_length=length_dist,
            background=background,
            seed=seed
        )
        
        self.add_motif(motif)
        self.set_planting(UniformPlanting(density=density, strand_bias=0.5))


class BenchmarkGenerator:
    """
    Generator for full benchmark suites with train/test splits.
    
    Organizes datasets into a directory structure:
    output_dir/
      train/
        sequences.fa
        truth.json
      test/
        sequences.fa
        truth.json
    """
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        
    def generate_split(self, dataset: SyntheticDataset, split_name: str,
                       filename: str = "sequences.fa"):
        """
        Generate a single split.
        
        Args:
            dataset: Configured SyntheticDataset instance.
            split_name: Name of split (e.g., 'train', 'test').
            filename: Output filename.
        """
        import os
        split_dir = os.path.join(self.output_dir, split_name)
        os.makedirs(split_dir, exist_ok=True)
        
        filepath = os.path.join(split_dir, filename)
        dataset.generate(filepath)
        
    def create_standard_suite(self, dataset_config: dict, 
                              train_size: int = 1000, test_size: int = 200,
                              benchmark_cls: type = None):
        """
        Create a standard train/test benchmark suit.
        
        Args:
            dataset_config: kwargs for benchmark_cls.
            train_size: Number of training sequences.
            test_size: Number of testing sequences.
            benchmark_cls: Dataset class to use (default: SimpleBenchmark).
        """
        if benchmark_cls is None:
            benchmark_cls = SimpleBenchmark
            
        # Train
        train_config = dataset_config.copy()
        train_config['n_sequences'] = train_size
        train_ds = benchmark_cls(**train_config)
        self.generate_split(train_ds, 'train')
        
        # Test (seed incremented)
        test_config = dataset_config.copy()
        test_config['n_sequences'] = test_size
        
        # Handle seed if present
        if 'seed' in test_config and test_config['seed'] is not None:
            test_config['seed'] += 1
            
        test_ds = benchmark_cls(**test_config)
        self.generate_split(test_ds, 'test')
