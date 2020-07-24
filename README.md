# mgxsim

mgxsim provides tools for accessing sets of genomes based on metadata, simulating reads or contigs from those genomes, and feeding the results into other tools. Our goals are minimal overhead and easy reproducibility.

TODOs:

1. Make subclassing and validation easier.
2. Remove the triple nested comprehensions necessary for simulation examples.

## Installation

Using pip, just run

```bash
pip install .
```

## Representing a Set of Genomes

We provide two core classes for manipulating sets of genomes. These are 
1. `GenomeIndex`: Provides a mapping from metadata to absolute path for every genome.
2. `Genome`: Represents a genome as a container of contigs.

## Adding a New Set of Genomes

To adapt mgxsim for your genome set, all you need to do is create a custom `GenomeIndex`. This is done by subclassing `GenomeIndex` and implementing `load_metadata`.

```python
import pandas as pd

from mgxsim import GenomeIndex


class TrivialIndex(GenomeIndex):

    def __init__(self, path_to_medata = "data/metadata.csv"):
        super().__init__()
        self.metadata = build_index(path_to_metadata)

    def load_metadata(self, path_to_metadata):
        # Assumes metadata.csv is already formatted with absolute paths.
        df = pd.read_csv(path_to_metadata)
        return df
```

Note that every subclass of `GenomeIndex` must create an attribute `self.metadata` and this DataFrame MUST have a column called "genome_path" containing the absolute path to a fasta file for that genome. Everything else is up to you.

For a more complicated example, check out `mgxsim/gtdb/gtdb_index.py` for a [GTDB](https://gtdb.ecogenomic.org/) index.


## Using a GenomeIndex

Once you've followed the steps above, you can subset, simulate, and more. We show a few common use cases below. Check out the `examples` directory for more!

### Subsetting

Any function from metadata to a list of genome paths can be used to load genomes. This is done with `GenomeIndex.query`, which returns a list of `Genomes`s.

```python
from typing import List

import matplotlib.pyplot as plt
import pandas as pd

# Subclass GenomeIndex in triv_idx.py
from triv_idx import TrivialIndex


def return_archaea(metadata: pd.DataFrame) -> List[str]:
    archaea_idx = metadata.domain == "archaea"
    all_archaea = metadata[archaea_idx]
    return all_archaea.genome_path.tolist()


my_index = TrivialIndex()
genome_subset = my_index.query(return_archaea)

# For example, plot histogram of lengths
lengths = [len(g) for g in genome_subset]
plt.hist(lengths)
plt.show()
```

### Simulating Reads

To simulate reads, we suggest using existing read simulators. Here's an example with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)

```python
from iss.abundance import lognormal
from iss.error_models.basic import BasicErrorModel
from iss.generator import simulate_read

from triv_idx import TrivialIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = TrivialIndex()
genomes = my_index.query(my_query)

# Simulate abundances per contig
num_reads = 1e9
abundances = lognormal(genomes.flattened_contigs())
error_model = BasicErrorModel()

# Nested comprehension says for all genomes, for all contigs in genome, simulate reads.
# We should write a util function to hide the triple comprehension.
reads = [
    [
        [
            simulate_read(contig, error_model) for _ in range(abundances[contig] * num_reads)
        ]
        for contig in genome
    ] 
    for genome in genomes
]

genome_5_contig_4_reads = reads[5][4]
```

### Simulating Contigs

We provide simple contig simulation functions directly. For example, the function `random_fixed_len_contig` samples a random position and returns a contig of fixed length.

```python
import random

from mgxsim.sample import random_fixed_length_contig

from triv_idx import TrivialIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = TrivialIndex()
genomes = my_index.query(my_query)

num_contigs = 1e5
counts = {c: int(random.random() * num_contigs) for c in genomes.flattened_contigs()}


# Again, we should write a util function to hide the triple comprehension.
sim_contigs = [
    [
        [
            random_fixed_len_contig(contig, 1000) for _ in range(counts[contig])
        ]
        for contig in genome
    ] 
    for genome in genomes
]
```

### Simulating Differential Coverage

Simulating correlated coverage for all contigs in a genome can be useful. We provide some distributions for handling this. The recommended one is `poisson_diffcov`. This generates a random lambda > 0 for every genome, then samples coverage for each contig from a Poisson with mean lambda * length of contig.

```python
from mgxsim.distributions import poisson_diffcov
from mgxsim.sample import random_fixed_length_contig

from triv_idx import TrivialIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = TrivialIndex()
genomes = my_index.query(my_query)

num_contigs = 1e5
abunds = poisson_diffcov(genomes)

# Continue simulating reads, contigs, what have you..
```

## What mgxsim is not for

This library is intended for relatively small-scale, self-contained projects. It is not a replacement for a proper database!
