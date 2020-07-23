# mgx-sim

mgx-sim provides tools for accessing sets of genomes based on metadata, simulating reads or contigs from those genomes, and feeding the results into other tools. We aim to add minimal overhead and enable easy reproducibility.

## Genome Libraries

We refer to a collection of genomes and metadata about those genomes as a "Genome Library". The central component of a Genome Library is a `GenomeLibraryIndex`, which contains the following data in a pandas DataFrame,

1. Absolute path to every genome fasta file in the library.
2. Any metadata associated with every genome.

## Creating a Genome Library

To create a Genome Library, all you need to do is create an index for it! This is done by subclassing `GenomeLibraryIndex` and implementing `build_index`. 

```python
import pandas as pd

from mgxsim import GenomeLibraryIndex


class MyLibraryIndex(GenomeLibraryIndex):

    def build_index(self, path_to_metadata = "data/metadata.csv"):
        # Assumes metadata.csv is already formatted with absolute paths.
        df = pd.read_csv(path_to_metadata)
        return df
```

The `build_index` method creates a pandas DataFrame containing all metadata of interest. Every `GenomeLibraryIndex`'s DataFrame MUST have a column called "genome_path" containing the absolute path to a fasta file for that genome. Everything else is up to you.

## Using a Genome Library

Once you've followed the steps above, you can subset, simulate, and more. The basic objects we provide are

- `GenomeLibraryIndex`: Wraps a `pandas.DataFrame` of all metadata. Has an absolute path for every genome.
- `Genome`: Wraps a list of `Bio.SeqRecord`s for each contig in a genome.
- `GenomeLibrary`: Wraps a list of `Genome`s as well as associated metadata in the associated `GenomeLibraryIndex`.

### Subsetting

Any function from the metadata DataFrame to a list of genomes can be used to extract genomes from a library. This is done with `GenomeLibraryIndex.query`, which returns a `GenomeLibrary`.

```python
from typing import List

import matplotlib.pyplot as plt
import pandas as pd

# Subclass GenomeLibraryIndex in my_lib_idx.py
from my_lib_idx import MyLibraryIndex


def return_archaea(metadata: pd.DataFrame) -> List[str]:
    archaea_idx = metadata.domain == "archaea"
    all_archaea = metadata[archaea_idx]
    return all_archaea.genome_path.tolist()


my_index = MyLibraryIndex()
genome_subset = my_index.query(return_archaea)

# For example, plot histogram of lengths
lengths = [len(g) for g in genome_subsets]
plt.hist(lengths)
plt.show()
```

### Simulating Reads

To simulate reads, we suggest using existing read simulators. Here's an example with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq)

```python
from iss.abundance import lognormal
from iss.error_models.basic import BasicErrorModel
from iss.generator import simulate_read

from mgxsim import MyLibraryIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = MyLibraryIndex()
genomes = my_index.query(list_of_paths)

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

We provide simple contig simulation utils directly. For example, the function `random_fixed_len_contig` samples a random position and returns a contig of fixed length.

```python
import random

from mgxsim.sim.contigs import random_fixed_length_contig

from my_lib_idx import MyLibraryIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = MyLibraryIndex()
genomes = my_index.query(list_of_paths)

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
from mgxsim.sim.distributions import poisson_diffcov
from mgxsim.sim.contigs import random_fixed_length_contig

from my_lib_idx import MyLibraryIndex


def my_query(metadata):
    # do something ...
    return list_of_paths


my_index = MyLibraryIndex()
genomes = my_index.query(list_of_paths)

num_contigs = 1e5
abunds = poisson_diffcov(genomes)

# Continue simulating reads, contigs, what have you..
```

## What mgx-sim is not for

This library is intended for relatively small-scale, self-contained projects. It is not a replacement for a proper database!
