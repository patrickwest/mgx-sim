from typing import Callable, List, Union

from abc import abstractmethod
from pathlib import Path

import pandas as pd

from .genome import Genome


class GenomeIndex:
    def __init__(self):
        self.df = None

    @abstractmethod
    def build_index(self) -> pd.DataFrame:
        pass

    def __getitem__(self, genome: Genome) -> pd.Series:
        return self.df[self.df.genome_path == genome.abs_path]

    def query(
        self, query_fn: Callable[[pd.DataFrame], List[Union[Path, str]]]
    ) -> List[Genome]:
        genome_paths = query_fn(self.df)
        genomes = [Genome(path) for path in genome_paths]
        return genomes

    @property
    def dataframe(self) -> pd.DataFrame:
        return self.df
