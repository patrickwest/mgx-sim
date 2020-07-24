from typing import Iterator

import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Genome:
    def __init__(self, abs_path):
        self._abs_path = str(abs_path)
        with gzip.open(self.abs_path, "rt") as handle:
            self.contigs = list(SeqIO.parse(handle, "fasta"))

    @property
    def abs_path(self) -> str:
        return self._abs_path

    @property
    def num_contigs(self) -> int:
        return len(self.contigs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)
