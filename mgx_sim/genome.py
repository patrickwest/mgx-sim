from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Genome:
    def __init__(self, abs_path):
        self._abs_path = str(abs_path)
        self.contigs = list(SeqIO.parse(self.abs_path, "fasta"))

    @property
    def abs_path(self) -> str:
        return self._abs_path

    @property
    def num_contigs(self) -> int:
        return len(self.contigs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)
