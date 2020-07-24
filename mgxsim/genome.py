from typing import Iterator

import binascii
import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def is_gz_file(filepath):
    with open(filepath, "rb") as f:
        return binascii.hexlify(f.read(2)) == b"1f8b"


class Genome:
    def __init__(self, abs_path):
        self._abs_path = str(abs_path)
        if is_gz_file(abs_path):
            with gzip.open(self.abs_path, "rt") as handle:
                self.contigs = list(SeqIO.parse(handle, "fasta"))
        else:
            self.contigs = list(SeqIO.parse(abs_path, "fasta"))

    @property
    def abs_path(self) -> str:
        return self._abs_path

    @property
    def num_contigs(self) -> int:
        return len(self.contigs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)
    
    def __len__(self) -> int:
        num = sum([len(c) for c in self])
        return num
