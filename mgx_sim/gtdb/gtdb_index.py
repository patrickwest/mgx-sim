from typing import List, Union

from pathlib import Path

import pandas as pd

from mgx_sim import GenomeIndex


class GTDBr95Index(GenomeIndex):
    def __init__(self, root_path: Union[str, Path]):
        super().__init__()
        self.df = self.build_index(root_path)

    def build_index(self, root_path: Union[str, Path]) -> pd.DataFrame:
        root_path = Path(root_path)
        if not root_path.exists():
            raise FileNotFoundError(root_path)
        metadata_path = root_path / Path("metadata") / Path("genome_metadata.tsv")
        metadata = pd.read_csv(metadata_path, sep="\t")

        tree_path = root_path / Path("taxonomy") / Path("gtdb_taxonomy.tsv")
        raw_tree = pd.read_csv(tree_path, sep="\t", header=None)

        tree = self._clean_raw_tree(raw_tree)

        # Join and create absolute path
        merged_data = pd.merge(metadata, tree, how="outer", on="accession")
        paths = [
            f"{root_path}/fastani/database/{accession}_genomic.fna.gz"
            for accession in merged_data.accession
        ]
        merged_data["genome_path"] = paths
        return merged_data

    def _clean_raw_tree(self, raw_tree: pd.DataFrame) -> pd.DataFrame:
        # Remove GB prefix in accessions and expand out all ranks into
        # separate columns.

        def split_tax_string(s: str) -> List[str]:
            separate = s.split(";")
            no_prefix = [x[3:] for x in separate]
            return no_prefix

        split_taxa_dict = raw_tree[1].apply(split_tax_string).to_dict()
        clean = pd.DataFrame(split_taxa_dict).T
        clean.columns = [
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]

        no_prefix_accession = [x[3:] for x in raw_tree[0]]
        clean["accession"] = no_prefix_accession
        return clean
