try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Tutorial requires matplotlib!")
    raise ImportError

from mgxsim.gtdb import GTDBr95Index

def extract_complete(metadata):
    idx = metadata.checkm_completeness >= 99.
    sub_df = metadata[idx]
    return sub_df.genome_path.tolist()[:100]


def extract_incomplete(metadata):
    idx = metadata.checkm_completeness < 55.
    sub_df = metadata[idx]
    return sub_df.genome_path.tolist()[:100]

gtdb_root = "/shared/db/gtdb/latest/"
index = GTDBr95Index(gtdb_root)
complete_genomes = index.query(extract_complete)
incomplete_genomes = index.query(extract_incomplete)

plt.figure(figsize=(10,10))
plt.hist([len(g) for g in complete_genomes], label="complete")
plt.hist([len(g) for g in incomplete_genomes], label="incomplete")
plt.xlabel("Length of genome")
plt.legend()
plt.savefig("example.png")
plt.close()
