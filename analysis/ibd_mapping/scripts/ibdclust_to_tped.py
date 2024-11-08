"""
    Export prominent IBD clusters as TPED / TFAM files for plink ... 
"""

import numpy as np
import gzip as gz
from tqdm import tqdm


def identify_ibd_clusters(genos, clique_size=5, max_clique_size=20, max_n_cliques=10):
    """Identifies the IBD clusters and converts them to alleles."""
    assert clique_size > 1
    geno = np.array([[int(i) for i in g.split("|")] for g in genos])
    ibd_clades, counts = np.unique(geno.flatten(), return_counts=True)
    tped_str = []
    idx = np.where(counts >= clique_size)[0]
    if idx.size <= max_n_cliques:
        for i in idx:
            if counts[i] <= max_clique_size:
                X = (geno == ibd_clades[i]).astype(np.uint8)
                new_geno = "\t".join([f"{x1+1}\t{x2+1}" for (x1, x2) in X])
                tped_str.append(new_geno)
    return tped_str


if __name__ == "__main__":
    with open(snakemake.output["tped"], "w+") as tped_out:
        with gz.open(snakemake.input["ibd_clust"], "rt") as f:
            for i, line in tqdm(enumerate(f)):
                if i == 0:
                    header = line.rstrip().split("\t")[3:]
                    with open(snakemake.output["tfam"], "w+") as tfam_out:
                        for x in header:
                            tfam_out.write(f"{x}\t{x}\t0\t0\t0\t0\n")
                else:
                    info = line.split("\t")[:3]
                    genos = line.rstrip().split("\t")[3:]
                    tped_str = identify_ibd_clusters(
                        genos, clique_size=int(snakemake.params["clique_size"])
                    )
                    for j, s in enumerate(tped_str):
                        assert len(s.split("\t")) == 2 * len(header)
                        tped_out.write(
                            f"{info[0][3:]}\t{info[0]}_{info[1]}_{j}\t{info[2]}\t{info[1]}\t{s}\n"
                        )
