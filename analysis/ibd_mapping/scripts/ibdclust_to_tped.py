"""
    Export prominent IBD clusters as TPED / TFAM files for plink ... 
"""

import numpy as np
import gzip as gz
from tqdm import tqdm


def identify_ibd_clusters(genos, clique_size=5):
    """Identifies the IBD clusters and converts them to alleles."""
    assert clique_size > 1
    geno = np.array([[int(i) for i in g.split("|")] for g in genos])
    ibd_clades, counts = np.unique(geno.flatten(), return_counts=True)
    tped_str = []
    for i in np.where(counts >= clique_size)[0]:
        X = (geno == ibd_clades[i]).astype(np.uint8)
        new_geno = [f"{x1}\t{x2}" for [x1, x2] in X]
        tped_str.append("\t".join(new_geno))
    return tped_str


if __name__ == "__main__":
    with open(snakemake.output["tped"], "w+") as tped_out:
        with gz.open(snakemake.input["ibd_clust"], "rt") as f:
            for i, line in tqdm(enumerate(f)):
                if i == 0:
                    header = line.split("\t")[3:]
                    with open(snakemake.output["tfam"], "w+") as tfam_out:
                        for x in header:
                            tfam_out.write(f"{x}\t{x}\t0\t0\t0\t0\n")
                else:
                    info = line.split("\t")[:3]
                    genos = line.split("\t")[3:]
                    tped_str = identify_ibd_clusters(
                        genos, clique_size=int(snakemake.params["clique_size"])
                    )
                    for i, s in enumerate(tped_str):
                        tped_out.write(
                            f"{info[0]}\t{info[0]}:{info[1]}-{i}\t{info[2]}\t{info[1]}\t{tped_str}"
                        )
