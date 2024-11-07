"""
    Export prominent IBD clusters as TPED / TFAM files for plink ... 
"""

import numpy as np
import gzip as gz


with gz.open(
    "../analysis/ibd_mapping/results/natera_parents.chr10.ibd_cluster.ibdclust.gz", "rt"
) as f:
    for i, line in enumerate(f):
        if i == 0:
            header = line.split("\t")[3:]
        else:
            info = line.split("\t")[:3]
            genos = line.split("\t")[3:]
            geno = np.array([[int(i) for i in g.split("|")] for g in genos])
            ibd_clades, counts = np.unique(geno.flatten(), return_counts=True)
            for i in np.where(counts > 5)[0]:
                print(ibd_clades[i], counts[i])
        if i > 30:
            break


def identify_ibd_clusters(genos, clique_size=5):
    """Identifies the IBD clusters and converts them to alleles."""
    assert clique_size > 1
    geno = np.array([[int(i) for i in g.split("|")] for g in genos])
    ibd_clades, counts = np.unique(geno.flatten(), return_counts=True)
    for i in np.where(counts > clique_size)[0]:
        X = geno == ibd_clades[i]
    pass


if __name__ == "__main__":
    with open(snakemake.output["tped"], "w+") as tped_out:
        with gz.open(snakemake.input["ibd_clust"], "rt") as f:
            for i, line in enumerate(f):
                if i == 0:
                    header = line.split("\t")[3:]
                else:
                    info = line.split("\t")[:3]
                    genos = line.split("\t")[3:]
                    X = identify_ibd_clusters(
                        genos, clique_size=int(snakemake.params["clique_size"])
                    )
