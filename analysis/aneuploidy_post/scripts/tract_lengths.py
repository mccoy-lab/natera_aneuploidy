import numpy as np
import pickle
import gzip as gz
import click
import pandas as pd
from karyohmm import MetaHMM

if __name__ == "__main__":
    # Read in the BAF data
    chroms = [f"chr{i}" for i in range(1, 23)]
    hmm = MetaHMM()
    # Read in the hmm data + centromere data
    hmm_traceback = np.load(
        gz.open(snakemake.input["hmm_traceback"], "r"), allow_pickle=True
    )
    n_tracts = []
    max_tract_frac = []
    for c in chroms:
        # NOTE: gammas are assumed to be in log-space here
        X, karyo = hmm.marginal_posterior_karyotypes(
            hmm_traceback[c]["gammas"], hmm_traceback[c]["karyotypes"]
        )
        map_path = np.argmax(X, axis=0)
        pos = hmm_data[c]["pos"]
        idx = map_path[1:] != map_path[:-1]
        dists = np.array([pos[0]] + [pos[i] for i in np.where(idx)[0]] + [pos[-1]])
        tract_fractions = (dists[1:] - dists[:-1]) / (pos[-1] - pos[0])
        n_tracts.append(tract_fractions.size)
        max_tract_frac.append(np.max(tract_fractions))
    res_dict = {
        "mother": np.repeat(snakemake.wildcards["mother"], len(chroms)),
        "father": np.repeat(snakemake.wildcards["father"], len(chroms)),
        "child": np.repeat(snakemake.wildcards["child"], len(chroms)),
        "chrom": chroms,
        "ntracts": n_tracts,
        "max_tract_frac": max_tract_frac,
    }
    # Convert to a pandas DataFrame
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(snakemake.output["tracts_tsv"], sep="\t", index=None)
