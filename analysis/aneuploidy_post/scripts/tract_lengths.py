import numpy as np
import pickle
import gzip as gz
import pandas as pd
from karyohmm import MetaHMM


def get_tract_dist(gammas, karyo):
    """Compute summaries of the chunk-length distribution for each call."""
    max_post_path = np.argmax(gammas, axis=0)
    changepts = np.where(np.diff(max_post_path) != 0)[0]
    if 0 not in changepts:
        changepts = np.insert(changepts, 0, 0)
    changepts = np.append(changepts, gammas.shape[1])
    assert changepts.size >= 2
    chunk_dict = {}
    for s, e in zip(changepts[:-1], changepts[1:]):
        k = karyo[max_post_path[s:e][-1]]
        nsnps = e - s
        if k not in chunk_dict:
            chunk_dict[k] = [nsnps]
        else:
            chunk_dict[k] = chunk_dict[k] + [nsnps]
    return changepts, chunk_dict


if __name__ == "__main__":
    # Read in the HMM data
    chroms = [f"chr{i}" for i in range(1, 23)]
    hmm = MetaHMM()
    # Read in the hmm data + centromere data
    hmm_traceback = np.load(
        gz.open(snakemake.input["hmm_traceback"], "r"), allow_pickle=True
    )
    n_snps = []
    n_chunks = []
    mean_chunk_size = []
    for c in chroms:
        # NOTE: gammas are assumed to be in log-space here
        X, karyo = hmm.marginal_posterior_karyotypes(
            hmm_traceback[c]["gammas"], hmm_traceback[c]["karyotypes"]
        )
        n_snps.append(X.shape[1])
        _, chunk_dict = get_tract_dist(X, karyo)
        cur_chunks = []
        cur_mean_chunk_size = []
        for k in karyo:
            if k in chunk_dict:
                cur_chunks.append(len(chunk_dict[k]))
                cur_mean_chunk_size.append(np.mean(chunk_dict[k]))
            else:
                cur_chunks.append(np.nan)
                cur_mean_chunk_size.append(np.nan)
        n_chunks.append(cur_chunks)
        mean_chunk_size.append(cur_mean_chunk_size)
    # Creating the full dictionary for this set of parents
    res_dict = {
        "mother": np.repeat(snakemake.wildcards["mother"], len(chroms)),
        "father": np.repeat(snakemake.wildcards["father"], len(chroms)),
        "child": np.repeat(snakemake.wildcards["child"], len(chroms)),
        "chrom": chroms,
        "nsnps": n_snps,
    }
    for i, k in enumerate(karyo):
        res_dict[f"n_tract_{k}"] = [n_chunks[j][i] for j in range(len(chroms))]
        res_dict[f"mean_tract_{k}"] = [
            mean_chunk_size[j][i] for j in range(len(chroms))
        ]
    # Convert to a pandas DataFrame
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(snakemake.output["tracts_tsv"], na_rep="NA", sep="\t", index=None)
