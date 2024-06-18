"""Running the Karyohmm MetaHMM for BAF in corresponding embryos."""

import sys
import pandas as pd
import numpy as np
import gzip as gz
import ruptures as rpt
from tqdm import tqdm
from karyohmm import MetaHMM

if __name__ == "__main__":
    hmm = MetaHMM()
    # 1. Read in the actual data
    data = np.load(gz.open(snakemake.input["hmm_trace"]), allow_pickle=True)
    intervals = []
    for c in tqdm(snakemake.params["chroms"]):
        if "pos" not in data[c]:
            intervals.append(
                [
                    snakemake.wildcards["mother"],
                    snakemake.wildcards["father"],
                    snakemake.wildcards["child"],
                    c,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            )
        else:
            pos = data[c]["pos"]
            gamma_rev, karyo = hmm.marginal_posterior_karyotypes(
                data[c]["gammas"], data[c]["karyotypes"]
            )
            pts = rpt.KernelCPD(min_size=snakemake.params["min_size"]).fit_predict(
                gamma_rev.T, pen=snakemake.params["penalty"]
            )
            # 2. Now for each segment we need to quantify the posterior distribution and the number of SNPs
            if len(pts) > 1:
                for s, e in zip(pts[:-1], pts[1:]):
                    cur_mean_gamma = gamma_rev[:, (s - 1) : (e - 1)].mean(axis=1)
                    # Only take the segments that are not euploidies ...
                    if karyo[np.argmax(cur_mean_gamma)] != "2":
                        intervals.append(
                            [
                                snakemake.wildcards["mother"],
                                snakemake.wildcards["father"],
                                snakemake.wildcards["child"],
                                c,
                                pos[s - 1],
                                pos[e - 1],
                                karyo[np.argmax(cur_mean_gamma)],
                                np.max(cur_mean_gamma),
                                e - s,
                            ]
                        )
    # 3. Write out the values to a dataframe
    if len(intervals) == 0:
        intervals.append(
            [
                snakemake.wildcards["mother"],
                snakemake.wildcards["father"],
                snakemake.wildcards["child"],
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ]
        )
    tot_df = pd.DataFrame(intervals)
    tot_df.columns = [
        "mother",
        "father",
        "child",
        "chrom",
        "start",
        "end",
        "karyotype",
        "mean_posterior",
        "nsnps",
    ]
    tot_df.to_csv(snakemake.output["changept_calls"], sep="\t", index=None)
