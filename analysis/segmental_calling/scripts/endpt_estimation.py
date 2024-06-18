import gzip as gz
from pathlib import Path

import numpy as np
import polars as pl
import scipy as sp
from karyohmm import MetaHMM
from tqdm import tqdm


def find_hmm_file(base_path="", mother=None, father=None, child=None, chrom=None):
    """Find the HMM traceback file."""
    fname = f"{base_path}/{mother}+{father}/{child}.hmm_model.pkl.gz"
    assert Path(fname).is_file()
    data = np.load(gz.open(fname), allow_pickle=True)
    assert chrom in data.keys()
    return data


def estimate_endpoints(
    gammas, pos, karyo, karyotype="3p", start=None, end=None, qs=[0.05, 0.5, 0.95]
):
    assert gammas.shape[1] == pos.size
    assert end >= start
    # Look at the starting point first.
    hmm = MetaHMM()
    g, x = hmm.marginal_posterior_karyotypes(gammas, karyo)
    assert karyotype in x
    start_idx = np.where(pos == start)[0][0]
    end_idx = np.where(pos == end)[0][0]
    start_i, start_e = int(np.maximum(start_idx - 100, 0)), int(
        np.minimum(start_idx + 100, pos.size)
    )
    end_i, end_e = int(np.maximum(end_idx - 100, 0)), int(
        np.minimum(end_idx + 100, pos.size)
    )
    post_path_start = g[np.where(x == karyotype)[0], start_i:start_e]
    post_path_end = g[np.where(x == karyotype)[0], end_i:end_e]
    start_pts = []
    end_pts = []
    for q in qs:
        idx_start = np.where(post_path_start >= q)[0]
        idx_end = np.where(post_path_end >= q)[0]
        x_start = np.arange(start_i, start_e)
        x_end = np.arange(end_i, end_e)
        if idx_start.size > 0:
            start_pts.append(pos[x_start[int(np.min(idx_start))]])
        else:
            start_pts.append(np.nan)
        if idx_end.size > 0:
            end_pts.append(pos[x_end[int(np.max(idx_end))]])
        else:
            end_pts.append(np.nan)
    return start_pts, end_pts


if __name__ == "__main__":
    # Read in the data frame of segmental calls
    seg_df = pl.read_csv(
        snakemake.input["segmental_calls"], null_values=["NA"], separator="\t"
    )
    hmm = MetaHMM()
    res_data = []
    for mother, father, child, chrom, start, end, karyotype in tqdm(
        seg_df[
            ["mother", "father", "child", "chrom", "start", "end", "karyotype"]
        ].iter_rows()
    ):
        hmm_file = find_hmm_file(
            base_path=snakemake.params["base_path"],
            mother=mother,
            father=father,
            child=child,
            chrom=chrom,
        )

        states = hmm.states
        gammas = hmm_file[chrom]["gammas"]
        karyo = hmm_file[chrom]["karyotypes"]

        pos = hmm_file[chrom]["pos"]
        start_epts, end_epts = estimate_endpoints(
            gammas,
            pos,
            karyo,
            karyotype=karyotype,
            start=start,
            end=end,
            qs=[0.05, 0.5, 0.95],
        )
        res_data.append(
            [mother, father, child, chrom, start, end, karyotype]
            + start_epts
            + end_epts
        )
    res_dict = {
        "mother": [x[0] for x in res_data],
        "father": [x[1] for x in res_data],
        "child": [x[2] for x in res_data],
        "chrom": [x[3] for x in res_data],
        "start": [x[4] for x in res_data],
        "end": [x[5] for x in res_data],
        "karyotype": [x[6] for x in res_data],
        "start_05": [x[7] for x in res_data],
        "start_50": [x[8] for x in res_data],
        "start_95": [x[9] for x in res_data],
        "end_05": [x[10] for x in res_data],
        "end_50": [x[11] for x in res_data],
        "end_95": [x[12] for x in res_data],
    }
    res_df = pl.from_dict(res_dict, strict=False)
    res_df.write_csv(
        snakemake.output["segmental_endpoints"], null_value="NA", separator="\t"
    )
