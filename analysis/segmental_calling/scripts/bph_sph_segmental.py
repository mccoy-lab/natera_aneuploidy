import gzip as gz
from pathlib import Path

import numpy as np
import polars as pl
from bph_vs_sph import bph_vs_sph_trisomy
from karyohmm import MetaHMM
from tqdm import tqdm
import pickle


def find_hmm_file(base_path="", mother=None, father=None, child=None, chrom=None):
    """Find the HMM traceback file."""
    fname = f"{base_path}/{mother}+{father}/{child}.hmm_model.pkl.gz"
    # assert Path(fname).is_file()
    data = pickle.load(gz.open(fname))
    # assert chrom in data.keys()
    return data


if __name__ == "__main__":
    # Read in the data frame
    seg_df = pl.read_csv(
        snakemake.input["segmental_calls"], null_values=["NA"], separator="\t"
    )
    trisomy_df = seg_df.filter(pl.col("karyotype").is_in(["3m", "3p"]))
    res_data = []
    states = MetaHMM().states
    for mother, father, child, chrom, start, end in tqdm(
        zip(
            trisomy_df["mother"],
            trisomy_df["father"],
            trisomy_df["child"],
            trisomy_df["chrom"],
            trisomy_df["start"],
            trisomy_df["end"],
        )
    ):
        hmm_file = find_hmm_file(
            base_path=snakemake.params["base_path"],
            mother=mother,
            father=father,
            child=child,
            chrom=chrom,
        )
        gammas = hmm_file[chrom]["gammas"]
        pos = hmm_file[chrom]["pos"]
        bf_bph_sph, posterior_bph, posterior_sph = bph_vs_sph_trisomy(
            gammas, states, pos=pos, start=start, end=end
        )
        res_data.append(
            [
                mother,
                father,
                child,
                chrom,
                start,
                end,
                float(bf_bph_sph[0]),
                float(posterior_bph),
                float(posterior_sph),
            ]
        )
    res_dict = {
        "mother": [x[0] for x in res_data],
        "father": [x[1] for x in res_data],
        "child": [x[2] for x in res_data],
        "chrom": [x[3] for x in res_data],
        "start": [x[4] for x in res_data],
        "end": [x[5] for x in res_data],
        "bf_bph_sph": [x[6] for x in res_data],
        "posterior_bph": [x[7] for x in res_data],
        "posterior_sph": [x[8] for x in res_data],
    }
    res_df = pl.from_dict(res_dict, strict=False)
    res_df.write_csv(
        snakemake.output["segmental_bph_sph"], null_value="NA", separator="\t"
    )
