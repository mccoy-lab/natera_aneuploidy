import sys
import numpy as np
import pickle
import gzip as gz
import pandas as pd
from karyohmm import MetaHMM


def convert_geno(gtypes):
    haps = []
    for g in gtypes:
        if g == "BB":
            haps.append([1, 1])
        elif g == "AB":
            haps.append([0, 1])
        elif g == "BA":
            haps.append([0, 1])
        elif g == "AA":
            haps.append([0, 0])
        else:
            haps.append([np.nan, np.nan])
    return np.array(haps).T


def process_ivan_data(csv_fp, chrom="1", sample_id=1):
    """Process data from Ivan."""
    df = pd.read_csv(
        csv_fp,
        dtype={
            "Chr": str,
            "mother_gtype": str,
            "father_gtype": str,
            "anonymized_sample_id": int,
            "anonymized_snp_id": int,
            "b_allele_freq": float,
            "category_num": int,
        },
    )
    cur_df = df[
        (df["Chr"] == chrom) & (df["anonymized_sample_id"] == sample_id)
    ].sort_values("anonymized_snp_id")
    baf = cur_df.b_allele_freq.values
    mat_haps = convert_geno(cur_df.mother_gtype.values)
    pat_haps = convert_geno(cur_df.father_gtype.values)
    # print(chrom, sample_id, baf, mat_haps, pat_haps)
    # Filter to only the indexes that are assessable in the HMM (e.g. no parental missing genotypes)
    if (mat_haps.ndim == 2) and (pat_haps.ndim == 2):
        assert baf.size == mat_haps.shape[1]
        assert baf.size == pat_haps.shape[1]
        idx1 = ~np.isnan(pat_haps[0, :])
        idx2 = ~np.isnan(mat_haps[0, :])
        idx = idx1 & idx2 & (~np.isnan(baf))
        vals, cnts = np.unique(cur_df.category_num.values, return_counts=True)
        prop = cnts / np.sum(cnts)
        categories = {v: p for (v, p) in zip(vals, prop)}
        max_cat = vals[np.argmax(prop)]
        return baf[idx], mat_haps[:, idx], pat_haps[:, idx], categories, max_cat
    else:
        return None, mat_haps, pat_haps, {}, "XX"


def run_meta_hmm_full(baf, mat_haps, pat_haps, unphased=True):
    hmm = MetaHMM()
    # Step 1. Parameter Inference
    if baf is None:
        res_dict = {
            "0": np.nan,
            "1m": np.nan,
            "1p": np.nan,
            "2m": np.nan,
            "2p": np.nan,
            "2": np.nan,
            "3m": np.nan,
            "3p": np.nan,
            "sigma_baf": np.nan,
            "pi0_baf": np.nan,
        }
    else:
        print("Running meta HMM parameter estimation ...", file=sys.stderr)
        pi0_est, sigma_est = hmm.est_sigma_pi0(
            bafs=baf[::5],
            mat_haps=mat_haps[:, ::5],
            pat_haps=pat_haps[:, ::5],
            unphased=unphased,
        )
        print("Finished meta HMM parameter estimation!", file=sys.stderr)
        gammas, states, karyotypes = hmm.forward_backward(
            bafs=baf,
            mat_haps=mat_haps,
            pat_haps=pat_haps,
            pi0=pi0_est,
            std_dev=sigma_est,
            unphased=unphased,
        )
        print(
            "Finished running meta HMM forward-backward algorithm ...",
            file=sys.stderr,
        )
        res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
        res_dict["sigma_baf"] = sigma_est
        res_dict["pi0_baf"] = pi0_est
        res_dict["karyotypes"] = karyotypes
    return res_dict


if __name__ == "__main__":
    # Read in the input data and params ...
    ivan_df = pd.read_csv(
        snakemake.input["ivan_csv"],
        dtype={
            "Chr": str,
            "mother_gtype": str,
            "father_gtype": str,
            "anonymized_sample_id": int,
            "anonymized_snp_id": int,
            "b_allele_freq": float,
            "category_num": int,
        },
    )
    n = int(snakemake.wildcards["n"])
    chroms = [c for c in ivan_df.Chr.unique() if c not in ["X", "Y"]]
    full_chrom_hmm_dict = {}
    for c in chroms:
        full_chrom_hmm_dict[c] = {}
        ind = int(snakemake.wildcards["i"])
        print(f"Processing chr{c} for individual {ind} from family {n}")
        baf, mat_haps, pat_haps, categories, max_cat = process_ivan_data(
            snakemake.input["ivan_csv"], chrom=c, sample_id=ind
        )
        res_dict = run_meta_hmm_full(baf, mat_haps, pat_haps, unphased=True)
        res_dict["ivan_calls"] = categories
        res_dict["ivan_max_cat"] = max_cat
        full_chrom_hmm_dict[c][ind] = res_dict
    pickle.dump(full_chrom_hmm_dict, gz.open(snakemake.output["hmm_pkl"], "wb"))
