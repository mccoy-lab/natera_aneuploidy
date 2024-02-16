import numpy as np
import pickle
import gzip as gz
import click
import pandas as pd
from karyohmm import MosaicEst


def est_gain_loss(df, mother, father, child, chrom):
    """Determine if this is likely a mosaic gain or loss."""
    assert "mother" in df.columns
    assert "father" in df.columns
    assert "child" in df.columns
    assert "chrom" in df.columns
    assert "3m" in df.columns
    assert "3p" in df.columns
    assert "1m" in df.columns
    assert "1p" in df.columns
    assert "2" in df.columns
    assert "0" in df.columns
    assert "sigma_baf" in df.columns
    aneu_cats = ["0", "1m", "1p", "2", "3m", "3p"]
    aneu_cat_posterior = df[
        (df.mother == mother)
        & (df.father == father)
        & (df.child == child)
        & (df.chrom == chrom)
    ][aneu_cats].values[0]
    assert np.isclose(np.sum(aneu_cat_posterior), 1.0)
    # NOTE: to determine gain vs. loss we look at the relative order
    idx = np.argsort(aneu_cat_posterior)
    sort_aneu = np.array(aneu_cats, dtype=str)[idx[::-1]]
    gain = None
    if sort_aneu[0] == "2":
        if (sort_aneu[1] == "1m") or (sort_aneu[1] == "1p"):
            gain = False
        if (sort_aneu[1] == "3m") or (sort_aneu[1] == "3p"):
            gain = True
    if (sort_aneu[0] == "1m") or (sort_aneu[0] == "1p"):
        if sort_aneu[1] == "2":
            gain = False
    if (sort_aneu[0] == "3m") or (sort_aneu[0] == "3p"):
        if sort_aneu[1] == "2":
            gain = True
    return gain


if __name__ == "__main__":
    # 1. Read in the mosaic dataframe
    mother = snakemake.wildcards["mother"]
    father = snakemake.wildcards["father"]
    child = snakemake.wildcards["child"]
    chrom = snakemake.wildcards["chrom"]
    mosaic_df = pd.read_csv(snakemake.input["mosaic_tsv"], sep="\t")
    gain = est_gain_loss(
        mosaic_df, mother=mother, father=father, child=child, chrom=chrom
    )
    try:
        if gain is not None:
            # 2. Read in the BAF data
            data = pickle.load(gz.open(snakemake.input["baf_pkl"], "r"))
            pos = data[chrom]["pos"]
            mat_haps = data[chrom]["mat_haps"]
            pat_haps = data[chrom]["pat_haps"]
            baf_embryo = data[chrom]["baf_embryo"]
            m_est = MosaicEst(
                mat_haps=mat_haps, pat_haps=pat_haps, bafs=baf_embryo, pos=pos
            )
            # 2a. Use the default parameter settings for mosaic estimation
            m_est.baf_hets()
            sigma_est = np.std(m_est.het_bafs)
            m_est.viterbi_hets(std_dev=sigma_est)
            m_est.create_transition_matrix()
            m_est.est_mle_theta(std_dev=sigma_est)
            if np.isnan(m_est.mle_theta):
                lrt_theta = np.nan
                ci_cf = [np.nan, np.nan, np.nan]
            else:
                ci_theta = m_est.ci_mle_theta(std_dev=sigma_est)
                lrt_theta = m_est.lrt_theta(std_dev=sigma_est)
                ci_cf = [
                    m_est.est_cf(theta=ci_theta[0], gain=gain),
                    m_est.est_cf(theta=ci_theta[1], gain=gain),
                    m_est.est_cf(theta=ci_theta[2], gain=gain),
                ]
            res_dict = {
                "mother": mother,
                "father": father,
                "child": child,
                "chrom": chrom,
                "gain": gain,
                "mosaic_sigma": sigma_est,
                "mle_theta": m_est.mle_theta,
                "lrt_theta": lrt_theta,
                "cellfrac_lower95": ci_cf[0],
                "cellfrac_mean": ci_cf[1],
                "cellfrac_upper95": ci_cf[2],
            }
        else:
            res_dict = {
                "mother": mother,
                "father": father,
                "child": child,
                "chrom": chrom,
                "gain": np.nan,
                "mosaic_sigma": np.nan,
                "mle_theta": np.nan,
                "lrt_theta": np.nan,
                "cellfrac_lower95": np.nan,
                "cellfrac_mean": np.nan,
                "cellfrac_upper95": np.nan,
            }
    except ValueError:
        res_dict = {
            "mother": mother,
            "father": father,
            "child": child,
            "chrom": chrom,
            "gain": np.nan,
            "mosaic_sigma": np.nan,
            "mle_theta": np.nan,
            "lrt_theta": np.nan,
            "cellfrac_lower95": np.nan,
            "cellfrac_mean": np.nan,
            "cellfrac_upper95": np.nan,
        }
    # Convert to a pandas DataFrame
    res_df = pd.DataFrame(res_dict, index=[0])
    res_df.to_csv(snakemake.output["mosaic_tsv"], sep="\t", index=None)
