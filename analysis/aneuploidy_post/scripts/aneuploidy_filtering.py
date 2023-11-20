import numpy as np
import pickle
import gzip as gz
import click
import pandas as pd
from scipy.stats import entropy
from sklearn.mixture import GaussianMixture


def create_entropy(df, categories=["0", "1m", "1p", "2", "3m", "3p"]):
    """Create a column for the posterior entropy in the overall aneuploidy data frame."""
    for x in categories:
        assert x in df.columns
    post_vals = df[categories].values
    df["post_entropy"] = entropy(post_vals, axis=1)
    return df


def max_posterior(df, categories=["0", "1m", "1p", "2", "3m", "3p"]):
    """Create a column for the maximum-posterior in the overall aneuploidy data frame."""
    for x in categories:
        assert x in df.columns
    post_vals = df[categories].values
    df["post_max"] = np.max(post_vals, axis=1)
    return df


def embryo_mean_noise(df):
    """Define mean embryo noise as a column in the aneuploidy dataframe."""
    assert "sigma_baf" in df.columns
    assert "mother" in df.columns
    assert "father" in df.columns
    assert "child" in df.columns
    # Aggregate the embryo noise
    embryo_noise_df = (
        df.groupby(["mother", "father", "child"])[["sigma_baf"]]
        .agg(np.nanmean)
        .reset_index()
        .rename(columns={"sigma_baf": "sigma_embryo_mean"})
    )
    df = df.merge(embryo_noise_df, on=["mother", "father", "child"], how="left")
    return df


def mark_day3(df, meta_df):
    """Mark all potential day3 embryos based on metadata criteria."""
    assert "sample_scale" in meta_df.columns
    assert "family_position" in meta_df.columns
    assert "array" in meta_df.columns
    day3_ids = meta_df[meta_df.sample_scale == "single_cell"]["array"].values
    df["day3_embryo"] = df.child.isin(day3_ids)
    return df


def cluster_mosaics(
    df,
    sd_filter=3,
    k=3,
    q=5,
    p_thresh=0.95,
    means_init=[[1.0], [0.75], [0.5]],
    names=["meiotic", "low-mosaic", "high-mosaic"],
):
    """Determine clustering on the mean-embryo noise vs. delta in posterior."""
    assert f"embryo_noise_{sd_filter}sd" in df.columns
    for a in [
        "2",
        "mother",
        "father",
        "child",
        "chrom",
        "sigma_embryo_mean",
        "post_max",
    ]:
        assert a in df.columns
    assert k >= 2
    assert q >= 2
    assert len(names) == k
    assert len(means_init) == k
    gmm_df = df[~df[f"embryo_noise_{sd_filter}sd"] & (df["2"] < p_thresh)][
        ["mother", "father", "child", "chrom", "sigma_embryo_mean", "post_max"]
    ]
    X = gmm_df[["sigma_embryo_mean", "post_max"]].values
    # Step 1: Create the quantiles that each entry belongs to
    qs = np.quantile(X[:, 0], np.linspace(0.0, 1.0, q))
    qs = np.insert(qs, 0, 0)
    qs = np.append(qs, 1.0)
    idx = np.digitize(X[:, 0], bins=qs)
    # Step 2: Run a GMM on each 'tile of embryo noise
    full_labels = np.zeros(idx.size)
    uniq_idx = np.unique(idx)
    gmm = GaussianMixture(n_components=k, covariance_type="tied", means_init=means_init)
    for i in uniq_idx:
        labels = gmm.fit_predict(X[idx == i, 1].reshape(-1, 1))
        full_labels[idx == i] = labels
    # Step 3: Create output that creates new named columns for each proposed category
    for l, n in zip(range(k), names):
        cur_x = np.zeros(X.shape[0])
        cur_x[full_labels == l] = 1
        gmm_df[n] = cur_x
    merged_df = df.merge(gmm_df, on=["mother", "father", "child", "chrom"], how="left")
    return merged_df


if __name__ == "__main__":
    # Define the core parameters for the filtering steps ...
    sd_filter = int(snakemake.params["sd"])
    n_clusters = int(snakemake.params["k"])
    quantiles = int(snakemake.params["q"])
    # 1. Read in the aneuploidy table as a dataframe
    aneuploidy_df = pd.read_csv(snakemake.input["aneuploidy_tsv"], sep="\t")
    meta_df = pd.read_csv(snakemake.input["meta_csv"])
    # 2a. Create the entropy column
    aneuploidy_df = create_entropy(aneuploidy_df)
    # 2b. Create the max posterior column
    aneuploidy_df = max_posterior(aneuploidy_df)
    # 3. Create the mean embryo-noise column
    aneuploidy_df = embryo_mean_noise(aneuploidy_df)
    # 4. Mark day3 embryos (not included in quantile calcs for noise filtering)
    aneuploidy_df = mark_day3(aneuploidy_df, meta_df)
    # 5. Define quantile for upper threshold & column indicating pass on embryo noise (day 5 embryos only!)
    sigma_embryo_means_day5 = (
        aneuploidy_df[~aneuploidy_df.day3_embryo]
        .groupby(["mother", "father", "child"])["sigma_embryo_mean"]
        .head(1)
        .reset_index()["sigma_embryo_mean"]
        .values
    )
    mean_noise = np.nanmean(sigma_embryo_means_day5)
    sd_noise = np.nanstd(sigma_embryo_means_day5)
    aneuploidy_df[f"embryo_noise_{sd_filter}sd"] = (
        aneuploidy_df.sigma_embryo_mean.values >= mean_noise + sd_filter * sd_noise
    ) | (aneuploidy_df.sigma_embryo_mean.values <= mean_noise - sd_filter * sd_noise)
    # 6. Determining the mosaic clusters (in tranches of sigma per-embryo)
    aneuploidy_df = cluster_mosaics(
        aneuploidy_df, sd_filter=sd_filter, k=n_clusters, q=10
    )
    # 7. Write the output to a gzipped TSV
    aneuploidy_df.to_csv(
        snakemake.output["filt_aneuploidy_tsv"], index=None, sep="\t", na_rep="NA"
    )
