import numpy as np
import pandas as pd
import warnings
import gzip as gz
from karyohmm import MetaHMM

# Defining a custom set of states for chrX for copying probabilities to run through karyohmm
chrX_states = [
    (-1, -1, -1, -1),
    (-1, -1, 0, -1),
    (0, -1, -1, -1),
    (1, -1, -1, -1),
    (0, -1, 0, -1),
    (1, -1, 0, -1),
    (0, 0, -1, -1),
    (1, 1, -1, -1),
    (0, 1, -1, -1),
    (0, 0, 0, -1),
    (1, 1, 0, -1),
    (0, 1, 0, -1),
]
chrX_state_names = [
    "x0",
    "x1p",
    "x1m",
    "x1m",
    "x2",
    "x2",
    "x2m",
    "x2m",
    "x2m",
    "x3",
    "x3",
    "x3",
]

chrY_states = [
    (-1, -1, -1, -1),
    (-1, -1, 0, -1),
]
chrY_state_names = ["y0", "y1"]


def bayes_factor(posteriors, priors=None):
    """Compute Bayes Factors for evidence of specific aneuploidy states."""
    if priors is None:
        priors = np.ones(posteriors.size) / posteriors.size

    assert posteriors.size == priors.size
    assert np.isclose(np.sum(priors), 1.0)
    bfs = np.zeros(posteriors.size)
    for i in range(posteriors.size):
        denom = np.sum(
            [posteriors[j] * priors[i] for j in range(posteriors.size) if j != i]
        )
        bfs[i] = posteriors[i] * (1 - priors[i]) / denom
    return bfs


def posterior_chrX_karyotype(bafs, mat_haps, pat_haps, **kwargs):
    """Compute the posterior karyotype"""
    assert bafs.ndim == 1
    assert mat_haps.ndim == 2
    assert pat_haps.ndim == 2
    assert bafs.size == mat_haps.shape[1]
    assert bafs.size == pat_haps.shape[1]
    # Check for heterozygotes in the male X
    if np.any(np.sum(pat_haps, axis=0) == 1):
        warnings.warn(
            "Heterozygotes observed in male chrX ... these are excluded from the overall likelihood."
        )
    # 1. Create the X-specific HMM
    # NOTE: we just make the states since karyohmm isn't tailored for this just yet.
    hmm = MetaHMM()
    hmm.states = chrX_states
    hmm.karyotypes = np.array(chrX_state_names, dtype=str)
    # 2: first estimate the HMM parameters
    pi0_est, sigma_est = hmm.est_sigma_pi0(bafs, mat_haps, pat_haps, r=1e-4)
    # 3: run fwd-bwd decoding for the HMM
    gammas, states, karyotypes = hmm.forward_backward(
        bafs=bafs,
        mat_haps=mat_haps,
        pat_haps=pat_haps,
        pi0=pi0_est,
        std_dev=sigma_est,
        unphased=True,
    )
    res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
    res_dict["pi0_est_chrX"] = pi0_est
    res_dict["sigma_est_chrX"] = sigma_est
    # 4. Obtain the maxBF and category for maxBF on chrX
    posteriors = np.array([res_dict[x] for x in karyotypes])
    bfs = bayes_factor(posteriors)
    res_dict["x_maxBF"] = np.max(bfs)
    res_dict["x_maxBF_cat"] = karyotypes[np.argmax(bfs)]
    return res_dict


def posterior_chrY_karyotype(bafs, pat_haps, **kwargs):
    """Compute the log-likelihood of the embryo BAF conditional on the number of Y-chromosome copies."""
    assert bafs.ndim == 1
    assert pat_haps.ndim == 2
    assert bafs.size == pat_haps.shape[1]
    if np.any(np.sum(pat_haps, axis=0) == 1):
        warnings.warn(
            "Heterozygotes observed on chrY...these are excluded from the likelihood."
        )
    # NOTE: just create a null maternal haplotype (which we never sample anyways ...)
    mat_haps = np.zeros(shape=(2, bafs.size))
    # 1. Setup HMM with altered state-space
    hmm = MetaHMM()
    hmm.states = chrY_states
    hmm.karyotypes = np.array(chrY_state_names, dtype=str)
    pi0_est, sigma_est = hmm.est_sigma_pi0(
        bafs, mat_haps, pat_haps, sigma_bounds=(1e-2, 0.4), r=1e-4
    )
    # 2. run fwd-bwd decoding for the HMM for state probabilities
    gammas, states, karyotypes = hmm.forward_backward(
        bafs=bafs,
        mat_haps=mat_haps,
        pat_haps=pat_haps,
        pi0=pi0_est,
        std_dev=sigma_est,
        unphased=True,
    )
    # 3. Aggregate probabilities across karyotypes
    res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
    res_dict["pi0_est_chrY"] = pi0_est
    res_dict["sigma_est_chrY"] = sigma_est
    # 4. Obtain the maxBF and category for maxBF for chrY
    posteriors = np.array([res_dict[x] for x in karyotypes])
    bfs = bayes_factor(posteriors)
    res_dict["y_maxBF"] = np.max(bfs)
    res_dict["y_maxBF_cat"] = karyotypes[np.argmax(bfs)]
    return res_dict


if __name__ == "__main__":
    """Run the full sex-assignment method and output a TSV."""
    # 1. Read in the BAF data for analysis
    data = pickle.load(gz.open(snakemake.input["baf_pkl"], "r"))
    assert "chrX" in data.keys()
    assert "chrY" in data.keys()
    baf_chrX_data = data["chrX"]
    baf_chrY_data = data["chrY"]

    # 2. Compute the posteriors/BF for chrX
    n01_chrX = np.nansum(
        (baf_chrX_data["baf_embryo"] == 1) | (baf_chrX_data["baf_embryo"] == 0)
    )
    m_chrX = baf_chrX_data["baf_embryo"].size
    if n01_chrX == m_chrX:
        print("Warning: all BAF values were either [0,1].", file=sys.stderr)
        res_dict_chrX = {}
        for k in chrX_state_names:
            res_dict_chrX[k] = np.nan
        res_dict_chrX["pi0_est_chrX"] = np.nan
        res_dict_chrX["sigma_est_chrX"] = np.nan
        res_dict_chrX["x_maxBF"] = np.nan
        res_dict_chrX["x_maxBF_cat"] = np.nan
    else:
        res_dict_chrX = posterior_chrX_karyotype(
            bafs=baf_chrX_data["baf_embryo"],
            mat_haps=baf_chrX_data["mat_haps"],
            pat_haps=baf_chrX_data["pat_haps"],
        )
    # 3. Compute the posteriors/BF for chrY
    n01_chrY = np.nansum(
        (baf_chrY_data["baf_embryo"] == 1) | (baf_chrY_data["baf_embryo"] == 0)
    )
    m_chrY = baf_chrY_data["baf_embryo"].size
    if n01_chrY == m_chrY:
        print("Warning: all BAF values were either [0,1].", file=sys.stderr)
        res_dict_chrY = {}
        for k in chrY_state_names:
            res_dict_chrY[k] = np.nan
        res_dict_chrY["pi0_est_chrY"] = np.nan
        res_dict_chrY["sigma_est_chrY"] = np.nan
        res_dict_chrY["y_maxBF"] = np.nan
        res_dict_chrY["y_maxBF_cat"] = np.nan
    else:
        res_dict_chrY = posterior_chrY_karyotype(
            bafs=baf_chrY_data["baf_embryo"], pat_haps=baf_chrY_data["pat_haps"]
        )
    # 4. Create the line for the primary output
    with open(snakemake.output["karyo_tsv"], "w+") as out:
        out.write(
            "mother\tfather\tchild\tpi0_x\tsigma_x\tx0\tx1p\tx1m\tx2\tx2m\tx3\tx_maxBF\tx_maxBFcat\tpi0_y\tsigma_y\ty0\ty1\ty_maxBF\ty_maxBFcat\n"
        )
        # NOTE: this is kind of long and gross - and can be probably made more clear in a separate function?
        line = f"{baf_chrX_data['mother_id']}\t{baf_chrX_data['father_id']}\t{baf_chrX_data['child_id']}\t{res_dict_chrX['pi0_est_chrX']}\t{res_dict_chrX['sigma_est_chrX']}\t{res_dict_chrX['x0']}\t{res_dict_chrX['x1p']}\t{res_dict_chrX['x1m']}\t{res_dict_chrX['x2']}\t{res_dict_chrX['x2m']}\t{res_dict_chrX['x3']}\t{res_dict_chrX['x_maxBF']}\t{res_dict_chrX['x_maxBF_cat']}\t{res_dict_chrY['pi0_est_chrY']}\t{res_dict_chrY['sigma_est_chrY']}\t{res_dict_chrY['y0']}\t{res_dict_chrY['y1']}\t{res_dict_chrY['y_maxBF']}\t{res_dict_chrY['y_maxBF_cat']}\n"
        out.write(line)
