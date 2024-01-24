import sys

import numpy as np
import pandas as pd
from karyohmm import MosaicEst

if __name__ == "__main__":
    # Read in the input data and params ...
    gain = True
    if int(snakemake.wildcards["p_mono"]) == 0:
        gain = True
    if int(snakemake.wildcards["p_tri"]) == 0:
        gain = False
    # Read in the data ...
    data = np.load(snakemake.input["baf"])
    m_est = MosaicEst(
        mat_haps=data["mat_haps"],
        pat_haps=data["pat_haps"],
        bafs=data["baf_embryo"],
        pos=data["pos"],
    )
    print(
        "Estimating Mosaic Cell Fraction using expected heterozygotes & simulated sigma!"
    )
    m_est.baf_hets()
    m_est.create_transition_matrix()
    sigma = data["std_dev"]
    pi0 = data["mix_prop"]
    seed = data["seed"]
    m_est.est_mle_theta(std_dev=sigma)
    n_hets = m_est.n_het
    ci_theta = m_est.ci_mle_theta(std_dev=sigma)
    sim_cf = np.sum(data["ploidies"] != 2) / data["ploidies"].size
    sigma = data["std_dev"]
    ci_cf = [
        m_est.est_cf(theta=ci_theta[0], gain=False),
        m_est.est_cf(theta=ci_theta[1], gain=False),
        m_est.est_cf(theta=ci_theta[2], gain=gain),
    ]
    exp_het_results = [
        sim_cf,
        ci_cf[0],
        ci_cf[1],
        ci_cf[2],
        m_est.mle_theta,
        sigma,
        sigma,
        pi0,
        n_hets,
        seed,
        "ExpHet",
        gain,
    ]
    print(
        "Estimating Mosaic Cell Fraction using expected heterozygotes & default sigma!"
    )
    m_est.est_mle_theta()
    n_hets = m_est.n_het
    ci_theta = m_est.ci_mle_theta()
    sim_cf = np.sum(data["ploidies"] != 2) / data["ploidies"].size
    sigma = data["std_dev"]
    ci_cf = [
        m_est.est_cf(theta=ci_theta[0], gain=False),
        m_est.est_cf(theta=ci_theta[1], gain=False),
        m_est.est_cf(theta=ci_theta[2], gain=gain),
    ]
    exp_het_results_default = [
        sim_cf,
        ci_cf[0],
        ci_cf[1],
        ci_cf[2],
        m_est.mle_theta,
        0.1,
        sigma,
        pi0,
        n_hets,
        seed,
        "ExpHet",
        gain,
    ]
    print("Estimating Mosaic Cell Fraction using viterbi-genotyping & simulated sigma!")
    m_est.viterbi_hets()
    m_est.create_transition_matrix()
    sigma = data["std_dev"]
    pi0 = data["mix_prop"]
    seed = data["seed"]
    m_est.est_mle_theta(std_dev=sigma)
    n_hets = m_est.n_het
    ci_theta = m_est.ci_mle_theta(std_dev=sigma)
    sim_cf = np.sum(data["ploidies"] != 2) / data["ploidies"].size
    sigma = data["std_dev"]
    ci_cf = [
        m_est.est_cf(theta=ci_theta[0], gain=False),
        m_est.est_cf(theta=ci_theta[1], gain=False),
        m_est.est_cf(theta=ci_theta[2], gain=gain),
    ]
    viterbi_het_results = [
        sim_cf,
        ci_cf[0],
        ci_cf[1],
        ci_cf[2],
        m_est.mle_theta,
        sigma,
        sigma,
        pi0,
        n_hets,
        seed,
        "ViterbiHet",
        gain,
    ]
    print("Estimating Mosaic Cell Fraction using viterbi-genotyping & default sigma!")
    m_est.est_mle_theta()
    n_hets = m_est.n_het
    ci_theta = m_est.ci_mle_theta()
    sim_cf = np.sum(data["ploidies"] != 2) / data["ploidies"].size
    sigma = data["std_dev"]
    ci_cf = [
        m_est.est_cf(theta=ci_theta[0], gain=False),
        m_est.est_cf(theta=ci_theta[1], gain=False),
        m_est.est_cf(theta=ci_theta[2], gain=gain),
    ]
    viterbi_het_results_default = [
        sim_cf,
        ci_cf[0],
        ci_cf[1],
        ci_cf[2],
        m_est.mle_theta,
        0.1,
        sigma,
        pi0,
        n_hets,
        seed,
        "ViterbiHet",
        gain,
    ]
    # Create an updated
    df = pd.DataFrame(
        [
            exp_het_results,
            exp_het_results_default,
            viterbi_het_results,
            viterbi_het_results_default,
        ]
    )
    df.columns = [
        "true_cf",
        "est_cf_lo_95",
        "est_cf_mle",
        "est_cf_hi_95",
        "mle_theta",
        "sigma",
        "sigma0",
        "pi0",
        "n",
        "seed",
        "method",
        "gain",
    ]
    df.to_csv(snakemake.output["mosaic_tsv"], index=None, sep="\t")
