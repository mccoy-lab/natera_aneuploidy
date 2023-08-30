import numpy as np 
import pandas as pd
import warnings
from cyvcf2 import VCF
from karyohmm import MetaHMM
from .preprocess_natera import obtain_parental_genotypes, obtain_child_data, valid_allele, complement
from scipy.special import logsumexp

def filter_parent_child_sex_chrom_data(child_df, mat_haps, pat_haps, rsids, pos, ref, alt):
    """Filter the resultant parent-child data for this chromosome."""
    assert rsids.size == ref.size
    assert ref.size == alt.size
    assert pos.size == alt.size
    assert (mat_haps.ndim == 2) and (pat_haps.ndim == 2)
    assert mat_haps.shape[1] == rsids.size
    assert pat_haps.shape[1] == rsids.size
    bafs = np.zeros(len(rsids))
    rsid_dict = {}
    for r, baf, B in tqdm(zip(child_df.rsid.values,  child_df.b.values, child_df.B.values)):
        rsid_dict[r] = (baf, B)
    for i, (r, rx, ax) in tqdm(enumerate(zip(rsids, ref, alt))):
        (cur_baf, b_allele) = rsid_dict[r]
        if (
            valid_allele(b_allele)
            and (np.sum(mat_haps[:, i]) in [0, 1, 2])
            and (np.sum(pat_haps[:, i]) in [0, 1, 2])
        ):
            if (b_allele == ax) | (b_allele == complement(ax)):
                bafs[i] = cur_baf
            elif (b_allele == rx) | (b_allele == complement(rx)):
                bafs[i] = 1.0 - cur_baf
            else:
                bafs[i] = np.nan
        else:
            bafs[i] = np.nan
    idx = ~np.isnan(bafs)
    bafs = bafs[idx]
    mat_haps = mat_haps[:, idx]
    pat_haps = pat_haps[:, idx]
    pos = pos[idx]
    ref = ref[idx]
    alt = alt[idx]
    rsids = rsids[idx]
    return bafs, mat_haps, pat_haps, rsids, pos, ref, alt


# Defining a custom set of states for chrX for copying probabilities to run through karyohmm
chrX_states = [
    (-1,-1,-1,-1),
    (-1, -1, 0, -1),
    (0, -1, -1, -1),
    (1, -1, -1, -1)
    (0, -1, 0, -1),
    (1, -1, 0, -1),
    (0, 0, -1, -1),
    (1, 1, -1, -1),
    (0, 1, -1, -1),
    (0, 0, 0, -1),
    (1, 1, 0, -1),
    (0, 1, 0, -1)
    ]
chrX_state_names = ["0", "1p", "1m", "1m", "2", "2", "2m", "2m", "2m", "3", "3", "3"]

chrY_states = [
    (-1,-1,-1,-1)
    (-1,-1,0,-1),
]
chrY_state_names = ["0", "1"]

def posterior_chrX_karyotype(bafs, mat_haps, pat_haps, **kwargs):
    """Compute the posterior karyotype """
    assert bafs.ndim == 1
    assert mat_haps.ndim == 2
    assert pat_haps.ndim == 2
    assert baf.size == mat_haps.shape[1]
    assert baf.size == pat_haps.shape[1]
    # Check for heterozygotes in the male X
    if np.any(np.sum(pat_haps, axis=0) == 1):
        warnings.warn("Heterozygotes observed in male chrX ... these are excluded from the overall likelihood.")
    hmm = MetaHMM()
    # NOTE: we just make the states since karyohmm isn't tailored for this just yet.
    hmm.states = chrX_states
    hmm.karyotypes = np.array(chrX_state_names, dtype=str)
    # Step 1: first estimate the HMM parameters
    pi0_est, sigma_est = hmm.est_sigma_pi0(bafs, mat_haps, pat_haps, r=1e-4)
    # Step 2: run fwd-bwd decoding for the HMM
    gammas, states, karyotypes = hmm.forward_backward(
        bafs=bafs,
        mat_haps=mat_haps,
        pat_haps=pat_haps,
        pi0=pi0_est,
        std_dev=sigma_est,
    )
    res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
    res_dict["pi0_est_chrX"] = pi0_est
    res_dict["sigma_est_chrX"] = sigma_est
    return res_dict

def loglik_chry(baf, pat_haps, **kwargs):
    """Compute the log-likelihood of the embryo BAF conditional on the number of Y-chromosome copies."""
    assert bafs.ndim == 1
    assert pat_haps.ndim == 2
    assert bafs.size == pat_haps.shape[1]
    if np.any(np.sum(pat_haps, axis=0) == 1):
        warnings.warn("Heterozygotes observed on chrY...these are excluded from the likelihood.")
    # NOTE: just create a null maternal haplotype (which we never sample anyways ... 
    mat_haps = np.zeros(shape=(2, bafs.size))
    # 1. Setup HMM with altered state-space 
    hmm = MetaHMM()
    hmm.states = chrY_states
    hmm.karyotypes = chrY_state_names
    pi0_est, sigma_est = hmm.est_sigma_pi0(bafs, mat_haps, pat_haps, r=1e-4)
    # 2. run fwd-bwd decoding for the HMM for state probabilities
    gammas, states, karyotypes = hmm.forward_backward(
        bafs=bafs,
        mat_haps=mat_haps,
        pat_haps=pat_haps,
        pi0=pi0_est,
        std_dev=sigma_est,
    )
    # 3. Aggregate probabilities across karyotypes 
    res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
    res_dict["pi0_est_chrY"] = pi0_est
    res_dict["sigma_est_chrY"] = sigma_est
    return res_dict


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

if __name__ == '__main__':
    """Run the full sex-assignment method."""
    #1. Read in the chrX and chrY files
    pass




