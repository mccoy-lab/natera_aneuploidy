import numpy as np 
import warnings
from karyohmm import MetaHMM
from scipy.special import logsumexp

# Define the copy-number for chrX and chrY that are possible 
# chrX_copy_num = [0,1,2,3]
# chrY_copy_num = [0,1]

# Defining a custom set of states for chrX for copying probabilities
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

def loglik_chrx(bafs, mat_haps, pat_haps, **kwargs):
    """Obtain the loglikelihood of the embryo BAF conditional on the number. """
    assert bafs.ndim == 1
    assert mat_haps.ndim == 2
    assert pat_haps.ndim == 2
    assert baf.size == mat_haps.shape[1]
    assert baf.size == pat_haps.shape[1]
    # Check for heterozygotes in the male X
    if np.any(np.sum(pat_haps, axis=0) == 1):
        warnings.warn("Heterozygotes observed in male chrX ... these are excluded from the overall likelihood.")
    hmm = MetaHMM()
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
    mat_haps = np.zeros(shape=(2, bafs.size))
    pass
    


if __name__ == '__main__':
    """Run the full sex-assignment methods."""
    pass
