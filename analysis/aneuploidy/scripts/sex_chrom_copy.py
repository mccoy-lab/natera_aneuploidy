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

def loglik_chrx(baf, mat_geno, pat_geno, **kwargs):
    """Obtain the loglikelihood of the embryo BAF conditional on the number. """
    assert baf.ndim == 1
    assert mat_geno.ndim == 2
    assert pat_geno.ndim == 2
    assert baf.size == mat_geno.shape[1]
    assert baf.size == pat_geno.shape[1]
    # Check for heterozygotes in the male X
    if np.any(np.sum(pat_geno, axis=0) == 1):
        warnings.warn("Heterozygotes observed in male chrX ... these are excluded from the overall likelihood.")
    hmm = MetaHMM()
    hmm.states = chrX_states
    hmm.karyotypes = np.array(chrX_state_names, dtype=str)
    # Step 1: first estimate the HMM parameters
    pi0_est, sigma_est = hmm.est_sigma_pi0()
    
    for i, b in enumerate(baf):
        # nullisomy for chrX
        loglik_x0 += emission_baf(b, m=-1, p=-1, k=2, **kwargs)
        # single chrX copyed from 
        loglik_x1 += emission_baf(b, m=mat_geno[0,i], p=-1, k=1, **kwargs)
        # two chrX copied from 
    pass


def loglik_chry(baf, pat_geno, **kwargs):
    """Compute the log-likelihood of the embryo BAF conditional on the number of Y-chromosome copies."""
    assert baf.ndim == 1
    assert pat_geno.ndim == 2
    if np.any(np.sum(pat_geno, axis=0) == 1):
        warnings.warn("Heterozygotes observed on chrY...these are excluded from the likelihood.")
    y0_loglik = 0.0
    y1_loglik = 0.0
    for i, b in enumerate(baf):
        if pat_geno[i,:].sum() != 1:
            # NOTE: this is just the full marginal likelihood. 
            y0_loglik += emission_baf(b, m=-1, p=-1, k=1, **kwargs)
            y1_loglik += emission_baf(b, m=-1, p=pat_geno[0,i], k=1, **kwargs)
    # NOTE: should we normalize here via logsumexp
    return {"y0": y0_loglik, "y1": y1_loglik} 


if __name__ == '__main__':
    """Run the full sex-assignment methods."""
    pass
