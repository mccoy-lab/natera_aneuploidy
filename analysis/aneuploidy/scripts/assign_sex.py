import numpy as np 
import warnings
from karyohmm_utils import emission_baf

# Define the copy-number for chrX and chrY that are possible 
chrX_copy_num = [0,1,2,3]
chrY_copy_num = [0,1]


def loglik_chrx(baf, mat_geno, pat_geno, **kwargs):
    """Obtain the loglikelihood of the embryo BAF conditional on the number. """
    assert baf.ndim == 1
    assert mat_geno.ndim == 2
    assert pat_geno.ndim == 2
    assert baf.size == mat_geno.shape[1]
    assert baf.size == pat_geno.shape[1]
    pass


def loglik_chry(baf, pat_geno, ncopy=0):
    """Compute the log-likelihood of the embryo BAF conditional on the number of Y-chromosome copies."""
    assert baf.ndim == 1
    assert pat_geno.ndim == 2
    if np.any(np.sum(pat_geno, axis=0) == 1):
        warnings.warn("Heterozygotes observed on chrY...be wary of genotyping noise!")
    loglik = 0.0
    for i, b in enumerate(baf):
        pass


if __name__ == '__main__':
    """Run the full sex-assignment methods."""
    pass
