import numpy as np 
import pickle
import gzip as gz
import click 
import pandas as pd


def bph(states):
    """Identify states that are BPH - both parental homologs."""
    idx = []
    for i,s in enumerate(states):
        assert len(s) == 4
        k = 0
        for j in range(4): 
            k += (s[i] >= 0)
        if k == 3:
            if s[1] != -1:
                if s[0] != s[1]:
                    # Both maternal homologs present 
                    idx.append(i)
            if s[3] != -1:
                if s[2] != s[3]:
                    # Both paternal homologs present 
                    idx.append(i)
    # Returns indices of both maternal & paternal BPH
    return idx


def sph(states):
    """Identify states that are SPH - single parental homolog."""
    idx = []
    for i,s in enumerate(states):
        assert len(s) == 4
        k = 0
        for j in range(4): 
            k += (s[i] >= 0)
        if k == 3:
            if s[1] != -1:
                if s[0] == s[1]:
                    # Both maternal homologs present 
                    idx.append(i)
            if s[3] != -1:
                if s[2] == s[3]:
                    # Both paternal homologs present 
                    idx.append(i)
    # Returns indices of both maternal & paternal SPH
    return idx

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



def bph_vs_sph_trisomy(gammas, states, pos, start=None, end=None):
    """Classify BPH vs SPH for a trisomy within a given start/end region of a chromosome."""
    assert gammas.shape[1] == pos.size
    assert len(states) == gammas.shape[0]
    if (start is not None) and (end is not None):
        assert start <= end
        assert start >= np.min(pos)
        assert end <= np.max(pos)
        pos_idx = np.where((pos <= end) & (pos >= start))[0]
    else:
        pos_idx = np.arange(pos.size)
    sph_idx = sph(states)
    bph_idx = bph(states)
    posterior_bph = np.sum(np.exp(gammas[bph_idx, pos_idx]))
    posterior_sph = np.sum(np.exp(gammas[sph_idx, pos_idx]))
    # scale both of the posteriors appropriately (and compute bayes factors...)
    posterior_bph_norm = posterior_bph / (posterior_bph + posterior_sph)
    posterior_sph_norm = posterior_sph / (posterior_bph + posterior_sph)
    bfs_bph_sph = bayes_factor(np.array([posterior_bph_norm, posterior_sph_norm]))
    return bfs_bph_sph, posterior_bph_norm, posterior_sph_norm


if __name__ == '__main__':
    pass

