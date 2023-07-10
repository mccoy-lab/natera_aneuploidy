import numpy as np 
import pickle
import gzip as gz
import click 
import pandas as pd
from karyohmm import MetaHMM

def bph(states):
    """Identify states that are BPH - both parental homologs."""
    idx = []
    for i,s in enumerate(states):
        assert len(s) == 4
        k = 0
        for j in range(4): 
            k += (s[j] >= 0)
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
            k += (s[j] >= 0)
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

def bph_vs_sph_trisomy(gammas, states, pos, start=None, end=None, exclude=False):
    """Classify BPH vs SPH for a trisomy within a given start/end region of a chromosome."""
    assert gammas.shape[1] == pos.size
    assert len(states) == gammas.shape[0]
    if (start is not None) and (end is not None):
        assert start <= end
        # assert start >= np.min(pos)
        # assert end <= np.max(pos)
        if exclude:
            pos_idx = np.where((pos >= end) | (pos <= start))[0]
        else:
            pos_idx = np.where((pos <= end) & (pos >= start))[0]
    else:
        pos_idx = np.arange(pos.size)
    sph_idx = sph(states)
    bph_idx = bph(states)
    posterior_bph = np.sum(np.exp(gammas[bph_idx,:][:, pos_idx]))
    posterior_sph = np.sum(np.exp(gammas[sph_idx,:][:, pos_idx]))
    # scale both of the posteriors appropriately (and compute bayes factors...)
    posterior_bph_norm = posterior_bph / (posterior_bph + posterior_sph)
    posterior_sph_norm = posterior_sph / (posterior_bph + posterior_sph)
    bfs_bph_sph = bayes_factor(np.array([posterior_bph_norm, posterior_sph_norm]))
    return bfs_bph_sph, posterior_bph_norm, posterior_sph_norm


if __name__ == '__main__':
    # Read in the BAF data
    data = pickle.load(gz.open(snakemake.input['baf_pkl'], 'r' ) )
    pos = data[snakemake.wildcards["chrom"]]["pos"]
    # Read in the hmm data + centromere data
    hmm_traceback = pickle.load(gz.open(snakemake.input['hmm_traceback'], 'r' ))
    centromere_df = pd.read_csv(snakemake.input["centromere_bed"], header=None, sep="\s+")
    centromere_df.columns = ["chrom", "start", "end", "feature"]
    bp_padding = snakemake.params["bp_padding"]
    min_centromere_pos = np.nanmin(centromere_df[centromere_df.chrom == snakemake.wildcards["chrom"]].start.values)
    max_centromere_pos = np.nanmax(centromere_df[centromere_df.chrom == snakemake.wildcards["chrom"]].end.values)
    # Obtain the appropriate values from the hmm results
    # NOTE: gammas are assumed to be in log-space here
    gammas = hmm_traceback[snakemake.wildcards['chrom']]['gammas']
    states = MetaHMM().states
    # Step 1. Get the centromeric regions associated appropriately here
    bfs_bph_sph_centromere, posterior_bph_centromere, posterior_sph_centromere = bph_vs_sph_trisomy(gammas, states, pos=pos, start=min_centromere_pos-bp_padding, end=min_centromere_pos+bp_padding)
    # Step 2. Get the centromeric regions associated appropriately here
    bfs_bph_sph_noncentro, posterior_bph_noncentro, posterior_sph_noncentro = bph_vs_sph_trisomy(gammas, states, pos=pos, start=min_centromere_pos-bp_padding, end=min_centromere_pos+bp_padding, exclude=True)
    res_dict = {
        "mother": snakemake.wildcards["mother"], 
        "father": snakemake.wildcards["father"],  
        "child": snakemake.wildcards["child"],  
        "chrom": snakemake.wildcards["chrom"],  
        "bf_bph_centro": bfs_bph_sph_centromere[0],
        "bf_sph_centro": bfs_bph_sph_centromere[1],
        "post_bph_centro": posterior_bph_centromere,
        "post_sph_centro": posterior_sph_centromere,
        "bf_bph_noncentro": bfs_bph_sph_noncentro[0],
        "bf_sph_noncentro": bfs_bph_sph_noncentro[1],
        "post_bph_noncentro": posterior_bph_noncentro,
        "post_sph_noncentro": posterior_sph_noncentro,
        }
    # Convert to a pandas DataFrame
    res_df = pd.DataFrame(res_dict, index=[0])
    res_df.to_csv(snakemake.output["bph_tsv"], sep="\t", index=None)



