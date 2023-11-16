import numpy as np 
import pickle
import gzip as gz
import click 
import pandas as pd
from scipy.stats import entropy

def create_entropy(df, categories = ['0', '1m', '1p', '2', '3m', '3p']):
    """Create a column for the posterior entropy in the overall aneuploidy data frame."""
    for x in categories:
        assert x in df.columns
    post_vals = df[categories].values
    df['post_entropy'] = entropy(post_vals, axis=1)
    return df

def max_posterior(df, categories = ['0', '1m', '1p', '2', '3m', '3p']):
    """Create a column for the maximum-posterior in the overall aneuploidy data frame."""
    for x in categories:
        assert x in df.columns
    post_vals = df[categories].values
    df['post_max'] = np.max(post_vals, axis=1)
    return df

def embryo_mean_noise(df):
    """Define mean embryo noise as a column in the aneuploidy dataframe."""
    assert 'sigma_baf' in df.columns
    assert 'mother' in df.columns
    assert 'father' in df.columns
    assert 'child' in df.columns
    # Aggregate the embryo noise 
    embryo_noise_df = df.groupby(['mother', 'father', 'child'])[['sigma_baf']].agg(
        np.nanmean
    ).reset_index().rename(columns={'sigma_baf': 'sigma_embryo_mean'})
    df = df.merge(embryo_noise_df, on=['mother', 'father', 'child'], how='left')
    return df

def mark_day3(df, meta_df):
    """Mark all potential day3 embryos based on metadata criteria."""
    assert 'sample_scale' in meta_df.columns
    assert 'family_position' in meta_df.columns
    assert 'array' in meta_df.columns
    day3_ids = meta_df[meta_df.sample_scale == 'single_cell']['array'].values
    df['day3_embryo'] = df.child.isin(day3_ids)
    return df
    

if __name__ == '__main__':
    # Define the core parameters for the filtering steps ...     
    sd_filter = int(snakemake.params['sd'])
    n_clusters = snakemake.params['k']
    # 1. Read in the aneuploidy table as a dataframe
    aneuploidy_df = pd.read_csv(snakemake.input['aneuploidy_tsv'], sep="\t")
    meta_df = pd.read_csv(snakemake.input['meta_csv'])
    # 2a. Create the entropy column 
    aneuploidy_df = create_entropy(aneuploidy_df)
    # 2b. Create the max posterior column
    aneuploidy_df = max_posterior(aneuploidy_df)
    # 3. Create the mean embryo-noise column
    aneuploidy_df = embryo_mean_noise(aneuploidy_df)
    # 4. Mark day3 embryos (not included in quantile calcs for noise filtering) 
    aneuploidy_df = mark_day3(aneuploidy_df, meta_df)
    # 5. Define quantile for upper threshold & column indicating pass on embryo noise (day 5 embryos only!)
    sigma_embryo_means_day5 = aneuploidy_df[~aneuploidy_df.day3_embryo].groupby(
        ['mother', 'father', 'child']
    )['sigma_embryo_mean'].head(1).reset_index()['sigma_embryo_mean'].values
    mean_noise = np.nanmean(sigma_embryo_means_day5)
    sd_noise = np.nanstd(sigma_embryo_means_day5)
    aneuploidy_df[f'embryo_noise_{sd_filter}sd'] =  (aneuploidy_df.sigma_embryo_mean.values >= mean_noise + sd_filter*sd_noise) | (aneuploidy_df.sigma_embryo_mean.values <= mean_noise - sd_filter*sd_noise)
    
