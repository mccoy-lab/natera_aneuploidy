import numpy as np 
import pickle
import gzip as gz
import click 
import pandas as pd
from scipy.stats import entropy

def create_entropy(df, categories = ['0', '1m', '1p', '2', '3m', '3p']):
    """Create a column for the posterior entropy."""
    for x in categories:
        assert x in df.columns
    post_vals = df[categories].values
    df['post_entropy'] = entropy(post_vals, axis=1)
    return df

def mean_embryo_noise(df):
    """Define mean embryo noise as a column."""
    assert 'sigma_baf' in df.columns
    assert 'mother' in df.columns
    assert 'father' in df.columns
    assert 'child' in df.columns
    # Aggregate the 
    embryo_noise_df = df.groupby(['mother', 'father', 'child'])[['sigma_baf']].agg(
        np.nanmean
    ).reset_index().rename(columns={'sigma_baf': 'sigma_embryo_mean'})
    df = df.merge(embryo_noise_df, on=['mother', 'father', 'child'], how='left')
    return df

if __name__ == '__main__':
    # Define the core parameters for the filtering steps ...     



