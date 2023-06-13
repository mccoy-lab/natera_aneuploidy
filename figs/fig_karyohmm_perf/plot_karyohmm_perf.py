import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from tqdm import tqdm

def bayes_factor(posteriors, priors=None):
    """Estimate the Bayes Factors """
    if priors is None:
        priors = np.ones(posteriors.size)/posteriors.size
    
    assert posteriors.size == priors.size
    assert np.isclose(np.sum(priors), 1.0)
    bfs = np.zeros(posteriors.size)
    for i in range(posteriors.size):
        denom = np.sum([posteriors[j]*priors[i] for j in range(posteriors.size) if j != i])
        bfs[i] = posteriors[i]*(1 - priors[i]) / denom
    return bfs
    
if __name__ == '__main__':
    sim_performance_df = pd.read_csv('../../data/simulations/total_hmm_ploidy.tsv.gz', sep="\t")
    filt_df = sim_performance_df[(~sim_performance_df.lrr) & (sim_performance_df.a == 0.3)]
    posterior_vals = filt_df[['0', '1m', '1p', '2', '3m', '3p']].values
    bayes_factors_sim = np.array([bayes_factor(posterior_vals[i]) for i in range(posterior_vals.shape[0])])
    
    filt_df['bf_0'] = bayes_factors_sim[:,0]
    filt_df['bf_1m'] = bayes_factors_sim[:,1]
    filt_df['bf_1p'] = bayes_factors_sim[:,2]
    filt_df['bf_2'] = bayes_factors_sim[:,3]
    filt_df['bf_3m'] = bayes_factors_sim[:,4]
    filt_df['bf_3p'] = bayes_factors_sim[:,5]
    filt_df['bf_max'] = filt_df[['bf_0','bf_1m', 'bf_1p', 'bf_2', 'bf_3m', 'bf_3p']].idxmax(axis=1)
    filt_df['bf_max_val'] = filt_df[['bf_0','bf_1m', 'bf_1p', 'bf_2', 'bf_3m', 'bf_3p']].max(axis=1)
    filt_df['post_cat_max'] = filt_df[['0','1m', '1p', '2', '3m', '3p']].idxmax(axis=1)
    filt_df['post_max'] = filt_df[['0','1m', '1p', '2', '3m', '3p']].max(axis=1)
    ploidy_vals = np.unique(filt_df.aploid.values)

    # Creating the plot here ...  
    fig, ax = plt.subplots(1, ploidy_vals.size, figsize=(6,2), sharey=True, sharex=True)
    for i, a in enumerate(ploidy_vals):
        sub_df = filt_df[(filt_df.aploid == a)]
        sigmas = np.unique(sub_df.sigma.values)
        pi0s = np.unique(sub_df.pi0.values)
        for p in pi0s:
            tprs = []
            for s in sigmas:
                true_pos_df = sub_df[(sub_df.sigma == s) & (sub_df.pi0 == p) & (sub_df.post_cat_max == a) & (sub_df.bf_max_val >= 5)]
                tprs.append(true_pos_df.shape[0] / sub_df[(sub_df.sigma == s) & (sub_df.pi0 == p)].shape[0])
            ax[i].scatter(sigmas, tprs, alpha=0.8, marker='+', label=f'$\pi_0$ = {p}')
        ax[i].set_xlabel(r'$\sigma$', fontsize=12)
        ax[i].set_title(f'{a}', fontsize=12)
        ax[i].axhline(0.8, linestyle='--', color='black')
    
    ax[0].legend(frameon=False, fontsize=8)        
    ax[0].set_ylabel(r'True Positive Rate', fontsize=12)
    plt.tight_layout()
    plt.savefig('karyohmm_perf.pdf', bbox_inches='tight')
    plt.savefig('karyohmm_perf.png', bbox_inches='tight', dpi=300)
    