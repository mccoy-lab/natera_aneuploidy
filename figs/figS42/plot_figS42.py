#!python3

"""Plots the suggested thresholds to maintain a calibrated rate of 99% power for aneuploidy calls."""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib

# Read in the data for creating the thresholds for posteriors
sim_df = pd.read_csv('../../data/simulations/fpr_sims_hmm_ploidy.tsv.gz', sep="\t")
ploids = np.unique(sim_df.aploid.values)

alpha = 1 #the lower percentile cutoff

fig, ax = plt.subplots(1, ploids.size,figsize=(ploids.size*2, 2.5), sharey=True)
thresholds = {}
for i,p in enumerate(ploids):
    data = sim_df[sim_df.aploid == p][p].values
    x = np.percentile(data, alpha)
    thresholds[p] = x
    ax[i].hist(data, density=True, bins=50)
    ax[i].axvline(x, color='orange', linestyle='--', label=f'{int(alpha)}st percentile: {x:.2f}')
    ax[i].set_title(p)
    ax[i].set_xlabel(r'Posterior Probability')
    ax[i].legend(fontsize=6, frameon=False)
plt.tight_layout()
thresholds['1p'] = thresholds['1m']
thresholds['3p'] = thresholds['3m']
print(thresholds)
# Create the folders if they are already not created
pathlib.Path("pdfs/").mkdir(parents=True, exist_ok=True)
plt.savefig('pdfs/figS37.pdf', bbox_inches='tight')
