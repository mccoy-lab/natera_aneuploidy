#!python3

"""Plots the suggested thresholds to maintain a calibrated rate of 99% power for aneuploidy calls."""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pathlib


df = pd.read_csv('../../data/simulations/total_hmm_ploidy.tsv.gz', sep="\t")

aploid = np.unique(df.aploid)

fig, ax = plt.subplots(1, aploid.size, figsize=(2*aploid.size,2), sharey=True)
for k,a in enumerate(aploid):
    cur_df = df[df.aploid == a]
    sigmas = np.unique(cur_df.sigma.values)
    pi0 = np.unique(cur_df.pi0.values)
    Z = np.zeros(shape=(sigmas.size, pi0.size))
    for i, s in enumerate(sigmas):
        for j,p in enumerate(pi0):
            Z[i,j] = np.median(cur_df[(cur_df.sigma == s) & (cur_df.pi0 == p)][a].values)
    im = ax[k].imshow(Z, extent=[np.min(sigmas), np.max(sigmas), np.min(pi0), np.max(pi0)], interpolation='lanczos', vmin=0.5, vmax=0.9)
    ax[k].set_xlabel(r'$\hat{\pi}_0$', fontsize=12)
    ax[k].set_title(f'{a}')
ax[0].set_ylabel(r'$\hat{\sigma}$', fontsize=12)
plt.tight_layout()
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.25, 0.01, 0.7])
fig.colorbar(im, cax=cbar_ax, label=r'Median Posterior')

# Create the folders if they are already not created
pathlib.Path("pdfs/").mkdir(parents=True, exist_ok=True)
plt.savefig('pdfs/figS43.pdf', bbox_inches='tight')
