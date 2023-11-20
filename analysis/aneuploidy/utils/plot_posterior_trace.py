"""Diagnostic to generate posterior traces in ploidy from HMM output."""

import numpy as np
import pickle
import click
from tqdm import tqdm
from karyohmm import MetaHMM
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gzip as gz


def plot_posterior_plots(baf_data, hmm_data, chrom="chr1"):
    """Draw a heatmap of the posteriors here."""
    pos = baf_data[chrom]["pos"]
    gammas = hmm_data[chrom]["gammas"]
    xhmm = MetaHMM()
    fig, ax = plt.subplots(1, 1, figsize=(7, 3), layout="tight")
    gamma_karyo = xhmm.marginal_posterior_karyotypes(gammas, xhmm.karyotypes)
    im = ax.imshow(
        gamma_karyo,
        interpolation="none",
        aspect="auto",
        extent=[np.min(pos), np.max(pos), 5, 0],
    )
    ax.set_yticks(range(6))
    ax.set_yticklabels(["0", "1m", "1p", "2", "3m", "3p"])
    ax.set_xlabel("Position (bp)")
    fig.colorbar(im, label="Posterior Probability")
    return fig, ax


@click.command()
@click.option(
    "-b", "--baf", required=True, help="karyoHMM input file (in NPZ numpy format)."
)
@click.option(
    "-h", "--hmm", required=True, help="karyoHMM input file (in NPZ numpy format)."
)
@click.option("-o", "--output", required=True, help="The output PDF file.")
def main(baf, hmm, output):
    chroms = [f"chr{i}" for i in range(1, 23)]
    with gz.open(baf, "r") as f:
        data = np.load(f, allow_pickle=True)
    with gz.open(hmm, "r") as g:
        hmm_data = np.load(g, allow_pickle=True)
    mother = data["mother"]
    father = data["father"]
    child = data["child"]
    with PdfPages(output) as pdf:
        for i, c in tqdm(enumerate(chroms)):
            fig, ax = plot_posterior_plots(baf_data=data, hmm_data=hmm_data, chrom=c)
            if i == 0:
                fig.suptitle(f"{mother}+{father}: {child}\n{c}")
            else:
                fig.suptitle(f"{c}")
            pdf.savefig()
            plt.close()


if __name__ == "__main__":
    main()
