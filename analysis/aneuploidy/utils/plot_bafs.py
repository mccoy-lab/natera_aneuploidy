"""Diagnostic to generate per-genotype BAF histograms per chromosome to diagnose tricky behavior."""

import numpy as np
import pickle
import click
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import gzip as gz


def plot_baf_plots(baf_data, chrom="chr1", full_data=False):
    """Draw BAF vs. Genotype histogram."""
    mat_haps = baf_data[chrom]["mat_haps"]
    pat_haps = baf_data[chrom]["pat_haps"]
    baf_embryo = baf_data[chrom]["baf_embryo"]
    mat_geno = np.sum(mat_haps, axis=0)
    pat_geno = np.sum(pat_haps, axis=0)

    fig, ax = plt.subplots(3, 3, figsize=(6, 5), sharex=True, layout="constrained")
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            idx = np.where((mat_geno == i) & (pat_geno == j))
            ax[i, j].hist(baf_embryo[idx], bins=50)
            ax[i, j].set_title(f"M={i},P={j}")
    for a in ax[2, :]:
        a.set_xlabel(r"BAF")
    if full_data:
        return fig, ax, mat_haps, pat_haps, baf_embryo
    return fig, ax


@click.command()
@click.option(
    "-b", "--baf", required=True, help="BAF input file (in NPZ numpy format)."
)
@click.option("-o", "--output", required=True, help="The output PDF file.")
def main(baf, output):
    chroms = [f"chr{i}" for i in range(1, 23)]
    with gz.open(baf, "r") as f:
        data = np.load(f, allow_pickle=True)
    mother = data["mother"]
    father = data["father"]
    child = data["child"]
    with PdfPages(output) as pdf:
        for i, c in tqdm(enumerate(chroms)):
            fig, ax = plot_baf_plots(baf_data=data, chrom=c)
            if i == 0:
                fig.suptitle(f"{mother}+{father}: {child}\n{c}")
            else:
                fig.suptitle(f"{c}")
            pdf.savefig()
            plt.close()


if __name__ == "__main__":
    main()
