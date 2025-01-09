#!python3

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

if __name__ == "__main__":
    meta_df = pd.read_csv(snakemake.input["natera_metadata"])
    meta_df["IID"] = meta_df.array.values
    pcs = pd.read_csv(snakemake.input["eigenvec"], sep="\t")
    pcs.rename(columns={"#IID": "IID"}, inplace=True)
    kg_phase3_samples = pd.read_csv(snakemake.input["meta_1kg_data"], sep="\t")
    kg_phase3_samples.rename(columns={"KGP_sample_id": "IID"}, inplace=True)
    merged_pc_df = pcs.merge(kg_phase3_samples, how="left")
    kg_phase3_df = merged_pc_df[~merged_pc_df.superpopulation.isna()]
    kg_phase3_pcs = kg_phase3_df[[f"PC{i}" for i in range(1, 21)]].values

    # Filter out 1KG individuals here ...
    natera_df = merged_pc_df[
        merged_pc_df.superpopulation.isna() & ~merged_pc_df.IID.str.contains("NA|HG")
    ]
    natera_pcs = natera_df[[f"PC{i}" for i in range(1, 21)]].values

    # Fit a nearest-neighbor to match the 5 nearest 1KG individuals
    nbrs = NearestNeighbors(n_neighbors=5, algorithm="ball_tree").fit(kg_phase3_pcs)
    distances, indices = nbrs.kneighbors(natera_pcs)

    pop_labels_natera = kg_phase3_df.superpopulation.values[indices]

    pop_labels = []
    for i in tqdm(range(pop_labels_natera.shape[0])):
        cats, counts = np.unique(pop_labels_natera[i, :], return_counts=True)
        pop_labels.append(cats[np.argmax(counts)])

    natera_df["inferred_pop"] = pop_labels
    natera_df.to_csv(snakemake.output["ancestry_table"], sep="\t", index=None)
