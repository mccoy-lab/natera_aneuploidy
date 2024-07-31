#!python3

import numpy as np
from karyohmm import PGTSim, PGTSimMosaic

if __name__ == "__main__":
    if snakemake.params["sfs"] != "None":
        afs = np.loadtxt(snakemake.params["sfs"])
    else:
        afs = None
    # Run the full simulation using the defined helper function
    if snakemake.params["mixed_ploidy"]:
        p_mono = snakemake.params["p_mono"]
        p_tri = snakemake.params["p_tri"]
        table_data = PGTSimMosaic().mixed_ploidy_sim(
            afs=afs,
            ploidies=np.array([0, 1, 2, 3]),
            props=np.array([0.0, p_mono, 1.0 - p_mono - p_tri, p_tri]),
            ncells=snakemake.params["n"],
            m=snakemake.params["m"],
            length=35e6,
            mat_skew=snakemake.params["mat_skew"],
            std_dev=snakemake.params["sigma"],
            mix_prop=snakemake.params["pi0"],
            alpha=snakemake.params["alpha"],
            seed=snakemake.params["seed"],
        )
    else:
        table_data = PGTSim().full_ploidy_sim(
            afs=afs,
            seed=snakemake.params["seed"],
            ploidy=snakemake.params["k"],
            m=snakemake.params["m"],
            length=35e6,
            std_dev=snakemake.params["sigma"],
            mix_prop=snakemake.params["pi0"],
            mat_skew=snakemake.params["mat_skew"],
            alpha=snakemake.params["alpha"],
            switch_err_rate=3e-2,
        )
    try:
        table_data["mother"] = snakemake.params["father_id"]
        table_data["father"] = snakemake.params["mother_id"]
        table_data["child"] = snakemake.params["child_id"]
    except KeyError:
        pass
    np.savez_compressed(snakemake.output["baf"], **table_data)
