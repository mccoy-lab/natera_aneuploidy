# natera_aneuploidy
GWAS + QuantGen of Aneuploidy in Natera PGT Data

## Description 

The code here is associated with the following [paper](). You can find each specific section of analysis within the `analysis` directory. 

### GWAS

### Phenotyping

### Metadata

### Aneuploidy

The directory `analysis/aneuploidy` contains code for calling whole-chromosome aneuploidies. You can start this section of the analysis using:

```
cd analysis/aneuploidy/
mamba env create -f env.yaml
conda activate natera_aneuploidy_calls
snakemake -s natera_data.smk -j 1 -p
```

which for the first time, will generate all the valid "trios" in the dataset in a file (`valid_trios.txt`). A subsequent run of `snakemake -s natera_data.smk -j 100 -p -n` will create the resulting tables for calling whole-chromosome aneuploidies.

### Simulations

The directory `analysis/simulations` contains code for establishing key benchmarks for the `karyohmm` method for calling aneuploidy. To reproduce the full results run the following: 

```
cd analysis/simulations/
mamba env create -f env.yaml
conda activate natera_aneuploidy_sims
snakemake -s sims.smk -p
```

The primary results files will be deposited as `.tsv` files. 

## Installation

TODO: describe conda environment creation ... 

## Contact 

* Sara Carioscia: @scarioscia
* Arjun Biddanda: @aabiddanda

## External Repositories

* Add in `karyohmm` when online 
