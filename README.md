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
snakemake -s natera_data.smk -j <njobs> -p -n 
```

which for the first time, will generate all the valid "trios" in the dataset in a file (`valid_trios.txt`). A subsequent run of `snakemake -s natera_data.smk -j 100 -p -n` will create the resulting tables for calling whole-chromosome aneuploidies. 

To concatenate the individual tsvs together that are the defined output of the pipeline we use: 

```
cd results/natera_inference/; 
find . -name "*.tsv" | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > natera_embryos.karyohmm_v14.070623.tsv.gz &
```

### Aneuploidy Post-Processing

We also have a post-processing workflow (`aneuploidy_post`) of the posterior tracebacks for defining some additional downstream features. Some of these analyses include:

* Dissection of SPH vs. BPH for trisomies in centromere-proximal and centromere-distal regions of chromosomes for determining trisomies originating at MI vs MII
* Coarse identification of segmental aneuploidies from a set of embryos.



### Simulations

The directory `analysis/simulations` contains code for establishing key benchmarks for the `karyohmm` method for calling aneuploidy. To reproduce the full results run the following: 

```
cd analysis/simulations/
mamba env create -f env.yaml
conda activate natera_aneuploidy_sims
snakemake -s sims.smk -p -j20
```

The primary results files will be deposited as `.tsv.gz` file under the `results` directory. 

### Misc

Miscellaneous analyses for the primary paper. Currently this directory runs: 

* PCA of all Natera parental samples jointly with 2504 1000 Genomes reference individuals. 


## Installation

TODO: describe conda environment creation ... 

## Contact 

* Sara Carioscia: @scarioscia
* Arjun Biddanda: @aabiddanda

## External Repositories

* Add in `karyohmm` when publicly available.
