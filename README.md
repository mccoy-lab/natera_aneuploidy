# natera_aneuploidy
GWAS + QuantGen of Aneuploidy in Natera PGT Data

## Description 

The code here is associated with the following [paper](https://www.medrxiv.org/content/10.1101/2025.04.02.25325097v1). You can find each  section of analysis in the `analysis` directory. 

### Aneuploidy

The directory `analysis/aneuploidy` contains code for calling whole-chromosome aneuploidies using [karyohmm](https://github.com/mccoy-lab/karyohmm). You can start this section of the analysis using:

```
cd analysis/aneuploidy/
mamba env create -f env.yaml
conda activate natera_aneuploidy_calls
pip install git+https://github.com/mccoy-lab/karyohmm
snakemake -s natera_data.smk -j <njobs> -p -n 
```

which for the first time, will generate all the valid "trios" in the dataset in a file (`valid_trios.txt`). A subsequent run of `snakemake -s natera_data.smk -j 100 -p -n` will create the resulting tables for calling whole-chromosome aneuploidies. 

To concatenate the individual TSV files that are the defined output of the pipeline, we use: 

```
cd results/natera_inference/; 
find . -name "*.tsv" | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > natera_embryos.karyohmm_v20.021024.tsv.gz &
```

### Aneuploidy Post-Processing

We also have a post-processing workflow (`aneuploidy_post`) of the posterior tracebacks for defining some additional downstream features. Some of these analyses include:

* Dissection of SPH vs. BPH for trisomies in centromere-proximal and centromere-distal regions of chromosomes for determining trisomies originating at MI vs MII

* Estimation of mosaic cell fraction for gains + losses that are putatively mitotic in nature (e.g. < 90% posterior probability of being an aneuploidy)

* Estimation of embryo-specific noise parameters for downstream filtering

To run this pipeline - which generates the final version of the autosomal aneuploidy calls: 

```
snakemake -s aneuploidy_post.smk -j 200  ...
```

NOTE: it will take a decent amount of time to build the full DAG for the pipeline. Make sure to have run the steps in the section above as well.  


### GWAS

This directory creates files necessary to run genome-wide association studies (GWAS) for aneuploidy phenotypes. The `gwas.smk` file here creates background files (discovery-test split, parental PCA), phenotype files (including a range of aneuploidy types), and GWAS summary statistic files. 

Specifically, the `scripts/phenotypes` directory within includes a script to create a table for each aneuploidy phenotype, sorted by each parent (e.g., number of embryos affected by maternal triploidy per mother). The `scripts/gwas` file includes a script to run a GWAS for a given phenotype intersected with the genotypes of each parent. 

### Sex embryo 

This directory infers the sex chromosome status for each embryo in the Natera data. First, it preprocesses the BAF data for each trio (mother-father-embryo). It then applies the sex chromosome-ploidy HMM to this BAF data and outputs a TSV with likelihood of each sex chromosome state (XY, XX, X0, XXY, 0, XXX, XXY, Y). 

To run the pipeline (we highly suggest running this on a computing cluster with SLURM): 

```
snakemake -s sex_embryos.smk -j 200 -p ... 
```

and run the following after the pipeline has finished in order to generate the table of sex-chromosome aneuploidy calls.
```
cd results/natera_inference/; 
find . -name "*.tsv" | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > natera_embryos.karyohmm_v20.sex_embryos.021024.tsv.gz &
```

### Simulations

The directory `analysis/simulations` contains code for establishing key benchmarks for the `karyohmm` method for calling aneuploidy. To reproduce the full results run the following: 

```
cd analysis/simulations/
mamba env create -f env.yaml
conda activate natera_aneuploidy_sims
pip install git+https://github.com/mccoy-lab/karyohmm
snakemake -s sims.smk -p -j20
```

The primary results files are  deposited as `.tsv.gz` file under the `results` directory. 

### Misc

Miscellaneous analyses for the primary paper. Currently this directory runs: 

* PCA of all Natera parental samples jointly with 2504 1000 Genomes reference individuals. 

## Contact 

* Sara Carioscia: @scarioscia
* Arjun Biddanda: @aabiddanda

## External Repositories

[karyohmm](https://github.com/mccoy-lab/karyohmm) - Software used in many of the scripts for aneuploidy-calling and generating posterior probability tables for downstream filtering. 

