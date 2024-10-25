"""
author: Margaret R. Starostik, Biology Dept., Johns Hopkins University
email: mstaros1@jhu.edu
last update: October 15, 2024
aim: 
"""

rule predict_expression:
    """
    
    """
    input:
        genotypes = GENOTYPES_DIR + "spectrum_imputed_chr{chromosome}_rehead_filter.cpra.vcf.gz",
        model = WORKING_DIR + "data/GTEx_v8/eqtl/mashr/mashr_{tissue}.db"
    output: 
        prediction = WORKING_DIR + "analysis/twas/predict_expression/{tissue}_chr{chromosome}_predict.txt",
        summary = WORKING_DIR + "analysis/twas/predict_expression/{tissue}_chr{chromosome}_summary.txt"
    log:
        WORKING_DIR + "analysis/twas/predict_expression/{tissue}_chr{chromosome}_log.txt"
    resources:
        partition = "shared",
        mem = "5G",
        time = "1:00:00"
    threads: 4
    shell:
        """
        ml anaconda
        conda activate metaxcan-env 
        python3 /home/mstaros1/tools/MetaXcan/software/Predict.py \
        --model_db_path {input.model} \
        --model_db_snp_key varID \
        --vcf_genotypes {input.genotypes} \
        --vcf_mode imputed \
        --on_the_fly_mapping METADATA "{{}}_{{}}_{{}}_{{}}_b38" \
        --prediction_output {output.prediction} \
        --prediction_summary_output {output.summary} \
        --verbosity 9 \
        --throw 2> {log}
        """        


rule per_tissue_association:
    """
    
    """
    input:
        prediction = WORKING_DIR + "analysis/twas/predict_expression/{tissue}_chr{chromosome}_predict.txt"
    output: 
        association = WORKING_DIR + "analysis/twas/association/single_tissue/{tissue}_chr{chromosome}_association.txt"
    log:
        WORKING_DIR + "analysis/twas/association/single_tissue/logs/{tissue}_chr{chromosome}_log.txt"
    resources:
        partition = "shared",
        mem = "5G",
        time = "2:00:00"
    threads: 4
    shell:
        """
        ml anaconda
        conda activate vscode-env 
        Rscript --vanilla /home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/code/test_association.R {input} {output} 2> {log}
        """ 