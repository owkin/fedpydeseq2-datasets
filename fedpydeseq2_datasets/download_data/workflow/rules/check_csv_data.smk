rule csvize_check_tcga_tumor_purity:
    input:
        "assets/ncomms9971-s2.xlsx",
    output:
        "results/ncomms9971-s2.csv",
        touch("results/ncomms9971-s2/csv.done")
    conda:
        "../envs/python.yaml"
    log:
        "logs/ncomms9971-s2/csv.log"
    script:
        "../scripts/csvize_check_tcga_tumor_purity.py"

rule csvize_check_tcga_cleaned_clinical:
    input:
        "results/1-s2.0-S0092867418302290-mmc1.xlsx",
    output:
        "results/1-s2.0-S0092867418302290-mmc1.csv",
        touch("results/1-s2.0-S0092867418302290-mmc1/csv.done")
    conda:
        "../envs/python.yaml"
    log:
        "logs/1-s2.0-S0092867418302290-mmc1/csv.log"
    script:
        "../scripts/csvize_check_tcga_clinical_data.py"

rule check_recount3_metadata_tcga:
    input:
        "results/{dataset}/metadata.tsv.gz",
    output:
        "results/{dataset}/checked/metadata.tsv.gz",
    log:
        "logs/{dataset}/check_recount3_metadata.log"
    script:
        "../scripts/check_recount3_metadata.py"
