rule download_recount3_data:
    output:
        "results/{dataset}/Counts_raw.tsv.gz",
        "results/{dataset}/metadata.tsv.gz",
        "results/{dataset}/gene_names.tsv.gz",
    retries: 3
    resources: download_concurrent=1
    conda:
        "../envs/recount3.yaml"
    log:
        "logs/{dataset}/download.log",
    script:
        "../scripts/download_recount3_cohort.R"

rule download_cleaned_tcga_clinical_data:
    output:
        "results/1-s2.0-S0092867418302290-mmc1.xlsx",
    retries: 3
    shell:
        "wget -P results/ https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx"
