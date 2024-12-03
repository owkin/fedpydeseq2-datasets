rule parquetize_raw_counts:
    input:
        "results/{dataset}/Counts_raw.tsv.gz",
        "results/{dataset}/gene_names.tsv.gz"
    output:
        "results/{dataset}/Counts_raw.parquet",
        touch("results/{dataset}/raw_parquet.done")
    conda:
        "../envs/python.yaml"
    log:
        "logs/{dataset}/raw_parquet.log"
    script:
        "../scripts/parquetize_check_data.py"
