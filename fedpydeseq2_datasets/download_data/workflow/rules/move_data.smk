output_path = config["output_path"]

rule move_recount3_metadata_tcga:
    input:
        metadata  = "results/{dataset}/checked/metadata.tsv.gz",
    output:
        f"{output_path}"+"/tcga/{dataset}/recount3_metadata.tsv.gz",
    shell:
        """
        mv {input.metadata} {output}
        """


rule move_clinical_data_tcga:
    input:
        clinical_metadata = "results/1-s2.0-S0092867418302290-mmc1.csv",
        clinical_metadata_done = "results/1-s2.0-S0092867418302290-mmc1/csv.done",
    output:
        f"{output_path}"+"/tcga/cleaned_clinical_metadata.csv",
    shell:
        """
        mv {input.clinical_metadata} {output}
        """

rule move_tumor_purity_data_tcga:
    input:
        tumor_purity_metadata ="results/ncomms9971-s2.csv",
        tumor_purity_metadata_done = "results/ncomms9971-s2/csv.done",
    output:
        f"{output_path}"+"/tcga/tumor_purity_metadata.csv",
    shell:
        """
        mv {input.tumor_purity_metadata} {output}
        """

rule move_counts_tcga:
    input:
        counts = "results/{dataset}/Counts_raw.parquet",
        parquet_done = "results/{dataset}/raw_parquet.done",
    output:
        f"{output_path}"+"/tcga/{dataset}/Counts_raw.parquet",
    shell:
        """
        mv {input.counts} {output}
        """

rule copy_centers_csv:
    output:
        f"{output_path}"+"/tcga/centers.csv"
    shell:
        """
        cp assets/centers.csv {output}
        """
