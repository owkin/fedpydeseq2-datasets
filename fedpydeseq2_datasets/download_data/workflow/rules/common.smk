from snakemake.utils import validate

validate(config, "../schemas/config.schema.yaml")

TCGA_DATASETS = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
"HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
"PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC",
"UCS", "UVM"]

def get_output(wildcards):
    files=[]
    datasets = config['datasets']
    output_path = config['output_path']
    for dataset in datasets :
        if dataset.upper() in TCGA_DATASETS :
            files += [
                f"{output_path}/tcga/{dataset}/recount3_metadata.tsv.gz",
                f"{output_path}/tcga/{dataset}/Counts_raw.parquet",
            ]
        else:
            raise ValueError(f"Config Error\n\tThe dataset '{dataset}' from the config is neither a TCGA nor a GTEx project.")

    files += [
        f"{output_path}/tcga/tumor_purity_metadata.csv",
        f"{output_path}/tcga/cleaned_clinical_metadata.csv",
        f"{output_path}/tcga/centers.csv"
    ]
    return files
