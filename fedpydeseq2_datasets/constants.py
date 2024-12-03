from typing import Literal
from typing import cast

TCGADatasetNames = Literal[
    "TCGA-LUAD",
    "TCGA-PAAD",
    "TCGA-BRCA",
    "TCGA-COAD",
    "TCGA-LUSC",
    "TCGA-READ",
    "TCGA-SKCM",
    "TCGA-PRAD",
    "TCGA-NSCLC",
    "TCGA-CRC",
]

TCGA_DATASET_NAMES = [
    cast(TCGADatasetNames, dataset)
    for dataset in [
        "TCGA-LUAD",
        "TCGA-PAAD",
        "TCGA-BRCA",
        "TCGA-COAD",
        "TCGA-LUSC",
        "TCGA-READ",
        "TCGA-SKCM",
        "TCGA-PRAD",
    ]
]
