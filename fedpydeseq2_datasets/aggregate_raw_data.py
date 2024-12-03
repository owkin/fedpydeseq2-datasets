from pathlib import Path

import pandas as pd

from fedpydeseq2_datasets.utils import tnm_to_series

LIST_COL_WITH_UNSPECIFIED_TYPE = [
    "tcga.xml_her2_and_centromere_17_positive_finding_other_measurement_scale_text",
    "tcga.xml_fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17"
    "_signal_result_range",
    "tcga.xml_metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_score",
    "tcga.xml_metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score",
    "tcga.xml_metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result",
    "tcga.xml_metastatic_breast_carcinoma_lab_proc_her2_neu_in_situ_hybridization"
    "_outcome_type",
]


def run_sanity_checks_raw_data(
    dataset_name: str,
    raw_data_path: str | Path,
):
    """
    Run sanity checks on the raw data.

    This function runs sanity checks on the raw data to ensure that the data is
    correctly formatted and that the data is not corrupted.

    It does so by checking the following conditions:
    - The Counts_raw.parquet for the given dataset exists and contains
        between 40 000 and 70 000 genes and between 10 and 1300 samples.
    - The recount3_metadata.tsv.gz file for the given dataset exists and contains
        the columns "external_id" and "tcga.tcga_barcode".
    - The cleaned_clinical_metadata.csv file for the given dataset exists and contains
        the columns "bcr_patient_barcode", "gender" and "ajcc_pathologic_tumor_stage".
    - The tumor_purity_metadata.csv file for the given dataset exists and contains
        the columns "Sample ID" and "CPE".
    - The centers.csv file for the given dataset exists and contains the columns
        "TSS Code" and "Region".

    Parameters
    ----------
    dataset_name : str
        The TCGA dataset name in the format of "tcga-cohort" (capitalized).
        For example, "TCGA-BRCA" for the breast cancer cohort.

    raw_data_path : str or Path
        The path to the raw data folder.
        This raw data folder is assumed to have the following structure and sub files
        <raw_data_path>
        ├── tcga
        │   ├── COHORT
        │   │   ├── Counts_raw.parquet
        │   │   └── recount3_metadata.tsv.gz
        │   ├── cleaned_clinical_metadata.csv
        │   ├── tumor_purity_metadata.csv
        │   └── centers.csv
        Note that the `centers` file is already in the repository. The rest of the
        files can be downloaded with the snakemake pipeline, and are
        already available for the LUAD dataset.

    """
    dataset, cohort = dataset_name.split("-")[:2]
    # Convert to path
    raw_data_path = Path(raw_data_path)
    assert dataset.lower() == "tcga"
    # Check that the centers file exists
    assert (raw_data_path / "tcga" / "centers.csv").exists()
    # Load it
    centers = pd.read_csv(
        raw_data_path / "tcga" / "centers.csv",
    )
    assert "TSS Code" in centers.columns
    assert "Region" in centers.columns
    # Check that the cleaned clinical metadata file exists

    assert (raw_data_path / "tcga" / "cleaned_clinical_metadata.csv").exists()
    # Load it
    cleaned_clinical = pd.read_csv(
        raw_data_path / "tcga" / "cleaned_clinical_metadata.csv"
    )
    # Check it contains the right columns
    assert "bcr_patient_barcode" in cleaned_clinical.columns
    assert "gender" in cleaned_clinical.columns
    assert "ajcc_pathologic_tumor_stage" in cleaned_clinical.columns
    # Check that the tumor purity metadata file exists
    assert (raw_data_path / "tcga" / "tumor_purity_metadata.csv").exists()
    # Load it
    tumor_purity = pd.read_csv(
        raw_data_path / "tcga" / "tumor_purity_metadata.csv",
    )
    assert "Sample ID" in tumor_purity.columns
    assert "CPE" in tumor_purity.columns
    # Check that the recount3 metadata file exists
    assert (raw_data_path / "tcga" / cohort / "recount3_metadata.tsv.gz").exists()
    # Load it

    # specify the dtype of the columns to avoid warnings
    dtype_dict = {col: "object" for col in LIST_COL_WITH_UNSPECIFIED_TYPE}
    recount3_metadata = pd.read_csv(
        raw_data_path / "tcga" / cohort / "recount3_metadata.tsv.gz",
        sep="\t",
        dtype=dtype_dict,
    )
    # Check the columns
    assert "external_id" in recount3_metadata.columns
    assert "tcga.tcga_barcode" in recount3_metadata.columns

    # Check that the counts file exists
    assert (raw_data_path / "tcga" / cohort / "Counts_raw.parquet").exists()
    # Load it
    counts = pd.read_parquet(raw_data_path / "tcga" / cohort / "Counts_raw.parquet")
    # Check that the number of genes is roughly between 40 000 and 70 000
    assert 40000 < counts.shape[0] < 70000
    # Check that the number of samples is greater than 10 and less than 1000
    assert 10 < counts.shape[1] < 1300
    return


def common_preprocessing_tcga(
    dataset_name: str,
    raw_data_path: str | Path,
    processed_data_path: str | Path,
    force: bool = False,
):
    """
    Preprocess the TCGA data and merge all different metadata files.

    This function preprocesses the TCGA data and merges all the different metadata
    files into a single clinical data file.

    It also indexes the count matrix by the barcodes and removes all genes ending with
    PAR_Y. It also removes the gene version by taking the first one.

    Parameters
    ----------
    dataset_name : str
        The TCGA dataset name in the format of "tcga-cohort" (capitalized).
        For example, "TCGA-BRCA" for the breast cancer cohort.

    raw_data_path : str or Path
        The path to the raw data folder.
        This raw data folder is assumed to have the following structure and sub files
        <raw_data_path>
        ├── tcga
        │   ├── COHORT
        │   │   ├── Counts_raw.parquet
        │   │   └── recount3_metadata.tsv.gz
        │   ├── cleaned_clinical_metadata.csv
        │   ├── tumor_purity_metadata.csv
        │   └── centers.csv
        Note that the `centers` file is already in the repository. The rest of the
        files can be downloaded with the snakemake pipeline, and are
        already available for the LUAD dataset.

    processed_data_path : str or Path
        The path to the processed data folder. This function will create the following
        files in this folder:
        <processed_data_path>
        ├── tcga
        │   └── COHORT
        │       ├── counts.parquet
        │       └── clinical_data.csv

    force : bool
        If True, the function will run the preprocessing even if the processed data
        already exists. Default is False.

    Raises
    ------
    ValueError
        If there are missing TSS codes in the centers.csv file.

    """
    if dataset_name in ["TCGA-NSCLC", "TCGA-CRC"]:
        if dataset_name == "TCGA-NSCLC":
            dataset_1, dataset_2 = "TCGA-LUAD", "TCGA-LUSC"
        elif dataset_name == "TCGA-CRC":
            dataset_1, dataset_2 = "TCGA-COAD", "TCGA-READ"

        common_preprocessing_tcga(
            dataset_name=dataset_1,
            raw_data_path=raw_data_path,
            processed_data_path=processed_data_path,
            force=force,
        )
        common_preprocessing_tcga(
            dataset_name=dataset_2,
            raw_data_path=raw_data_path,
            processed_data_path=processed_data_path,
            force=force,
        )
        return

    run_sanity_checks_raw_data(dataset_name, raw_data_path)
    dataset, cohort = dataset_name.split("-")[:2]
    # Convert to path
    raw_data_path = Path(raw_data_path)
    processed_data_path = Path(processed_data_path)
    assert dataset.lower() == "tcga"
    # Load the data from the cohort
    # Check if the outputs exist. If they both do, return
    if (
        (processed_data_path / "tcga" / cohort / "counts.parquet").exists()
        and (processed_data_path / "tcga" / cohort / "clinical_data.csv").exists()
        and not force
    ):
        return

    counts = pd.read_parquet(raw_data_path / "tcga" / cohort / "Counts_raw.parquet")
    # specify the dtype of the columns to avoid warnings
    dtype_dict = {col: "object" for col in LIST_COL_WITH_UNSPECIFIED_TYPE}
    recount3_metadata = pd.read_csv(
        raw_data_path / "tcga" / cohort / "recount3_metadata.tsv.gz",
        sep="\t",
        dtype=dtype_dict,
    )
    # Load only the barcodes and the external ids in the recount3 data

    recount3_metadata = recount3_metadata[
        ["external_id", "tcga.tcga_barcode", "tcga.xml_stage_event_tnm_categories"]
    ]
    # Now create a mapping between the external id and the barcodes
    recount3_metadata = recount3_metadata.set_index("external_id")
    # apply this to the columns
    counts.columns = counts.columns.map(recount3_metadata["tcga.tcga_barcode"])
    # remove the gene_name index
    counts = counts.droplevel("gene_name")
    # Remove all genes ending with PAR_Y
    counts = counts.loc[~counts.index.str.endswith("PAR_Y")]
    counts = counts.T
    counts.index.name = "barcode"
    # Now we filter the gene version by taking the first one.
    counts.columns = counts.columns.str.split(".").str[0]
    counts = counts.loc[:, ~counts.columns.duplicated()]
    counts = counts[~counts.index.duplicated(keep="first")]
    # We now have a clean count matrix
    # we create a column indicating if normal tissue
    is_normal_tissue = (
        recount3_metadata["tcga.tcga_barcode"]
        .map(lambda x: x[:15])
        .map(lambda x: 10 <= int(x[-2:]) <= 29)
    )
    is_normal_tissue.index = recount3_metadata["tcga.tcga_barcode"]
    is_normal_tissue = is_normal_tissue[
        ~is_normal_tissue.index.duplicated(keep="first")
    ]
    # We create T, N, M columns
    tnm_columns = (
        recount3_metadata["tcga.xml_stage_event_tnm_categories"]
        .apply(tnm_to_series)
        .astype("Int8")
    )
    tnm_columns.index = recount3_metadata["tcga.tcga_barcode"]
    tnm_columns = tnm_columns[~tnm_columns.index.duplicated(keep="first")]

    # Now we load the clinical data
    cleaned_clinical = pd.read_csv(
        raw_data_path / "tcga" / "cleaned_clinical_metadata.csv"
    )
    # We filter 3 columns
    cleaned_clinical = cleaned_clinical[
        ["bcr_patient_barcode", "ajcc_pathologic_tumor_stage", "gender"]
    ]
    # We create a correspondance between sample and patient
    cleaned_clinical = cleaned_clinical.set_index("bcr_patient_barcode")

    # Filter the patient that are not in cleaned_clinical
    counts = counts.loc[counts.index.map(lambda x: x[:12] in cleaned_clinical.index)]
    is_normal_tissue = is_normal_tissue.loc[counts.index]
    # create a dataframe

    barcode_to_patient = pd.DataFrame(
        index=counts.index,
        data=counts.index.str[:12].to_list(),
        columns=["bcr_patient_barcode"],
    )

    cleaned_clinical = cleaned_clinical.loc[barcode_to_patient["bcr_patient_barcode"]]
    cleaned_clinical.index = counts.index
    # lower case the gender
    cleaned_clinical.loc[:, "gender"] = cleaned_clinical["gender"].str.lower()

    # Now we define the stage in 1,2,3,4 or NA if it is normal tissue
    def process_stage(stage_str):
        if stage_str.startswith("Stage IV"):
            return 4
        elif stage_str.startswith("Stage III"):
            return 3
        elif stage_str.startswith("Stage II"):
            return 2
        elif stage_str.startswith("Stage I"):
            return 1
        else:
            return pd.NA

    cleaned_clinical.loc[:, "ajcc_pathologic_tumor_stage"] = (
        cleaned_clinical["ajcc_pathologic_tumor_stage"]
        .apply(process_stage)
        .astype("Int16")
    )
    # Add a is_normal_tissue column
    cleaned_clinical.loc[:, "is_normal_tissue"] = is_normal_tissue
    # Set the tumor stage to NA if normal tissue
    cleaned_clinical.loc[is_normal_tissue, "ajcc_pathologic_tumor_stage"] = pd.NA
    # Add the T, N, M columns
    cleaned_clinical = cleaned_clinical.join(tnm_columns)
    # Set T, N, M columns to pd.NA if normal tissue
    cleaned_clinical.loc[is_normal_tissue, ["T", "N", "M"]] = pd.NA
    # now load the tumor purity metadata
    tumor_purity = pd.read_csv(
        raw_data_path / "tcga" / "tumor_purity_metadata.csv",
    )
    barcode_to_sample = pd.DataFrame(
        index=counts.index,
        data=counts.index.str[:16].to_list(),
        columns=["Sample ID"],
    )
    tumor_purity = tumor_purity.set_index("Sample ID")
    # Add nan rows for all sample ids which are in the counts but not in
    # the tumor purity
    missing_sample_ids = list(
        set(barcode_to_sample["Sample ID"]) - set(tumor_purity.index)
    )
    missing_samples = pd.DataFrame(
        index=missing_sample_ids,
        columns=tumor_purity.columns,
    ).astype(tumor_purity.dtypes)
    tumor_purity = pd.concat([tumor_purity, missing_samples], axis=0)

    tumor_purity = tumor_purity.loc[barcode_to_sample["Sample ID"]]
    tumor_purity.index = counts.index
    # We add the CPE column of the tumor purity to the clinical data
    cleaned_clinical.loc[:, "CPE"] = tumor_purity["CPE"]
    # Getting the centers info
    centers_metadata = pd.read_csv(
        raw_data_path / "tcga" / "centers.csv",
    )
    centers_metadata = centers_metadata[["TSS Code", "Region"]]
    centers_metadata["TSS Code"] = centers_metadata["TSS Code"].apply(
        lambda x: x.lstrip("0")
    )
    centers_metadata = centers_metadata.set_index("TSS Code")
    # Get the TSS code
    cleaned_clinical["TSS"] = cleaned_clinical.index.to_series().apply(
        lambda x: x.split("-")[1]
    )
    # Remove leading zeros in both indexes
    cleaned_clinical["TSS"] = cleaned_clinical["TSS"].apply(lambda x: x.lstrip("0"))
    diff = set(cleaned_clinical["TSS"]) - set(centers_metadata.index)
    if len(diff) > 0:
        raise ValueError(
            f"Missing TSS codes {diff} in the centers.csv file."
            f"Please add these missing codes."
        )
    cleaned_clinical = cleaned_clinical.join(centers_metadata, on="TSS")
    cleaned_clinical.drop(columns="TSS", inplace=True)
    # encode the center
    cleaned_clinical["center_id"] = pd.Categorical(cleaned_clinical.Region).codes

    cleaned_clinical.drop(columns="Region", inplace=True)
    # Rename ajcc_pathologic_tumor_stage to stage
    cleaned_clinical.rename(
        columns={"ajcc_pathologic_tumor_stage": "stage"}, inplace=True
    )
    # Now we save the data

    output_folder = processed_data_path / "tcga" / cohort
    output_folder.mkdir(parents=True, exist_ok=True)
    counts.to_parquet(output_folder / "counts.parquet")
    cleaned_clinical.to_csv(output_folder / "clinical_data.csv")
    return
