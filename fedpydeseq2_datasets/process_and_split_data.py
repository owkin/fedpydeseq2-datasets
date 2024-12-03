"""Generate the centers' TCGA datasets.

The function is used to preprocess and split into centers the tcga dataset.
"""
import shutil
from functools import partial
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
from loguru import logger

from fedpydeseq2_datasets.aggregate_raw_data import common_preprocessing_tcga
from fedpydeseq2_datasets.constants import TCGADatasetNames
from fedpydeseq2_datasets.utils import get_experiment_id
from fedpydeseq2_datasets.utils import mix_centers

ALLOWED_DESIGN_FACTORS_TCGA = {"stage", "gender", "CPE"}
ALLOWED_CONTINUOUS_FACTORS_TCGA = {"CPE"}


def setup_tcga_dataset(
    raw_data_path: str | Path,
    processed_data_path: str | Path,
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
    small_samples: bool = False,
    small_genes: bool = False,
    only_two_centers: bool = False,
    design_factors: str | list[str] = "stage",
    continuous_factors: list[str] | None = None,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
    force: bool = False,
    **pydeseq2_kwargs,
):
    """Load, clean and split a TCGA dataset into clients.

    This function is given the path to the raw data.
    Afterward, the natural split (centers) of the dataset is used
    to separate the data into different dataframes, one for each center.
    Finally the function saves those dataframe on the disk as a csv file for each
    center. Those datasets would correspond to center's dataset in substrafl.

    The file structure expected as an input is the following:

    ```
    <raw_data_path>
    ├── tcga
    │   ├── COHORT
    │   │   ├── Counts_raw.parquet
    │   │   └── recount3_metadata.tsv.gz
    │   ├── centers.csv
    │   ├── tumor_purity_metadata.csv
    │   └── cleaned_clinical_metadata.csv
    └──
    ```

    The processed data path will be built/filled/ checked with the following structure.
    The cohort is the tcga cohort, which we collect from the dataset name.
    The true dataset name is generated from the dataset name, the design factors,
    the continuous factors, the small samples, the small genes and the only two centers
    parameters.

    ```
    <processed_data_path>
    ├── tcga
    │   └── COHORT
    │       ├── counts.parquet
    │       └── clinical_data.csv
    ├── centers_data
    │   └── tcga
    │       └── <experiment_id>
    │           └── center_0
    │               ├── counts_data.csv
    │               └── metadata.csv
    ├── pooled_data
    │   └── tcga
    │       └── <experiment_id>
    │           ├── counts_data.csv
    │           └── metadata.csv
    └──
    ```

    Parameters
    ----------
    raw_data_path : str or Path
        Path to raw data.
    processed_data_path : str or Path
        Path to processed data.
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
        The dataset to preprocess, by default "TCGA-LUAD".
    small_samples : bool
        If True, only preprocess a small subset of the data, by default False.
        This small subset is composed of 10 samples per center (or the number
        of samples in the center if this number is inferior).
    small_genes : bool
        If True, only preprocess a small subset of the data features (genes)
        , by default False. This small subset is composed of 100 genes.
    only_two_centers : bool
        If True, split the data in two centers only, by default False.
    design_factors : str or list
        The design factors.
    continuous_factors : list[str] or None
        The continuous design factors. Factors not in this list will be considered
        as categorical.
    heterogeneity_method : str or None
        The method to use to generate heterogeneity in the dataset. If None, no
        heterogeneity is generated. Default is None.
    heterogeneity_method_param : float or None
        The parameter of the heterogeneity method. Default is None.
    force : bool
        If True, force the setup of the dataset even if the metadata already exists.
        Default is False.
    **pydeseq2_kwargs
        Additional arguments to pass to the pydeseq2 and fedpydeseq2 classes and
        strategies.
    """
    raw_data_path = Path(raw_data_path).resolve()
    processed_data_path = Path(processed_data_path).resolve()

    common_preprocessing_tcga(
        dataset_name=dataset_name,
        raw_data_path=raw_data_path,
        processed_data_path=processed_data_path,
        force=force,
    )

    if isinstance(design_factors, str):
        design_factors = [design_factors]

    experiment_id = get_experiment_id(
        dataset_name,
        small_samples,
        small_genes,
        only_two_centers,
        design_factors,
        continuous_factors,
        heterogeneity_method=heterogeneity_method,
        heterogeneity_method_param=heterogeneity_method_param,
        **pydeseq2_kwargs,
    )

    logger.info(f"Setting up TCGA dataset: {experiment_id}")

    center_data_path = processed_data_path / "centers_data" / "tcga" / experiment_id
    first_center_metadata_path = center_data_path / "center_0" / "metadata.csv"

    if not first_center_metadata_path.exists() or force:
        logger.info(
            "First center metadata does not exist or force=True. Setting up the "
            "dataset."
        )
        return _setup_tcga_dataset(
            processed_data_path,
            dataset_name,
            small_samples,
            small_genes,
            only_two_centers,
            design_factors,
            continuous_factors,
            heterogeneity_method,
            heterogeneity_method_param,
            **pydeseq2_kwargs,
        )
    # Check if the metadata contains all the design factors
    logger.info(
        f"First center metadata exists at {first_center_metadata_path}. "
        f"Checking if all design factors are present."
    )
    metadata = pd.read_csv(first_center_metadata_path, index_col=0)
    for design_factor in design_factors:
        if design_factor not in metadata.columns:
            logger.info(
                f"Design factor {design_factor} not present in the metadata."
                f" Setting up the dataset."
            )
            return _setup_tcga_dataset(
                processed_data_path,
                dataset_name,
                small_samples,
                small_genes,
                only_two_centers,
                design_factors,
                continuous_factors,
                heterogeneity_method,
                heterogeneity_method_param,
                **pydeseq2_kwargs,
            )


def _setup_tcga_dataset(
    processed_data_path: str | Path,
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
    small_samples=False,
    small_genes=False,
    only_two_centers=False,
    design_factors: str | list[str] = "stage",
    continuous_factors: list[str] | None = None,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
    **pydeseq2_kwargs,
):
    """Load, clean and split a TCGA dataset into clients.

    This is the main function to perform these operations.

    This function is given the path to the raw data.
    Afterward, the natural split (centers) of the dataset is used
    to separate the data into different dataframes, one for each center.
    Finally the function saves those dataframe on the disk as a csv
    file for each center.
    Those datasets would correspond to center's dataset in substrafl.

    For the file structure description, see the `setup_tcga_dataset` function.

    Parameters
    ----------
    processed_data_path : str or Path
        This path will contain a folder 'centers_data' and 'pooled_data'
        to store the preprocessed data.
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
        The dataset to preprocess, by default "TCGA-LUAD".
    small_samples : bool
        If True, only preprocess a small subset of the data, by default False.
        This small subset is composed of 10 samples per center.
    small_genes : bool
        If True, only preprocess a small subset of the data features (genes)
        , by default False. This small subset is composed of 100 genes.
    only_two_centers : bool
        If True, split the data in two centers only, by default False.
    design_factors : str or list
        The design factors.
    continuous_factors : list[str] or None
        The continuous design factors. Factors not in this list will be considered
        as categorical.
    heterogeneity_method : str or None
        The method to use to generate heterogeneity in the dataset. If None, no
        heterogeneity is generated. Default is None.
    heterogeneity_method_param : float or None
        The parameter of the heterogeneity method. Default is None.
    **pydeseq2_kwargs
        Additional arguments to pass to the pydeseq2 and fedpydeseq2 classes and
        strategies.

    """
    # Check that the design factors and continuous factors are allowed
    if isinstance(design_factors, str):
        design_factors = [design_factors]
    if continuous_factors is None:
        continuous_factors = []
    assert set(design_factors).issubset(
        ALLOWED_DESIGN_FACTORS_TCGA
    ), f"Design factors should be in {ALLOWED_DESIGN_FACTORS_TCGA}"
    assert set(continuous_factors).issubset(
        ALLOWED_CONTINUOUS_FACTORS_TCGA
    ), f"Continuous factors should be in {ALLOWED_CONTINUOUS_FACTORS_TCGA}"

    processed_data_path = Path(processed_data_path).resolve()
    sampling_random_generator = np.random.default_rng(42)

    if small_genes:
        n_genes = 100
    if small_samples:
        n_samples = 10

    experiment_id = get_experiment_id(
        dataset_name,
        small_samples,
        small_genes,
        only_two_centers,
        design_factors,
        continuous_factors,
        heterogeneity_method,
        heterogeneity_method_param,
        **pydeseq2_kwargs,
    )

    center_data_path = processed_data_path / "centers_data" / "tcga" / experiment_id

    # -- Process the data
    logger.info(f"Processing the data for the TCGA dataset: {experiment_id}")
    counts_data, metadata = preprocess_tcga(
        processed_data_path=processed_data_path,
        dataset_name=dataset_name,
        only_two_centers=only_two_centers,
        design_factors=design_factors,
        heterogeneity_method=heterogeneity_method,
        heterogeneity_method_param=heterogeneity_method_param,
    )

    if small_genes:
        counts_data = counts_data.sample(
            n_genes, axis=1, random_state=sampling_random_generator
        )

    if small_samples:
        new_counts_data_list: list[pd.DataFrame] = []
        new_metadata_list: list[pd.DataFrame] = []

    # -- Split the data
    logger.info(f"Saving the data for each center {center_data_path}")
    for center_id in metadata.center_id.unique():
        counts_dataframe = counts_data.loc[metadata.center_id == center_id]
        metadata_dataframe = metadata.loc[metadata.center_id == center_id]
        if small_samples:
            categorical_factors = list(set(design_factors) - set(continuous_factors))
            n_levels = len(metadata_dataframe[categorical_factors].drop_duplicates())
            n_samples_per_level = max(n_samples // n_levels, 1)

            def _sampling_function(df_, n_samples_, sampling_rng_):
                return df_.sample(
                    min(len(df_), n_samples_),
                    random_state=sampling_rng_,
                    replace=False,
                )

            _partialized_sampling_function = partial(
                _sampling_function,
                n_samples_=n_samples_per_level,
                sampling_rng_=sampling_random_generator,
            )

            metadata_dataframe = (
                metadata_dataframe.groupby(categorical_factors, dropna=False)
                .apply(
                    _partialized_sampling_function,
                    include_groups=False,
                )
                .reset_index(categorical_factors)
            )

            counts_dataframe = counts_dataframe.loc[metadata_dataframe.index]
            new_counts_data_list.append(counts_dataframe)
            new_metadata_list.append(metadata_dataframe)
        path = center_data_path / f"center_{center_id}"

        # delete the folder if it exists
        # This avoids having old files in the folder when registering the data
        # in substra
        if path.exists():
            for file in path.iterdir():
                if file.is_dir():
                    if file.name == ".ipynb_checkpoints":
                        logger.info(f"Removing {file}")
                        shutil.rmtree(file)
                    else:
                        raise ValueError(
                            f"Unexpected directory in the center folder: {file}"
                        )
                else:
                    file.unlink()
            path.rmdir()
        path.mkdir(parents=True)

        counts_dataframe.to_csv(path / "counts_data.csv")
        metadata_dataframe.to_csv(path / "metadata.csv")

    if small_samples:
        counts_data = pd.concat(new_counts_data_list)
        metadata = pd.concat(new_metadata_list)
    # Now save to the pooled data folder
    pooled_data_path = processed_data_path / "pooled_data" / "tcga" / experiment_id
    logger.info(f"Saving the pooled data at {pooled_data_path}")
    pooled_data_path.mkdir(parents=True, exist_ok=True)
    counts_data.to_csv(pooled_data_path / "counts_data.csv")
    metadata.to_csv(pooled_data_path / "metadata.csv")


def preprocess_tcga(
    processed_data_path: Path,
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
    only_two_centers=False,
    design_factors: str | list[str] = "stage",
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Preprocess the TCGA dataset.

    If the `stage` design factor is used, the function will binarize the `stage` design
    into two categories: `Advanced` and `Non-advanced`.
    `Advanced` corresponds to stage `IV`,
    and `Non-advanced` corresponds to stages `I`, `II` and `III`.
    For the TCGA-PRAD cohort, we do not have the stage information, but we infer the
    stage from the `T`, `N` and `M` columns. If the `N` or `M` columns are > 0,
    the stage is IV according to:
    https://www.cancer.org/cancer/types/prostate-cancer\
    /detection-diagnosis-staging/staging.html
    and hence the `Advanced` stage. Otherwise, it is `Non-advanced`.


    Parameters
    ----------
    processed_data_path : Path
        Where to find the folder tcga with raw data.
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
        The dataset to preprocess, by default "TCGA-LUAD".
    only_two_centers : bool
        If True, split the data in two centers only, by default False.
    design_factors : str or list
        The design factors.
    heterogeneity_method : str or None
        The method to use to generate heterogeneity in the dataset. If None, no
        heterogeneity is generated. Default is None.
    heterogeneity_method_param : float or None
        The parameter of the heterogeneity method. Default is None.

    Returns
    -------
    counts_data: pd.DataFrame
        Processed pooled counts dataset.
    metadata: pd.DataFrame
        Processed pooled metadata dataset.

    """
    _, cohort = dataset_name.split("-")[:2]
    if cohort in ["NSCLC", "CRC"]:
        if cohort == "NSCLC":
            dataset_1, dataset_2 = cast(TCGADatasetNames, "TCGA-LUAD"), cast(
                TCGADatasetNames, "TCGA-LUSC"
            )
        elif cohort == "CRC":
            dataset_1, dataset_2 = cast(TCGADatasetNames, "TCGA-COAD"), cast(
                TCGADatasetNames, "TCGA-READ"
            )
        counts_data_1, metadata_1 = preprocess_tcga(
            processed_data_path=processed_data_path,
            dataset_name=dataset_1,
            only_two_centers=only_two_centers,
            design_factors=design_factors,
        )
        metadata_1["center_id"] = 0
        counts_data_2, metadata_2 = preprocess_tcga(
            processed_data_path=processed_data_path,
            dataset_name=dataset_2,
            only_two_centers=only_two_centers,
            design_factors=design_factors,
        )
        metadata_2["center_id"] = 1
        counts_data = pd.concat([counts_data_1, counts_data_2])
        metadata = pd.concat([metadata_1, metadata_2])

        metadata = mix_centers(
            metadata, heterogeneity_method, heterogeneity_method_param
        )
        return counts_data, metadata

    data_path = processed_data_path / "tcga" / cohort
    path_to_counts = data_path / "counts.parquet"
    path_to_metadata = data_path / "clinical_data.csv"
    # -- Load the data
    counts_data = pd.read_parquet(path_to_counts)

    metadata = pd.read_csv(path_to_metadata, index_col=0)

    # -- Process the metadata

    # If the cohort is PRAD, we need to add the stage information
    if dataset_name == "TCGA-PRAD":

        def binarize_stage_prad(metadata_row):
            if (
                pd.isna(metadata_row["M"])
                and pd.isna(metadata_row["N"])
                and pd.isna(metadata_row["T"])
            ):
                return pd.NA
            if (not pd.isna(metadata_row["N"])) and (metadata_row["N"] > 0):
                return "Advanced"
            elif (not pd.isna(metadata_row["M"])) and (metadata_row["M"] > 0):
                return "Advanced"
            else:
                return "Non-advanced"

        metadata["stage"] = metadata.apply(binarize_stage_prad, axis=1)

    else:
        # Binarize the stage
        def binarize_stage(stage):
            if pd.isna(stage):
                return pd.NA
            elif stage == 4:
                return "Advanced"
            else:
                return "Non-advanced"

        metadata["stage"] = metadata["stage"].apply(binarize_stage)

    # -- Process the data
    cols_to_keep = (
        design_factors.copy() if isinstance(design_factors, list) else [design_factors]
    )
    cols_to_keep.append("center_id")
    metadata = metadata[cols_to_keep]

    # remove samples with NaN design factor
    metadata.dropna(subset=design_factors, inplace=True)

    if only_two_centers:
        metadata = _merge_centers(metadata)

    counts_data = counts_data.loc[metadata.index]

    return counts_data, metadata


def _merge_centers(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Merge the centers into two centers.

    Parameters
    ----------
    metadata : pd.DataFrame
        Metadata

    Returns
    -------
    metadata : pd.DataFrame
        Metadata
    """
    list_center_ids = np.sort(metadata.center_id.unique())
    assert len(list_center_ids) > 2, "The dataset has less than 2 centers"
    dict_mapping = {list_center_ids[k]: 0 for k in range(len(list_center_ids) // 2)}
    dict_mapping.update(
        {
            list_center_ids[k]: 1
            for k in range(len(list_center_ids) // 2, len(list_center_ids))
        }
    )
    logger.info(f"Merging centers using the mapping: {dict_mapping}")
    metadata["center_id"] = metadata["center_id"].map(dict_mapping)

    return metadata
