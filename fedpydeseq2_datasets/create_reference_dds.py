import pickle as pkl
from inspect import signature
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger
from pydeseq2.dds import DeseqDataSet

from fedpydeseq2_datasets.constants import TCGADatasetNames
from fedpydeseq2_datasets.utils import get_experiment_id
from fedpydeseq2_datasets.utils import get_ground_truth_dds_name

ALLOWED_DESIGN_FACTORS_TCGA = {"stage", "gender", "CPE"}
ALLOWED_CONTINUOUS_FACTORS_TCGA = {"CPE"}


def setup_tcga_ground_truth_dds(
    processed_data_path: str | Path,
    dataset_name: TCGADatasetNames = "TCGA-LUAD",
    small_samples: bool = False,
    small_genes: bool = False,
    only_two_centers: bool = False,
    design_factors: str | list[str] = "stage",
    continuous_factors: list[str] | None = None,
    reference_dds_ref_level: tuple[str, ...] | None = ("stage", "Advanced"),
    force: bool = False,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
    pooled: bool = True,
    default_refit_cooks: bool = False,
    **pydeseq2_kwargs: Any,
):
    """Set the ground truth DeseqDataSet for the TCGA dataset.

    This function is given the path to the processed data.
    Then it preprocesses the data and initializes the DeseqDataSet.
    WARNING: by default, the cooks outliers are NOT refitted.
    Afterward, it performs the Deseq2 pipeline.
    Finally, it saves the pooled_dds and all the local_dds on the disk.
    The local dds is simply the restriction of the global dds to the center id, for
    each center.


    The file structure for the processed data is the following (result
    from the `setup_tcga_dataset` function):

    ```
    <processed_data_path>
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

    ```

    The file structure output for the processed data is the following:

    ```
    <processed_data_path>
    ├── centers_data
    │   └── tcga
    │       └── <experiment_id>
    │           └── center_0
    │               ├── counts_data.csv
    │               ├── metadata.csv
    │               └── ground_truth_dds.pkl # if pooled
    │               └── ground_truth_dds-center.pkl # if not pooled

    ├── pooled_data
    │   └── tcga
    │       └── <experiment_id>
    │           ├── counts_data.csv
    │           ├── metadata.csv
    │           └── ground_truth_dds.pkl # if pooled
    └──

    ```

    Parameters
    ----------
    processed_data_path : str or Path
        Path to processed data.
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
        The continuous design factors.
    reference_dds_ref_level : tuple or None
        The reference level for the design factor.
    force : bool
        If True, force the setup of the dataset even if the metadata already exists.
        Default is False.
    heterogeneity_method : str or None
        The method to use to generate heterogeneity in the dataset. If None, no
        heterogeneity is generated. Default is None.
    heterogeneity_method_param : float or None
        The parameter of the heterogeneity method. Default is None.
    pooled : bool
        If True, We compute the pooled dds, and then restrict it to the
        center to get the local dds. If False, we only compute the local dds
        from the center data only. Default is True.
    default_refit_cooks : bool
        If True, refit the cooks outliers. Default is False.
    **pydeseq2_kwargs
        Additional arguments to pass to the pydeseq2 and fedpydeseq2 classes and
        strategies.

    """
    processed_data_path = Path(processed_data_path).resolve()
    logger.info(f"Design factors in setup_tcga_ground_truth_dds: {design_factors}")

    refit_cooks = pydeseq2_kwargs.get("refit_cooks", default_refit_cooks)

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

    pooled_data_path = processed_data_path / "pooled_data" / "tcga" / experiment_id

    ground_truth_dds_name = get_ground_truth_dds_name(
        reference_dds_ref_level, refit_cooks=refit_cooks, pooled=pooled
    )

    pooled_dds_file_path = pooled_data_path / f"{ground_truth_dds_name}.pkl"

    # Check that all centers dds files were generated
    metadata = pd.read_csv(
        processed_data_path / "pooled_data" / "tcga" / experiment_id / "metadata.csv",
        index_col=0,
    )

    if pooled and (not pooled_dds_file_path.exists() or force):
        # -- Process the data and initialize the DeseqDataSet
        dds_file_name = get_ground_truth_dds_name(
            reference_dds_ref_level, refit_cooks=refit_cooks, pooled=True
        )

        pooled_dds_file_path = pooled_data_path / f"{dds_file_name}.pkl"

        # We pass the default refit_cooks to the setup_ground_truth_dds function
        # It is overwritten by the refit_cooks parameter if it is passed
        # as a pydesq2_kwargs
        setup_ground_truth_dds_kwargs = {"refit_cooks": refit_cooks, **pydeseq2_kwargs}

        _setup_ground_truth_dds(
            pooled_data_path,
            pooled_dds_file_path,
            design_factors,
            continuous_factors,
            reference_dds_ref_level,
            **setup_ground_truth_dds_kwargs,
        )

    for center_id in metadata.center_id.unique():
        center_data_path = (
            processed_data_path
            / "centers_data"
            / "tcga"
            / experiment_id
            / f"center_{center_id}"
        )
        center_dds_file_path = center_data_path / f"{ground_truth_dds_name}.pkl"
        if not center_dds_file_path.exists() or force:
            if pooled:
                # In that case we need to reprocess the data
                _setup_local_ground_truth_dds(
                    pooled_dds_file_path,
                    center_data_path.parent,
                    metadata,
                )
                break
            # Else, build the local dds from the center data only
            try:
                _setup_ground_truth_dds(
                    center_data_path,
                    center_dds_file_path,
                    design_factors,
                    continuous_factors,
                    reference_dds_ref_level,
                    **{"refit_cooks": refit_cooks, **pydeseq2_kwargs},
                )
            except ValueError as e:
                logger.warning(
                    f"Error while setting up the local dds for center {center_id}: "
                    f"{e}, will set None for this center"
                )
                with open(center_dds_file_path, "wb") as file:
                    pkl.dump(None, file)


def _setup_local_ground_truth_dds(
    pooled_dds_file_path: Path,
    center_data_path: Path,
    metadata: pd.DataFrame,
):
    """Set the local ground truth DeseqDataSet for the TCGA dataset.

    This function is given the path to the pooled_dds and the path to the center data.
    Then it loads the pooled_dds and initializes the local DeseqDataSet for each center.
    Finally, it saves the local_dds on the disk.

    For the file structure description, see the `setup_tcga_ground_truth_dds` function.

    Parameters
    ----------
    pooled_dds_file_path : Path
        Path to the pooled_dds.
    center_data_path : Path
        Path to the center data.
    metadata : pd.DataFrame
        Metadata of the pooled data.
    """
    with open(pooled_dds_file_path, "rb") as file:
        pooled_dds = pkl.load(file)

    for k in metadata.center_id.unique():
        local_reference_dds = pooled_dds[pooled_dds.obs.center_id == k].copy()
        path = center_data_path / f"center_{k}"
        path.mkdir(parents=True, exist_ok=True)
        with open(path / pooled_dds_file_path.name, "wb") as file:
            pkl.dump(local_reference_dds, file)


def _setup_ground_truth_dds(
    data_path: Path,
    output_file_path: Path,
    design_factors: str | list[str],
    continuous_factors: list[str] | None,
    reference_dds_ref_level: tuple[str, ...] | None,
    **pydeseq2_kwargs,
):
    """Process the data and initialize the DeseqDataSet.

    This function is given the path to the pooled data.
    Then it loads the data and initializes the DeseqDataSet.
    Afterward, it performs the Deseq2 pipeline.
    Finally, it saves the pooled_dds on the disk.

    For the file structure description, see the `setup_tcga_ground_truth_dds` function.

    Parameters
    ----------
    data_path : Path
        Path to the data necessary to initialize the DeseqDataSet.
        We require that the data path be a directory containing the following files:
        - counts_data.csv
        - metadata.csv

    output_file_path : Path
        Path to save the DeseqDataSet.
        It is expected to be a pkl file.

    design_factors : str or list
        The design factors.

    continuous_factors : list[str] or None
        The list of continuous factors. Factors not in this list will be considered
        as categorical.

    reference_dds_ref_level : tuple or None
        The reference level for the design factor.

    **pydeseq2_kwargs
        Additional arguments to pass to the pydeseq2 and fedpydeseq2 classes and
        strategies.

    """
    counts_data = pd.read_csv(
        data_path / "counts_data.csv",
        index_col=0,
    )
    metadata = pd.read_csv(
        data_path / "metadata.csv",
        index_col=0,
    )
    if "center_id" in counts_data.columns:
        counts_data.drop(columns="center_id", inplace=True)

    dds_kwargs = {
        parameter_name: parameter_value
        for parameter_name, parameter_value in pydeseq2_kwargs.items()
        if parameter_name in signature(DeseqDataSet).parameters
    }

    ref_level = (
        list(reference_dds_ref_level) if reference_dds_ref_level is not None else None
    )
    pooled_dds = DeseqDataSet(
        counts=counts_data,
        metadata=metadata,
        design_factors=design_factors,
        ref_level=ref_level,
        continuous_factors=continuous_factors,
        **dds_kwargs,
    )

    # log the dds parameters
    logger.info("Creating a DeseqDataSet with the following parameters:")
    logger.info(f"Counts path { data_path / 'counts_data.csv'}")
    logger.info(f"Metadata path { data_path / 'metadata.csv'}")
    logger.info(f"Design factors {design_factors}")
    logger.info(f"Continuous factors {continuous_factors}")
    logger.info(f"Reference level {ref_level}")
    logger.info(f"PyDESeq2 kwargs {pydeseq2_kwargs}")

    # -- Perform the Deseq2 pipeline.

    # Compute DESeq2 normalization factors using the Median-of-ratios method
    pooled_dds.fit_size_factors()
    # Fit an independent negative binomial model per gene
    pooled_dds.fit_genewise_dispersions()
    # Fit a parameterized trend curve for dispersions, of the form
    # f(\mu) = \alpha_1/\mu + a_0
    pooled_dds.fit_dispersion_trend()
    disp_function_type = pooled_dds.uns["disp_function_type"]
    logger.info(f"Dispersion function type: {disp_function_type}")
    # Compute prior dispersion variance
    pooled_dds.fit_dispersion_prior()
    # Refit genewise dispersions a posteriori (shrinks estimates towards trend curve)
    pooled_dds.fit_MAP_dispersions()
    # Fit log-fold changes (in natural log scale)
    pooled_dds.fit_LFC()

    pooled_dds.calculate_cooks()

    if pooled_dds.refit_cooks:
        logger.info("DDSEq2: Refitting cooks outliers")
        pooled_dds.refit()
        # Here, we must change the replaced genes to a boolean array
        # As AnnData does not support the series version of replaced
        # defined in deseq2 (issues when copying)
        if pooled_dds.obsm["replaceable"].sum() == 0:
            pooled_dds.varm["replaced"] = np.zeros((pooled_dds.n_vars,), dtype=bool)

    # -- Save the pooled_dds

    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Saving the pooled dds at {output_file_path}")
    with open(output_file_path, "wb") as file:
        pkl.dump(pooled_dds, file)
