import copy
import pickle
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from fedpydeseq2_datasets.constants import TCGADatasetNames

IDENTIFYING_PARAMETERS = [
    "contrast",
    "refit_cooks",
    "alt_hypothesis",
    "lfc_null",
]


def tnm_to_series(tnm_string: str | None) -> pd.Series:
    """Convert a TNM string into a pandas Series with T, N, and M categories.

    If there are multiple categories for the same letter, the maximum number is used.
    A typical TNM string looks like 'T1bN2M0'.

    Parameters
    ----------
    tnm_string : str or None
        The TNM string to convert.

    Returns
    -------
    pd.Series
        A pandas Series with the T, N, and M categories.
    """
    if pd.isna(tnm_string):
        return pd.Series({"T": pd.NA, "N": pd.NA, "M": pd.NA})

    assert isinstance(tnm_string, str), "The TNM string must be a string."
    # Split the string into separate TNM categories
    tnm_categories = re.findall(r"[TNM]\d+[a-z]*", tnm_string)

    # Initialize the maximum numbers for T, N, and M as None
    max_t = max_n = max_m = None

    # Iterate over each category
    for category in tnm_categories:
        # Extract the category letter and number
        letter = category[0]
        number = int(re.search(r"\d+", category).group())  # type: ignore

        # Update the maximum number for the corresponding category
        if letter == "T":
            max_t = max(max_t, number) if max_t is not None else number  # type: ignore
        elif letter == "N":
            max_n = max(max_n, number) if max_n is not None else number  # type: ignore
        elif letter == "M":
            max_m = max(max_m, number) if max_m is not None else number  # type: ignore

    # Convert the maximum numbers into a pandas Series
    series = pd.Series(
        {
            "T": max_t if max_t is not None else pd.NA,
            "N": max_n if max_n is not None else pd.NA,
            "M": max_m if max_m is not None else pd.NA,
        }
    )

    return series


def get_experiment_id(
    dataset_name: TCGADatasetNames | list[TCGADatasetNames],
    small_samples: bool,
    small_genes: bool,
    only_two_centers: bool,
    design_factors: str | list[str] = "stage",
    continuous_factors: list[str] | None = None,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
    **pydeseq2_kwargs: Any,
):
    """
    Generate the experiment id.

    Parameters
    ----------
    dataset_name : TCGADatasetNames or list[TCGADatasetNames]
        Dataset name
    small_samples : bool
        If True, only preprocess a small subset of the data, by default False.
        This small subset is composed of 10 samples per center.
    small_genes : bool
        If True, only preprocess a small subset of the data features (genes)
    only_two_centers : bool
        If True, split the data in two centers only.
    design_factors : str or list
        The design factors.
    continuous_factors : list or None
        The continuous design factors. Factors not in this list will be considered
        as categorical.
    heterogeneity_method : str or None
        The method used to generate heterogeneity in the data.
    heterogeneity_method_param : float or None
        The parameter for the heterogeneity method.
    **pydeseq2_kwargs : Any
        Additional arguments to pass to pydeseq2.

    Returns
    -------
    experiment_id : str
        The true dataset name
    """
    if isinstance(dataset_name, list):
        # Create a string with all the dataset names
        # Get the tcga cohorts
        all_cohorts = "-".join(
            sorted([dataset_name.split("-")[1] for dataset_name in dataset_name])
        )
        full_dataset_name = f"TCGA-{all_cohorts}"
    else:
        full_dataset_name = dataset_name
    if isinstance(design_factors, str):
        design_factors = [design_factors]
    if continuous_factors is None:
        continuous_factors = []
    design_factor_str = "_".join(design_factors)
    continuous_factor_str = "_".join(continuous_factors)
    small_genes_postfix = "-small-genes" if small_genes else ""
    small_samples_postfix = "-small-samples" if small_samples else ""
    two_centers_postfix = "-two-centers" if only_two_centers else ""

    pydeseq2_kwargs = copy.deepcopy(pydeseq2_kwargs)
    if len(design_factors) > 1:
        # Add the contrast to the pydeseq2_kwargs
        pydeseq2_kwargs["contrast"] = pydeseq2_kwargs.get("contrast", None)

    # Sort the kwargs alphabetically
    parameter_names = sorted(pydeseq2_kwargs.keys())
    parameter_str = ""
    for parameter_name in parameter_names:
        if parameter_name not in IDENTIFYING_PARAMETERS:
            continue
        parameter_value = pydeseq2_kwargs[parameter_name]
        if parameter_value is None:
            parameter_str += f"-default_{parameter_name}"
        elif isinstance(parameter_value, list):
            parameter_str += f"-{parameter_name}-{'_'.join(parameter_value)}"
        else:
            parameter_str += f"-{parameter_name}-{parameter_value}"

    heterogeneity_postfix = (
        f"heterogeneity-{heterogeneity_method}"
        if heterogeneity_method is not None
        else ""
    )
    heterogeneity_method_param_postfix = (
        f"-{heterogeneity_method_param}"
        if heterogeneity_method_param is not None
        else ""
    )

    experiment_id = (
        f"{full_dataset_name}-{design_factor_str}-{continuous_factor_str}"
        f"{small_genes_postfix}{small_samples_postfix}{two_centers_postfix}"
        f"{heterogeneity_postfix}{heterogeneity_method_param_postfix}{parameter_str}"
    )
    experiment_id = experiment_id.strip("-_")  # Remove trailing '-' and '_'

    return experiment_id


def get_ground_truth_dds_name(
    reference_dds_ref_level: tuple[str, ...] | None = ("stage", "Advanced"),
    refit_cooks: bool = False,
    pooled: bool = True,
) -> str:
    """
    Generate the ground truth dds name.

    Parameters
    ----------
    reference_dds_ref_level : tuple or None
        The reference level for the design factor.
    refit_cooks : bool
        If True, refit the genes with cooks outliers.
        TODO this is now obsolete, we should remove it at the end
        TODO (but for compatibility reasons we keep it for now).
    pooled : bool
        If True, we compute the pooled dds, and then restrict it to the
        center to get the local dds. If False, we only compute the local dds
        from the center data only.
        If not pooled, we add the suffix '-center' to the name.

    Returns
    -------
    ground_truth_dds_name : str
        The ground truth dds name
    """
    ref_level_str = (
        "_".join(reference_dds_ref_level)
        if reference_dds_ref_level is not None
        else "None"
    )
    ground_truth_dds_name = f"ground_truth_dds-{ref_level_str}"
    if refit_cooks:
        ground_truth_dds_name += "-refit_cooks"
    if not pooled:
        ground_truth_dds_name += "-center"
    return ground_truth_dds_name


def get_n_centers(
    processed_data_path: str | Path,
    dataset_name: TCGADatasetNames,
    small_samples: bool = False,
    small_genes: bool = False,
    only_two_centers: bool = False,
    design_factors: str | list[str] = "stage",
    continuous_factors: list[str] | None = None,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
    **pydeseq2_kwargs: Any,
) -> int:
    """
    Get the number of centers in the dataset.

    To do so, we open the file containing the metadata and count the number
    of unique center ids and count them.

    Parameters
    ----------
    processed_data_path : str or Path
        The path to the processed data.

    dataset_name : TCGADatasetNames
        The name of the dataset to use.

    small_samples : bool
        Whether to use a small number of samples. Default is False.

    small_genes : bool
        Whether to use a small number of genes. Default is False.

    only_two_centers : bool
        Whether to use only two centers. Default is False.

    design_factors : str or list[str]
        The design factors to use. Default is "stage".

    continuous_factors : list[str] or None
        The continuous factors to use. Default is None.

    heterogeneity_method : str or None
        The method to used to define the heterogeneity
        of the center's attribution.

    heterogeneity_method_param : float or None
        The parameter of the heterogeneity method.

    **pydeseq2_kwargs : Any
        Additional arguments to pass to pydeseq2 and fedpydeseq2.
        For example, the contrast.

    Returns
    -------
    n_centers : int
        The number of centers in the dataset

    """
    # Get the number of centers
    processed_data_path = Path(processed_data_path)
    experiment_id = get_experiment_id(
        dataset_name=dataset_name,
        small_samples=small_samples,
        small_genes=small_genes,
        only_two_centers=only_two_centers,
        design_factors=design_factors,
        continuous_factors=continuous_factors,
        heterogeneity_method=heterogeneity_method,
        heterogeneity_method_param=heterogeneity_method_param,
        **pydeseq2_kwargs,
    )
    n_centers = len(
        np.unique(
            pd.read_csv(
                processed_data_path
                / "pooled_data"
                / "tcga"
                / experiment_id
                / "metadata.csv"
            )["center_id"]
        )
    )
    return n_centers


def get_n_centers_from_subfolders(centers_path: str | Path) -> int:
    """
    Get the number of centers from a folder.

    Parameters
    ----------
    centers_path : str or Path
        The path to the folder which must contain
        center_{i} subfolders

    Returns
    -------
    int
        The number of centers
    """
    centers_path = Path(centers_path)
    n_centers = len([x for x in centers_path.iterdir() if x.name.startswith("center_")])
    # check that the names are correct
    assert all((centers_path / f"center_{i}").exists() for i in range(n_centers))
    return n_centers


def get_valid_centers_from_subfolders_file(
    centers_path: str | Path, filename: str, pkl: bool = False
) -> tuple[int, list[int]]:
    """
    Get the number of centers from a folder.

    Parameters
    ----------
    centers_path : str or Path
        The path to the folder which must contain
        center_{i} subfolders

    filename : str
        The name of the file to check for in the center_{i} subfolders.
        If the file is not found in a center_{i} subfolder, the center is
        not considered valid (in practice, this will mean that a DeSEQ2 analysis
        will not have been run on the center for lack of samples of all
        design factors). If pkl is True, the file is assumed to be a pickle file,
        and the file is loaded with pickle.load() to check if it is None.
        If it is None, the center is not considered valid.

    pkl : bool
        If True, the file is assumed to be a pickle file, and the file is loaded
        with pickle.load() to check if it is None. If it is None, the center is
        not considered valid.

    Returns
    -------
    int
        The number of centers
    """
    centers_path = Path(centers_path)
    n_centers = get_n_centers_from_subfolders(centers_path)
    valid_centers = []
    for i in range(n_centers):
        if (centers_path / f"center_{i}" / filename).exists():
            if pkl:
                with open(centers_path / f"center_{i}" / filename, "rb") as f:
                    data = pickle.load(f)
                if data is not None:
                    valid_centers.append(i)
            else:
                valid_centers.append(i)
    return n_centers, valid_centers


def mix_centers(
    metadata: pd.DataFrame,
    heterogeneity_method: str | None = None,
    heterogeneity_method_param: float | None = None,
):
    """
    Mix the centers in the metadata.

    It uses the heterogeneity method and the heterogeneity
    method parameter.

    Parameters
    ----------
    metadata : pd.DataFrame
        The metadata to mix the centers
    heterogeneity_method : Optional[str]
        The method used to generate heterogeneity in the data. The only option
        supported is 'binomial'.
    heterogeneity_method_param : Optional[float]
        The parameter for the heterogeneity method. Should be between 0 and 1. If
        the value is 0, the centers are not mixed. If the value is 1, the centers
        are completely mixed.

    Returns
    -------
    metadata : pd.DataFrame
        The metadata with the centers mixed

    """
    if heterogeneity_method is None:
        logger.info("No heterogeneity method specified. Keeping the centers as is.")
        return metadata
    assert (
        heterogeneity_method_param is not None
    ), "The heterogeneity method parameter must be specified."
    assert (
        0 <= heterogeneity_method_param <= 1
    ), "The heterogeneity method parameter must be between 0 and 1."
    np.random.seed(42)
    assert set(metadata.center_id.unique()) == {
        0,
        1,
    }, "The metadata must contain only two centers."
    if heterogeneity_method == "binomial":
        assert heterogeneity_method_param is not None
        mask_center_0 = metadata.center_id == 0
        mask_center_1 = metadata.center_id == 1
        random_binomial_0 = (
            np.random.binomial(1, heterogeneity_method_param / 2.0, mask_center_0.sum())
        ).astype(int)
        random_binomial_1 = (
            np.random.binomial(
                1, 1.0 - heterogeneity_method_param / 2.0, mask_center_1.sum()
            )
        ).astype(int)
        metadata["new_center_id"] = np.zeros(len(metadata), dtype=int)
        metadata.loc[mask_center_0, "new_center_id"] = random_binomial_0
        metadata.loc[mask_center_1, "new_center_id"] = random_binomial_1
        cross_table = pd.crosstab(metadata["center_id"], metadata["new_center_id"])
        logger.info("Mixing centers results")
        logger.info(cross_table)
        metadata.drop(columns=["center_id"], inplace=True)
        metadata.rename(columns={"new_center_id": "center_id"}, inplace=True)

    else:
        raise ValueError(
            f"Unknown heterogeneity method parameter: {heterogeneity_method_param}"
        )
    return metadata
