import pytest

from fedpydeseq2_datasets.create_reference_dds import setup_tcga_ground_truth_dds
from fedpydeseq2_datasets.process_and_split_data import setup_tcga_dataset


@pytest.mark.parametrize(
    "dataset_name",
    [
        "TCGA-LUAD",
    ],
)
@pytest.mark.usefixtures("raw_data_path", "tmp_processed_data_path")
def test_tcga_preprocessing_and_pooled_pydeseq2_all_indications_small(
    raw_data_path, tmp_processed_data_path, dataset_name
):
    """Build all necessary quantities for all tests.

    Builds ground truths.
    """

    setup_tcga_dataset(
        raw_data_path,
        tmp_processed_data_path,
        dataset_name=dataset_name,
        small_samples=True,
        small_genes=True,
        only_two_centers=False,
        design_factors="stage",
        continuous_factors=None,
        refit_cooks=True,
        force=True,
    )
    print("Setting up TCGA ground truth DESeq2 datasets...")
    setup_tcga_ground_truth_dds(
        tmp_processed_data_path,
        dataset_name=dataset_name,
        small_samples=True,
        small_genes=True,
        only_two_centers=False,
        design_factors="stage",
        continuous_factors=None,
        reference_dds_ref_level=("stage", "Advanced"),
        refit_cooks=True,
        force=True,
    )


@pytest.mark.parametrize(
    "dataset_name",
    [
        "TCGA-LUAD",
    ],
)
@pytest.mark.usefixtures("raw_data_path", "tmp_processed_data_path")
def test_tcga_preprocessing_and_per_center_pydeseq2_all_indications_small(
    raw_data_path, tmp_processed_data_path, dataset_name
):
    """Build all necessary quantities for all tests.

    Builds ground truths.
    """

    setup_tcga_dataset(
        raw_data_path,
        tmp_processed_data_path,
        dataset_name=dataset_name,
        small_samples=True,
        small_genes=True,
        only_two_centers=False,
        design_factors="stage",
        continuous_factors=None,
        refit_cooks=True,
        force=True,
    )
    print("Setting up TCGA ground truth DESeq2 datasets...")
    setup_tcga_ground_truth_dds(
        tmp_processed_data_path,
        dataset_name=dataset_name,
        small_samples=True,
        small_genes=True,
        only_two_centers=False,
        design_factors="stage",
        continuous_factors=None,
        reference_dds_ref_level=("stage", "Advanced"),
        refit_cooks=True,
        force=True,
        pooled=False,
    )
