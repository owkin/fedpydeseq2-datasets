"""Module to preprocess the TCGA data for fedpydeseq2."""

from fedpydeseq2_datasets.process_and_split_data import setup_tcga_dataset
from fedpydeseq2_datasets.create_reference_dds import setup_tcga_ground_truth_dds
from fedpydeseq2_datasets.utils import get_experiment_id
from fedpydeseq2_datasets.utils import get_ground_truth_dds_name
