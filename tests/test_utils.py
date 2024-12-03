import numpy as np
import pandas as pd
import pytest

from fedpydeseq2_datasets.utils import mix_centers


@pytest.mark.parametrize(
    "heterogeneity_method, heterogeneity_method_param",
    [
        ("binomial", 0.0),
        ("binomial", 0.2),
        ("binomial", 0.5),
        ("binomial", 0.8),
        ("binomial", 1),
    ],
)
def test_mix_centers(heterogeneity_method, heterogeneity_method_param):
    # Create a sample metadata DataFrame
    np.random.seed(0)
    metadata = pd.DataFrame(
        {
            "SampleID": [f"sample_{i}" for i in range(1000)],
            "center_id": np.random.choice([0, 1], size=1000),
        }
    )
    old_metadata = metadata.copy()
    # Call the mix_centers function
    mix_centers(metadata, heterogeneity_method, heterogeneity_method_param)

    cross_table = pd.crosstab(metadata["center_id"], old_metadata["center_id"])

    # Expected proportions
    expected_cross_table = np.zeros((2, 2))
    center_0_count = sum(old_metadata["center_id"] == 0)
    center_1_count = sum(old_metadata["center_id"] == 1)

    # Calculate the expected values based on the binomial distribution
    expected_cross_table[0, 0] = center_0_count * (1 - heterogeneity_method_param / 2.0)
    expected_cross_table[0, 1] = center_1_count * (heterogeneity_method_param / 2.0)
    expected_cross_table[1, 0] = center_0_count * (heterogeneity_method_param / 2.0)
    expected_cross_table[1, 1] = center_1_count * (1 - heterogeneity_method_param / 2.0)

    # Allow some tolerance due to randomness
    tolerance = 0.1 * (center_0_count + center_1_count)

    # Check that the observed values are close to the expected values
    assert np.allclose(cross_table.values, expected_cross_table, atol=tolerance)
