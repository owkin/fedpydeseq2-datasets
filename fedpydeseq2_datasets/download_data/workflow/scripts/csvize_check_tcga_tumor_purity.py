# %%
import sys

import pandas as pd


def csvize_and_check_tumor_purity(input_file: str, output_file: str):
    """
    Convert tumor purity metadata to a CSV file and check its content.

    Parameters
    ----------
    input_file : str
        Path to the input file, containing metadata.
    output_file : str
        Path to the output file, containing the metadata in CSV format.

    """
    df = pd.read_excel(input_file)
    # Remove first two rows
    df = df.iloc[2:]
    # Remove last column
    df = df.iloc[:, :-1]
    # Set first line as header
    df.columns = df.iloc[0]
    # Remove first line
    df = df.iloc[1:]
    # Check that the columns are the expected ones
    assert "Sample ID" in df.columns
    assert "CPE" in df.columns

    df.to_csv(output_file, index=False)
    return


with open(snakemake.log[0], "w") as f:  # type: ignore # noqa: F821
    sys.stderr = sys.stdout = f
    try:
        input_file = snakemake.input[0]  # type: ignore  # noqa: F821
        output_file = snakemake.output[0]  # type: ignore # noqa: F821
        csvize_and_check_tumor_purity(input_file, output_file)
    except Exception as e:  # noqa: BLE001
        print(e)
        sys.exit(1)

# %%
