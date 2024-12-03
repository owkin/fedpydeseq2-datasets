# %%
import sys

import pandas as pd


def csvize_and_check_clinical(input_file: str, output_file: str):
    """
    Convert the TCGA clinical data to a CSV file and check its content.

    Parameters
    ----------
    input_file : str
        Path to the input file, containing metadata.
    output_file : str
        Path to the output file, containing the metadata in CSV format.

    """
    df = pd.read_excel(input_file, index_col=0)

    # Check it contains the right columns
    assert "bcr_patient_barcode" in df.columns
    assert "gender" in df.columns
    assert "ajcc_pathologic_tumor_stage" in df.columns

    df.to_csv(output_file, index=False)

    return


with open(snakemake.log[0], "w") as f:  # type: ignore # noqa: F821
    sys.stderr = sys.stdout = f
    try:
        input_file = snakemake.input[0]  # type: ignore  # noqa: F821
        output_file = snakemake.output[0]  # type: ignore # noqa: F821
        csvize_and_check_clinical(input_file, output_file)
    except Exception as e:  # noqa: BLE001
        print(e)
        sys.exit(1)

# %%
