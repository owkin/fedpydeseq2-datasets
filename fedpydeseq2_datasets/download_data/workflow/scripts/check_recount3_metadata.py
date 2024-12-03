# %%
import sys

import pandas as pd


def check_recount3_metadata(input_file: str, output_file: str):
    """
    Convert the TCGA clinical data to a CSV file and check its content.

    Parameters
    ----------
    input_file : str
        Path to the input file, containing metadata.
    output_file : str
        Path to the output file, containing the metadata in CSV format.

    """
    df = pd.read_csv(input_file, sep="\t")
    # Check the columns
    assert "external_id" in df.columns
    assert "tcga.tcga_barcode" in df.columns

    # Save the dataframe to a tsv file
    df.to_csv(output_file, index=False, sep="\t", compression="gzip")

    return


with open(snakemake.log[0], "w") as f:  # type: ignore # noqa: F821
    sys.stderr = sys.stdout = f
    try:
        input_file = snakemake.input[0]  # type: ignore  # noqa: F821
        output_file = snakemake.output[0]  # type: ignore # noqa: F821
        check_recount3_metadata(input_file, output_file)
    except Exception as e:  # noqa: BLE001
        print(e)
        sys.exit(1)

# %%
