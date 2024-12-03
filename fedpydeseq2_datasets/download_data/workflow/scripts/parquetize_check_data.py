# %%
import sys

import pandas as pd


def parquetize_and_check(input_file: str, gene_names: str, output_file: str):
    """
    Convert a tabular file to a parquet file and check its content.

    Parameters
    ----------
    input_file : str
        Path to the input file, containing the counts.
    gene_names : str
        Path to the file containing the gene names, with two columns:
        'gene_name' corresponding to HGNC and 'gene_id' corresponding to ENSEMBL.
    output_file : str
        Path to the output file.

    """
    df = pd.read_table(input_file, index_col=0)
    # Use multiindex to use both HGNC and Ensembl gene names
    df.index = pd.MultiIndex.from_frame(pd.read_table(gene_names))
    df = df.astype("int32")
    # Check that the number of genes is roughly between 40 000 and 70 000
    assert 40000 < df.shape[0] < 70000
    # Check that the number of samples is greater than 10 and less than 1000
    assert 10 < df.shape[1] < 1300
    df.to_parquet(output_file)
    return


with open(snakemake.log[0], "w") as f:  # type: ignore # noqa: F821
    sys.stderr = sys.stdout = f
    try:
        input_file = snakemake.input[0]  # type: ignore  # noqa: F821
        gene_names = snakemake.input[1]  # type: ignore # noqa: F821
        output_file = snakemake.output[0]  # type: ignore # noqa: F821
        parquetize_and_check(input_file, gene_names, output_file)
    except Exception as e:  # noqa: BLE001
        print(e)
        sys.exit(1)

# %%
