# Data download

Repository for downloading the raw data necessary to run the tests and experiments.

This directory contains a [snakemake](https://snakemake.readthedocs.io/) pipeline for downloading RNA-seq data from [RECOUNT3](https://rna.recount.bio/),
tumor purity metadata for TCGA from the [Systematic pan-cancer analysis of tumour purity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/) paper,
and cleaned clinical metadata for TCGA from the [An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics](https://www.sciencedirect.com/science/article/pii/S0092867418302290#app2)
paper.


## Setup
The only configuration needed is the list of project to be downloaded, specified in [config/config.yaml](config/config.yaml).
The gene nomenclature is "ENSEMBL".
By default the pipeline downloads the data necessary to run the experiments of the paper,
that is the following TCGA cohorts:
```yaml
datasets :
  - LUAD
  - LUSC
  - PAAD
  - COAD
  - BRCA
  - PRAD
  - READ
  - SKCM
```


## Execution
For execution, you will need to clone this repository,
and go to the root of this README, that is `download_data`

Note that this pipeline requires a working installation of `conda`,
and requires access to certain `R` packages through `curl`.

### Directly running the full pipeline  with the `fedpydeseq2-download-data` command
If you want to run the pipeline directly, you can use the script which is available in the distribution: `fedpydeseq2-download-data`



```bash
fedpydeseq2-download-data
```

By default, this script download the data in the `data/raw` directory at the root of the github repo.

To change the location of the raw data download, add the following option:
```bash
fedpydeseq2-download-data --raw_data_output_path <path>
```

If you only want the LUAD dataset, add the `--only_luad` flag.

You can pass the `conda` activation path as an argument as well, for example:

```bash
fedpydeseq2-download-data --raw_data_output_path <path> --conda_activate_path /opt/miniconda/bin/activate
```




If you encounter errors, we recommend you use the step by step protocol below.

### A step by step run of the pipeline


To run the pipeline you will need to create a snakemake conda environment :
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge mamba
```

To launch the pipeline :
```bash
snakemake --cores all --resources download_concurrent=3 --use-conda
```
- `--cores all` specifies to use all available cores for parallelization.
- `--resouces download_concurrent=3` means that only 3 download scripts can be run in parallel per physical host. A download script has to be run for each project in the config. It is limited due to some weird behavior of the recount3 R package, and because the download speed is the bottleneck for this rule.
- `--use-conda` means that each rule will use the corresponding conda env, as specified in [workflow/envs/](workflow/envs/). The envs will be created on the first pipeline run, which can take some time.

### Options
See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for all useful command line arguments.

## Details on the `fedpydeseq2-download-data` script

### Overview

This script, `download_data.py`, is part of the `fedpydeseq2_datasets` package and is designed to download data using Snakemake workflows. It sets up a temporary environment, configures the necessary paths, and runs Snakemake to download the specified datasets.

### Prerequisites

- Python 3.7+
- Conda
- Snakemake

### Usage

#### Command Line Arguments

- `--only_luad`: Optional flag to download only the LUAD dataset.
- `--raw_data_output_path`: Optional argument to specify the path to the raw data output directory.
- `--conda_activate_path`: Optional argument to specify the path to the Conda activate script.

#### Example Commands

1. **Download all datasets**:
   ```sh
   fedpydeseq2-download-data
   ```

2. **Download only the LUAD dataset**:
   ```sh
   fedpydeseq2-download-data --only_luad
   ```

3. **Specify a custom raw data output path**:
   ```sh
   fedpydeseq2-download-data --raw_data_output_path /path/to/raw_data
   ```

4. **Specify a custom Conda activate script path**:
   ```sh
   fedpydeseq2-download-data --conda_activate_path /path/to/conda_activate.sh
   ```

### Script Details


#### `download_data` Function

This function handles the main logic for downloading data.

##### Parameters

- `config_path`: The path to the configuration file.
- `download_data_directory`: The path to the download data directory.
- `raw_data_output_path`: The path to the raw data output directory.
- `snakemake_env_name`: The name of the Snakemake environment.
- `conda_activate_path`: The path to the Conda activate script (optional).

#### `create_conda_env` Function

This function creates a Conda environment based on the provided environment file.

##### Parameters

- `env_file`: The path to the environment file.
- `env_prefix`: The prefix (location) where the Conda environment will be created.

#### `main` Function

This function parses command line arguments and calls the `download_data` function with the appropriate parameters.

### Configuration Files

- `config/config.yaml`: General configuration file for downloading all datasets.
- `config/config_luad.yaml`: Configuration file for downloading only the LUAD dataset.

### Example Workflow

1. **Set up the environment**:
   - The script creates a temporary directory and copies the specified download data directory to it.
   - It reads the configuration file and updates it with the raw data output path.

2. **Create Conda environment**:
   - The script creates a Conda environment using the `snakemake_env.yaml` file.

3. **Run Snakemake**:
   - The script runs Snakemake to download the data, using the specified number of cores and resources.

4. **Clean up**:
   - The script removes the temporary directory after the download is complete.

### Error Handling

- The script prints an error message if any exception occurs during the execution.
- It ensures that the temporary directory is cleaned up even if an error occurs.

### Notes

- Ensure that Conda is installed and properly configured on your system.
- The script assumes that Snakemake is available in the Conda environment specified by `snakemake_env.yaml`.


## References

The data downloaded here has mainly been obtained from TCGA and processed by the following
works.

[1] Aran D, Sirota M, Butte AJ.
        Systematic pan-cancer analysis of tumour purity.
        Nat Commun. 2015 Dec 4;6:8971.
        doi: 10.1038/ncomms9971.
        Erratum in: Nat Commun. 2016 Feb 05;7:10707.
        doi: 10.1038/ncomms10707.
        PMID: 26634437; PMCID: PMC4671203.
        <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/>

[2] Jianfang Liu, Tara Lichtenberg, Katherine A. Hoadley, Laila M. Poisson, Alexander J. Lazar, Andrew D. Cherniack, Albert J. Kovatich, Christopher C. Benz, Douglas A. Levine, Adrian V. Lee, Larsson Omberg, Denise M. Wolf, Craig D. Shriver, Vesteinn Thorsson et al.
    An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics,
    Cell,
    Volume 173, Issue 2, 2018, Pages 400-416.e11,
    ISSN 0092-8674,
    <https://doi.org/10.1016/j.cell.2018.02.05>
    <https://www.sciencedirect.com/science/article/pii/S0092867418302290>

[3] Wilks C, Zheng SC, Chen FY, Charles R, Solomon B, Ling JP, Imada EL,
        Zhang D, Joseph L, Leek JT, Jaffe AE, Nellore A, Collado-Torres L,
        Hansen KD, Langmead B (2021).
        "recount3: summaries and queries for
        large-scale RNA-seq expression and splicing."
        _Genome Biol_.
        doi:10.1186/s13059-021-02533-6
        <https://doi.org/10.1186/s13059-021-02533-6>,
        <https://doi.org/10.1186/s13059-021-02533-6>.
