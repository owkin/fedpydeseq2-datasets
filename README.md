# Datasets organisation

This directory contains the data, assets and scripts necessary to:
- download the raw data necessary to run the tests and experiments, when not
    available in the repository, in the `download_data` directory;
- open the data when performing a Substra experiment in the `assets` directory;
- store the data in the `data` directory.

## Data download

For a detailed description of the data download process, please refer to the
[README](fedpydeseq2_datasets/download_data/README.md).

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


**Origin of the data**
- The `Counts_raw.parquet` and `recount3_metadata.csv` files are downloaded from
    the [RECOUNT3](https://rna.recount.bio/) database.
- The `tumor_purity_metadata.csv` file is downloaded from the
    [Systematic pan-cancer analysis of tumour purity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/) paper.
- The `cleaned_clinical_metadata.csv` file is downloaded from the
    [An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics](https://www.sciencedirect.com/science/article/pii/S0092867418302290#app2) paper.

For more detailed references, see the [References](#references) section.
## Assets

The `assets` directory contains a TCGA opener necessary to open the data on each center
when performing a federated experiment with [Substra](https://docs.substra.org/en/stable/).

In particular, the `fedpydeseq2_datasets/assets/tcga` directory contains the following files:
```
assets/tcga
├── description.md
├── opener.py
```
The opener is a Python script that opens the data and makes it available to the
Substra platform. The `description.md` file contains a description of the data.

For more details on how the opener works, please refer to
the [Substra documentation](https://docs.substra.org/en/stable/).

## Raw data organisation

The `data` directory contains the raw data.
The `raw` directory contains the data downloaded from the original sources,
with the `download_data` scripts.

It is organized as follows:
```
data
├── raw
│   └── tcga
│       ├── <COHORT NAME>
│       │   ├── Counts_raw.parquet
│       │   └── recount3_metadata.csv
│       ├── centers.csv
│       ├── cleaned_clinical_metadata.csv
│       └── tumor_purity_metadata.csv

```

## Data preprocessing

This module not only provides the raw data on which to test `fedpydeseq2`; it also provides the necessary
preprocessing functions, to organise the data according to their center of origin, and to aggregate the raw
data into metadata and counts acceptable to run pydeseq2 or fedpydeseq2.

These preprocessing function usually create the preprocessed data in a `processed_data_path` directory, with the following
structure (the files shown below are created by different preprocessing functions).

```
└── <processed_data_path>
    ├── tcga
    │   ├── <COHORT NAME>
    │   │   ├── counts.parquet
    │   │   └── clinical_data.csv
    ├── centers_data
    │   └── tcga
    │       ├── <dataset_and_dge_experiment_id>
    │       │   ├── <center_i>
    │       │   │   ├── counts_data.csv
    │       │   │   ├── metadata.csv
    │       │   │   └── ground_truth_dds.pkl
    └── pooled_data
        └── tcga
            ├── <dataset_and_dge_experiment_id>
            │   ├── counts_data.csv
            │   ├── metadata.csv
            │   └── ground_truth_dds.pkl
```

These files are automatically generated, if they are not
already present from the raw files or if the `force` option is on.


Note that the centers are always indexed by an integer, starting from 0. For example, one would
have `center_0,...,center_3` if there are
4 centers in the experiment.

The `<dataset_and_dge_experiment_id>` is an identifier of a differential gene expression task (and its specific hyperparameters)
and a TCGA dataset. `fedpydeseq2` or `deseq2` can then be run on the data corresponding to that DGE task.

#### Details on the processed data

In this repository, we study the following cofactors:
- the `gender` of the patients, which is obtained from the `cleaned_clinical_metadata.csv` file;
- the `CPE` of the samples, which is obtained from the `tumor_purity_metadata.csv` file;
- the `stage` of the patients, which is obtained from the `cleaned_clinical_metadata.csv` file.
this stage is originally a stage between `I` and `IV`, but we have grouped them into `I-II-III`
(`Non-advanced`) and `IV` (`Advanced`) stages, to have a binary covariate.
- the `center_id` of the samples, which is obtained from the `centers.csv` file and used
to create natural centers for the federated experiments.


The processing is done by functions in the `fedpydeseq2_datasets` directory. There are three main functions.

- the `common_preprocessing_tcga` function in the `fedpydeseq2_datasets/common_preprocessing.py`
file;
- the `setup_tcga_dataset` function in the `fedpydeseq2_datasets/process_and_split_data.py` file;
- the `setup_tcga_ground_truth_dds` function in the `tcga_preprocessing/create_reference_dds.py` file.


The role of the `common_preprocessing_tcga` function is to generate counts and processed
clinical data for a given cohort (e.g. `LUAD`), from the raw data.
```
└── processed
    ├── tcga
    │   ├── <COHORT NAME>
    │   │   ├── counts.parquet
    │   │   └── clinical_data.csv
```
The `counts.parquet` file contains the counts data, indexed by TCGA sample barcode,
and with columns corresponding to the gene_id in ENSEMBL convention.
Note that we filter out the `PAR_Y` genes, as they are not common to all patients.
The `clinical_data.csv` file aggregates the different metadata from the different sources
described above in a per-cohort fashion. This `csv` is indexed by the TCGA sample barcode.
It contains the following columns:
- `gender`: the gender of the patient;
- `CPE`: CPE stands for "consensus measurement of purity estimations", and is an
aggregate of different purity estimations for the sample;
- `stage`: the stage of the patient, as an integer between 1 and 4
- `center_id`: the center id of the sample, as an integer
- `is_normal_tissue`: a boolean indicating if the sample is a normal tissue or not.
- `T` : the tumor grade of the patient, as an integer between 1 and 4
- `N` : the nodal status of the patient, as an integer between 0 and 3
- `M` : the metastasis status of the patient, as an integer between 0 and 1


The role of the `setup_tcga_dataset` function and the `setup_tcga_ground_truth_dds` function
is to generate the data necessary for the federated AND corresponding pooled experiments, creating
this part of the arborescence:
```
└── processed
    ├── centers_data
    │   └── tcga
    │       ├── <dataset_and_dge_experiment_id>
    │       │   ├── <center_i>
    │       │   │   ├── counts_data.csv
    │       │   │   ├── metadata.csv
    │       │   │   └── ground_truth_dds.pkl
    └── pooled_data
        └── tcga
            ├── <dataset_and_dge_experiment_id>
            │   ├── counts_data.csv
            │   ├── metadata.csv
            │   └── ground_truth_dds.pkl
```
The `<dataset_and_dge_experiment_id>` identifies an experiment. It concatenates
not only the dataset name (TCGA cohort), but also the design factors, continuous factors
as well as other parameters used to filter the data.
The `counts_data.csv` file contains the counts data, indexed by TCGA sample barcode,
and with columns corresponding to the gene_id in ENSEMBL convention.
The `metadata.csv` file contains the clinical data, indexed by the TCGA sample barcode, and
containing only the columns corresponding to a design factor.
The `ground_truth_dds.pkl` file contains the ground truth for the differential expression
analysis, as a `dds` object from the `DESeq2` package.

For more details on these functions, please refer to their respective documentations.

> **Note**: the `setup_tcga_dataset` function
> will binarize the `stage` into two categories: `Advanced` and `Non-advanced`.
> `Advanced` corresponds to stage `IV`, and `Non-advanced` corresponds to stages `I`, `II` and `III`.
> For the TCGA-PRAD cohort, we do not have the stage information, but we infer the stage
> from the `T`, `N` and `M` columns. If the `N` or `M` columns are > 0, the stage is IV (see the
> following [reference](https://www.cancer.org/cancer/types/prostate-cancer/detection-diagnosis-staging/staging.html))
> and hence the `Advanced` stage. Otherwise, it is `Non-advanced`.
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
