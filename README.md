# MorphNet Datasets

This repository contains datasets and codes used in [MorphNet]().

## Getting Started

To use the Patch-seq dataset with MorphGAN, please download the processed zip files from [Google Drive](). The directory structure should look like below, with (*) indicating required folders:

```
patchseq
├── data (* required)
│   ├── processed (* required)
│   │   ├── latent_gene_expressions (* required)
│   │   ├── metadata (* required)
│   │   ├── morphological_images (* required)
│   └── raw
├── notebooks
├── src
├── __init__.py (* required)
└── patchseq.py (* required)
```
## Advanced Usage
It is possible to preprocess the Patch-seq dataset from raw data. First follow the [guide on downloading raw data](#downloading-raw-data) to obtain the raw dataset. We uploaded the code to preprocess the raw dataset under the `src/` folder.

## Downloading Raw Data
#### Tolias Dataset
[Link to Paper](https://www.nature.com/articles/s41586-020-2907-3#data-availability)

##### Downloading Raw Data
Link to download the `tolias` dataset can be found at [berenslab/mini-atlas](https://github.com/berenslab/mini-atlas) Github repo. For MorphGAN, only the following data needs to be downloaded:

1. Raw morphological reconstructions in SWC file format
2. Gene Experssion dataset
    - mini-atlas/data/m1_patchseq_exon_counts.csv.gz
    - mini-atlas/data/m1-patchseq_intron_counts.csv.gz
3. Metadata
    - mini-atlas/data/m1_patchseq_meta_data.csv

#### Allen Dataset

##### Downloading Raw Data

The morphological reconstruction files are hosted in the Brain Image Library (BIL), which may be downloaded using HTTP, SFTP, or Globus. Further information on each method is available at https://www.brainimagelibrary.org/download.html.

- The link to morphological reconstruction files is: ftp://download.brainlib.org:8811/biccn/zeng/pseq/morph/200526/

```
wget -r ftp://download.brainlib.org:8811/biccn/zeng/pseq/morph/200526/
```

Transcriptomics and metadata are hosted at the Neuroscience Multi-omic Data Archive (NeMO).

- The link to gene expression and metadata is: http://data.nemoarchive.org/other/grant/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/

#### Macosko 10x Nuclei V3 Dataset

This folder have been downloaded from UMICH Box owned by Dr. Joshua Welch.

Path: `IterativeRefinement > macosko_10x_nuclei_v3`
Files:
- cluster.annotaiton.csv
- cluster.membership.csv
- matrix.RDS
- neurons_passing_qc.RDS
