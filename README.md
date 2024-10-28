
# Gene Expression Prediction with Pre-Trained Models

This repository provides an R script to predict gene expression levels using pre-trained Random Forest (CRE-RF) models and Python script to predict gene expression levels using pre-trained CNN (Binned-CNN) models. The models are stored on Zenodo, and users can download, unzip, and use them for prediction on their own data.

## Table of Contents
- [CRE-RF](#CRE-RF)
  - [Overview](#overview)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Usage](#usage)
    - [1. Download Pre-Trained Models from Zenodo](#1-download-pre-trained-models-from-zenodo)
    - [2. Prepare Input Data](#2-prepare-input-data)
    - [3. Run Predictions](#3-run-predictions)
  - [Notes](#notes)
  - [License](#license)
- [Binned-CNN](#Binned-CNN)

## CRE-RF

## Overview

This repository allows users to:
- Download pre-trained gene expression models from Zenodo.
- Load your own data for prediction.
- Predict gene expression using the pre-trained models.
- Generate results, including Mean Squared Error (MSE) and correlation, for the predicted and actual values.

The models are built using Random Forest (RF) for ~ 28,000 genes, and the predictions are based on a user's dataset containing the same features.

## Requirements

- R version 4.0 or higher
- R libraries: `randomForest`, `dplyr`, `utils`

Install these libraries in R if you don’t have them already:

```r
install.packages(c("randomForest", "dplyr"))
```

## Installation

Download the script files from GitHub.

## Usage

### 1. Download Pre-Trained Models from Zenodo

Before running the predictions, you need to download the pre-trained models from Zenodo. The models are provided as a zip file containing 28,000 Random Forest models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/records/13945441).
- Download the zip file (`models.zip`) containing the pre-trained models.
- Unzip the file into a directory on your local machine.

Alternatively, you can automate the download and unzip process by running the script in R. Update the `zenodo_url`, `download_path`, `unzip_path` in the script:

```r
zenodo_url <- "https://zenodo.org/record/your_record_id/files/models.zip"  # Replace with actual Zenodo URL
download_path <- "/path/to/download/directory"
unzip_path <- "/path/to/unzip/directory"
download_models_from_zenodo(zenodo_url, download_path, unzip_path)
```


### 2. Prepare Input Data

To get gene expression predictions on your own data, we provide scripts to generate the input files
for the models in the right format, both for the CRE- and the Binned-feature setup. Both setups are based on a matrix of samples*genomic regions filled with the H3K7ac ChiP-seq signal. The CRE setup uses the ENCODE CREs within a 1 MB window around a gene's 5'TSS as genomic regions. The Binned setup splits the 1 MB window into consecutive bins of size 100 bp. For your own data, you'll only need H3K27ac
ChIP-seq BigWig-files that should contain the fold-change over the control. The other required data can be downloaded or
are provided by us. You will have to create a JSON-file with the paths as shown below. You can find an [example JSON file here](https://github.com/SchulzLab/ExpressionPredictionModels/blob/main/GenerateInput/Example_InputRun.JSON):

| Parameter       | Description                                                                                                                                                                                                                                                                                                                                                                          |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `bigwigs`       | Bigwigs from which to take the average signal in the region features from.Can either be a list (e.g. ["path1.bw", "path2.bw"]), a folder from which all .bw and .bigwig files will be taken (e.g. "bigwig_folder/"), or a path pattern where all files matching the pattern will be used and the string at the asterisk will be used as sample ID (e.g. "bigwig_folder/sample*.bw"). |
| `mode`          | Which output files to produce. Either 'all' or 'CRE,Binned' to get both feature types, 'CRE' for only CRE-based or 'Binned' for only Binned.                                                                                                                                                                                                                                         |
| `out_folder`    | Output folder to which the files will be written to.                                                                                                                                                                                                                                                                                                                                 |
| `gene_file`     | A file with Ensembl IDs (one per line) for which the output will be generated.                                                                                                                                                                                                                                                                                                       |
| `Gene_CRE_file` | Path to the mapping file "CRE_GenePeakMap_Example.txt.gz", available in the Zenodo repository. Only required for the CRE mode.                                                                                                                                                                                                                                                       |
| `gtf_file`      | Path to the hg38 gtf-file, available at [GENCODE](https://www.gencodegenes.org/human/release_38.html). It has to be the comprehensive gene annotation with reference chromosomes only (gencode.v38.annotation.gtf.gz) to ensure the same regions per gene. Required for the Binned mode only.                                                                                        |
| `cores`         | Number of cores to use for steps that are parallelized (Default 1).                                                                                                                                                                                                                                                                                                                  |

To test the input generation and see the output files' format, have a look at the [GenerateInput folder](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/GenerateInput).
There, we provide the scripts and a miniature examples of the files. Inside the folder you can call:

```
python3 BigWigToInput.py Example_InputRun.JSON
```

The script BigWigToInput.py writes the output based on the information given by Example_InputRun.JSON, or your own JSON-file. The output produced 
by this example run is also stored in [ExampleOutput](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/GenerateInput/ExampleOutput).

Example file structure for gene expression data:

```
input_data/
├── ENSG00000123456.txt.gz
├── ENSG00000234567.txt.gz
├── ENSG00000345678.txt.gz
└── ...  # One file per gene
```

Each file should contain tab-delimited data with at least two columns:
- `Sample`: The sample IDs.
- Other columns: H3k27ac signal values for the regions.

### 3. Run Predictions

Now you can run the script "CRE_RF_Prediction.R" to predict gene expression using the pre-trained models. You can predict for a single gene or multiple genes.

#### Predict for a Single Gene

To predict the expression of a single gene, use the following function:

```r
predict_gene_expression(
  gene_name = "ENSG00000123456",  # Example gene name
  test_samples_file = "/path/to/test_samples.csv",
  input_data_dir = "/path/to/input_data",  # Directory containing the gene files
  models_dir = "/path/to/unzipped/models",  # Directory containing unzipped models
  output_dir = "/path/to/output"  # Directory to save the results
)
```

#### Predict for Multiple Genes

To predict for multiple genes listed in a file, use the `predict_multiple_genes` function:

```r
predict_multiple_genes(
  gene_list_file = "/path/to/gene_list.txt",  # File with a list of gene names (one gene per line)
  test_samples_file = "/path/to/test_samples.csv",
  input_data_dir = "/path/to/input_data",  # Directory containing the gene files
  models_dir = "/path/to/unzipped/models",  # Directory containing unzipped models
  output_dir = "/path/to/output"  # Directory to save the results
)
```

### Output

For each gene, the script will generate a CSV file with the prediction results. It includes:
- Mean Squared Error (MSE) between the predicted and actual values.
- Correlation between the predicted and actual values.
- Back-scaled MSE after reversing normalization.

The results for multiple genes will be aggregated into a single CSV file.

## Notes

- Ensure that the gene expression files and the model names match. For example, if your gene file is named `ENSG00000123456.txt.gz`, the corresponding model should be `ENSG00000123456.RDS`.
- The input files should be properly formatted with correct sample IDs and H3k27ac values.


## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

## Binned-CNN
This section describes the Binned-CNN approach. Updating!
