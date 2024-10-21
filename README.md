
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
  - [File Structure](#file-structure)
  - [Notes](#notes)
  - [License](#license)
- [Binned-CNN](#Binned-CNN)

## CRE-RF

## Overview

This repository allows users to:
- Download pre-trained gene expression models from Zenodo.
- Load their own data for prediction.
- Predict gene expression using the pre-trained models.
- Generate results, including Mean Squared Error (MSE) and correlation, for the predicted and actual values.

The models are built using Random Forest (RF) for 28,000 genes, and the predictions are based on a user's dataset containing the same features.

## Requirements

- R version 4.0 or higher
- R libraries: `randomForest`, `dplyr`, `utils`

Install these libraries in R if you don’t have them already:

```r
install.packages(c("randomForest", "dplyr"))
```

## Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/your-username/gene-expression-prediction.git
cd gene-expression-prediction
```

Alternatively, download the script files from GitHub.

## Usage

### 1. Download Pre-Trained Models from Zenodo

Before running the predictions, you need to download the pre-trained models from Zenodo. The models are provided as a zip file containing 28,000 Random Forest models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/records/13945441).
- Download the zip file (`models.zip`) containing the pre-trained models.
- Unzip the file into a directory on your local machine.

Alternatively, you can automate the download and unzip process by running the script in R. Update the `zenodo_url`, `download_path`, and `unzip_path` in the script:

```r
zenodo_url <- "https://zenodo.org/record/your_record_id/files/models.zip"  # Replace with actual Zenodo URL
download_path <- "/path/to/download/directory"
unzip_path <- "/path/to/unzip/directory"

download_models_from_zenodo(zenodo_url, download_path, unzip_path)
```

### 2. Prepare Input Data

You need to prepare two main input files:
1. **Gene Expression Data**: A file for each gene containing log-transformed expression levels, stored in a directory.
2. **Train/Test Sample Files**: CSV files containing lists of sample IDs for training and testing.

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
- Other columns: Gene expression values for that gene.

### 3. Run Predictions

Now you can run the script "CRE_RF_Prediction.R" to predict gene expression using the pre-trained models. You can predict for a single gene or multiple genes.

#### Predict for a Single Gene

To predict the expression of a single gene, use the following function:

```r
predict_gene_expression(
  gene_name = "ENSG00000123456",  # Example gene name
  train_samples_file = "/path/to/train_samples.csv",
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
  train_samples_file = "/path/to/train_samples.csv",
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

## File Structure

After downloading and unzipping the models, your directory structure should look something like this:

```
gene-expression-prediction/
├── models/  # Directory where models are stored after unzipping
├── input_data/  # Directory with gene expression data
├── train_samples.csv  # List of training sample IDs
├── test_samples.csv  # List of testing sample IDs
└── gene_list.txt  # List of gene names (one per line)
```

## Notes

- Ensure that the gene expression files and the model names match. For example, if your gene file is named `ENSG00000123456.txt.gz`, the corresponding model should be `ENSG00000123456.RDS`.
- The input files should be properly formatted with correct sample IDs and gene expression values.

## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

## Binned-CNN
This section describes the Binned-CNN approach.
