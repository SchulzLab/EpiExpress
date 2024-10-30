
# Gene Expression Prediction with Pre-Trained Models

This repository provides an R script to predict gene expression levels using pre-trained Random Forest (CRE-RF) models and Python script to predict gene expression levels using pre-trained CNN (Binned-CNN) models. The models are stored on Zenodo, and users can download, unzip, and use them for prediction on their own data.

## Table of Contents
  - [Overview](#overview)
  - [Requirements](#requirements)
  - [Usage](#usage)
    - [1. Download Pre-Trained Models from Zenodo](#1-download-pre-trained-models-from-zenodo)
    - [2. Prepare Input Data](#2-prepare-input-data)
    - [3. Clone the repository](#3-Clone-the-repository)
    - [4. Prediction](#4-prediction)
  - [Notes](#notes)
  - [License](#license)




## Overview

This repository allows users to:
- Download pre-trained models from Zenodo.
- Build your own input data for prediction.
- Predict gene expression using the pre-trained models.
- Generating the results (predicted expression).

The CRE-RF models are built using Random Forest (RF) in R for ~ 28,000 genes, and the predictions are based on a user's dataset containing the same features.

The Binned-CNN models are built using Keras in python for ~ 28,000 genes, and the predictions are based on a user's dataset containing the same feature setup (Fatemeh).

## Requirements 
#### CRE-RF Requirements
 - R version 4.0 or higher
 - R libraries: `randomForest`, `dplyr`, `utils`, `jsonlite`

 Install these libraries in R if you don’t have them already:

 ```r
 install.packages(c("randomForest", "dplyr", "utils", "jsonlite"))
 ```
#### Binned-CNN Requirements
- Fatemeh will add this part

## Usage

### 1. Download Pre-Trained Models from Zenodo

Before running the predictions, you need to download the pre-trained models from Zenodo. The models for each methods are provided as a tar file containing ~ 28,000 models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/uploads/13992024).
  - If you would like to use CRE-RF method, download CRE_RF_models.tar and RF_min_max.tar files.
  - If you would like to use Binned-CNN method, download Binned_CNN_models.tar and CNN_min_max.tar files.
- Uncompress the files into one directory on your local machine.
- you will have to add this path in the JSON file in the next steps.


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

### 3. Clone the repository 

```bash

git clone https://github.com/SchulzLab/ExpressionPredictionModels.git

```

### 4. Prediction
#### Run the script by specifying the path to your configuration JSON file:
To predict the expression by CRE-RF, run this command:

```bash
Rscript CRE-RF/CRE_Predict_RF.R path/to/config.json

```

To predict the expression by Binned-CNN, run this command:
```bash
python Binned_CNN_prediction/CNN_prediction.py path/to/config.json

```


### Output

For each gene, the script will generate a CSV file with the prediction results. 

## Notes

- Ensure that the gene expression files and the model names match. For example, if your gene file is named `ENSG00000123456.txt.gz`, the corresponding model should be `ENSG00000123456.RDS`.


## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

