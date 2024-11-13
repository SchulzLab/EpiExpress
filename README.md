
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

The Binned-CNN models are built using Keras in python for ~ 28,000 genes, and the predictions are based on a user's dataset containing the same feature setup (need to be checked).

## Requirements 
#### CRE-RF Requirements
 - R version 4.0 or higher
 - R libraries: `randomForest`, `dplyr`, `utils`, `jsonlite`

 Install these libraries in R if you don’t have them already:

 ```r
 install.packages(c("randomForest", "dplyr", "utils", "jsonlite"))
 ```
#### Binned-CNN Requirements
- (need to be added)

## Usage

### 1. Download Pre-Trained Models from Zenodo

Before running the predictions, you need to download the pre-trained models from Zenodo. The models for each methods are provided as a tar file containing ~ 28,000 models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/uploads/13992024).
  - If you would like to use CRE-RF method, download CRE_RF_models.tar and RF_min_max.tar files.
  - If you would like to use Binned-CNN method, download CNN_best_models.tar.gz and CNN_min_max.tar files.
- Uncompress the files into one directory on your local machine.
- you will have to add this path in the JSON file in the next steps.


### 2. Generate Input Data

To get gene expression predictions on your own data, you do not only need to download the pre-trained models (step above), but also generate the input matrices in the right format. For that, you will need:
- BigWig-files of H3K27ac ChIP-seq in hg38 for each sample that contain the fold-change over the control.
- To download the folder with the other required data we provide on Zenodo (TODO LINK).
- Prepare a JSON file with the paths and options as explained below. You can find an [example JSON file here](https://github.com/SchulzLab/ExpressionPredictionModels/blob/main/GenerateInput/Example_InputRun.JSON).

| Parameter            | Description                                                                                                                                                                                                                                                                                                                                                                           |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `bigwigs`            | Bigwigs from which to take the average signal in the region features from. Can either be a list (e.g. ["path1.bw", "path2.bw"]), a folder from which all .bw and .bigwig files will be taken (e.g. "bigwig_folder/"), or a path pattern where all files matching the pattern will be used and the string at the asterisk will be used as sample ID (e.g. "bigwig_folder/sample*.bw"). |
| `mode`               | Which output files to produce. Either 'all' or 'CRE-RF,Binned-CNN' to get both feature types, 'CRE-RF' for only CRE-based or 'Binned-CNN' for only Binned.                                                                                                                                                                                                                            |
| `out_folder`         | Output folder to which the files will be written to. A subfolder will be created for each mode.                                                                                                                                                                                                                                                                                       |
| `gene_file`          | A file with Ensembl IDs (one per line) for which the output will be generated.                                                                                                                                                                                                                                                                                                        |
| `model_folder`       | Path to where the trained models were downloaded in the step above.                                                                                                                                                                                                                                                                                                                   |
| `provided_input`     | Path to the folder 'ProvidedInput' that holds additional files necessary for generating the input. It is available on Zenodo (TODO LINK). Please do not change the file names, otherwise the scripts will crash.                                                                                                                                                                      |
| `cores`              | Number of cores to use for steps that are parallelized (Default 1).                                                                                                                                                                                                                                                                                                                   |
 | `correlation_cutoff` | Only models with a Pearson correlation coefficient between predicted and test data above the cutoff will be used (Default 0).                                                                                                                                                                                                                                                        |

We provide scripts to generate the input files for the models in the right format, both for the CRE- and the Binned-feature setup. Both setups are based on a matrix of samples*genomic regions filled with the H3K27ac ChiP-seq signal. The CRE setup uses the ENCODE CREs within a 1 MB window around a gene's 5'TSS as genomic regions. The Binned setup splits the 1 MB window into consecutive bins of size 100 bp. To test the input generation and see the output files' format, have a look at the [GenerateInput folder](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/GenerateInput).
There, we provide the scripts and a miniature examples of the files. Inside the folder you can call:

```
python3 BigWigToInput.py Example_InputRun.JSON
```

The script BigWigToInput.py writes the output based on the information given by Example_InputRun.JSON, or your own JSON-file. The output produced 
by this example run is also stored in [ExampleOutput](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/GenerateInput/ExampleOutput).
Once the input matrices are generated, you can proceed with the next step and get the models' expression prediction. If the input matrix for a gene for a feature setup was not written, it is listed in the file _FailedGenes.txt_ with a note saying why.

### 3. Clone the repository 

```bash

git clone https://github.com/SchulzLab/ExpressionPredictionModels.git

```

### 4. Prediction
#### Run the script by specifying the path to your configuration JSON file:
 Make sure that you provided the paths to all nedded files on the JSON file:
 
 - "model_folder" should be the path to the model files that you have downloaded from Zenodo.
 - "provided_input" should contain the min-max files that you have downloaded from Zenodo.
 - "out_folder" should be the path to your input data (the one you created in previous step). The final prediction files will be stored in a subfolder named "result".

To predict the expression by CRE-RF, run this command:

```bash
Rscript ExpressionPredictionModels/Prediction/Prediction_CRE_RF/Prediction_CRE_RF.R path/to/config.json

```

To predict the expression by Binned-CNN, run this command:
```bash
python ExpressionPredictionModels/Prediction/Binned_CNN_prediction/Prediction_Binned_CNN.py path/to/config.json

```


### Output

For each gene, the script will generate a CSV file with the prediction results. 

## Notes

- Ensure that the input files and the model names match. For example, if your gene file is named `ENSG00000123456.txt.gz`, the corresponding model should be `ENSG00000123456.RDS`.


## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

