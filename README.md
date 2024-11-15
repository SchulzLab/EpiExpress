
# Gene Expression Prediction with Pre-Trained Models

This repository provides two main folders:

1. **`Prediction`**: Contains the code for users to predict gene expression levels using pre-trained Random Forest (CRE-RF) models and Python scripts for Binned-CNN models. Detailed instructions for how users can predict gene expression with their own data are provided below.

2. **`manuscript_code`**: Contains the code and analysis we performed during the development of the project and the manuscript. This includes training models, generating plots, and conducting the analysis presented in the paper. The code in this folder was used for internal purposes, such as in silico perturbation (ISP), figure generation, and model training.

 For users interested in running predictions with pre-trained models, please focus on the **`Prediction`** folder. For those looking for the code related to model development and the 
 analysis presented in the manuscript, refer to the **`manuscript_code`** folder, which has subfolders like `ISP`, `Plots`, and `Training` that organize our project scripts.


 This repository provides an R script to predict gene expression levels using pre-trained Random Forest (CRE-RF) models and Python script to predict gene expression levels using pre- 
 trained CNN (Binned-CNN) models. The models are stored on Zenodo, and users can download, unzip, and use them for prediction on their own data.

## Table of Contents
  - [Overview](#Gene Expression Prediction with Pre-Trained Models)
  - [Requirements](#requirements)
  - [Usage](#usage)
    - [1. Download Pre-Trained Models from Zenodo](#1-download-pre-trained-models-from-zenodo)
    - [2. Prepare Input Data](#2-prepare-input-data)
    - [3. Clone the repository](#3-Clone-the-repository)
    - [4. Prediction](#4-prediction)
  - [Notes](#notes)
  - [License](#license)




## Gene Expression Prediction with Pre-Trained Models

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

Before running the predictions, you need to download the pre-trained models from Zenodo. The models for each method are provided as a tar file containing ~ 28,000 models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/uploads/13992024).
  - If you would like to use the CRE-RF models, download CRE_RF_models.tar.
  - If you would like to use the Binned-CNN method, download CNN_best_models.tar.gz.
- Uncompress the files into one directory on your local machine.
- You will have to add this path in the JSON file in the next steps.

### 2. Clone the repository 

```bash
git clone https://github.com/SchulzLab/ExpressionPredictionModels.git
```

### 3. Generate Input Data

To get gene expression predictions on your own data, you do not only need to download the pre-trained models (step above), but also generate the input matrices in the right format. For that, you will need:
- BigWig-files for each sample with H3K27ac ChIP-seq in hg38 that contain the fold-change over the control.
- To download the folder with the other required data we provide on Zenodo (TODO LINK).
- Prepare a JSON file with the paths and options as explained below. You can find an [example JSON file here](https://github.com/SchulzLab/ExpressionPredictionModels/blob/main/Prediction/Example_Run.JSON).

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

We provide scripts to generate the input files for the models in the right format, both for the CRE- and the Binned-feature setup. 
Both setups are based on a matrix of samples*genomic regions filled with the H3K27ac ChiP-seq signal. The CRE setup uses the 
ENCODE CREs within a 1 MB window around a gene's 5'TSS as genomic regions. The Binned setup splits the 1 MB window 
into consecutive bins of size 100 bp. To test the input generation and the subsequent expression prediction, 
you can follow the commands explained here, which will use a small example data set. The entire output of this example
run can also be found in the [ExampleOutput folder](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/Prediction/ExampleOutput).
For your own data, exchange the example JSON file with your own. If you cloned the directory, move into the Prediction/ folder. Inside the folder you can call:

```bash
python3 GenerateInput/BigWigToInput.py Example_Run.JSON
```

The script BigWigToInput.py writes the output based on the information given by Example_Run.JSON, or your own JSON-file.
Once the input matrices are generated, you can proceed with the next step and get the models' expression prediction. 
If the input matrix for a gene for a feature setup was not written, it is listed in the file _FailedGenes.txt_ with a note saying why.


### 4. Prediction
#### Run the script by specifying the path to your configuration JSON file:
This step uses the same JSON file as the input generation.
To predict the expression by CRE-RF, run this command:

```bash
Rscript Prediction_CRE_RF/Prediction_CRE_RF.R Example_Run.JSON
```

This will generate a sub-folder CRE_RF_result inside the out_folder you declared in the JSON file. For each gene, 
one csv-file will be written with a predicted expression value for each sample found in the bigwigs folder.

In the same manner, you can get the predictions from the Binned-CNN:
```bash
python3 Prediction_Binned_CNN/Prediction_Binned_CNN.py Example_Run.JSON
```


## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

