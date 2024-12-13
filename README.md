
# Epiexpress (Gene Expression Prediction)

1. **`Prediction`**: Contains the code for users to predict gene expression levels using pre-trained Random Forest (CRE-RF) models and Python scripts for Binned-CNN models. Detailed instructions for how users can predict gene expression with their own data are provided below.

## Prediction: Table of Contents
  - [Gene Expression Prediction with Pre-Trained Models](#gene-expression-prediction-with-pre-trained-models)
  - [Usage](#usage)
    - [1. Download Pre-Trained Models from Zenodo](#1-download-pre-trained-models-from-zenodo)
    - [2. Software and packages](#2-software-and-packages)
    - [3. Generate Input Data](#3-generate-input-data)
    - [4. Expression Prediction](#4-expression-prediction)
  - [License](#license)


## Gene Expression Prediction with Pre-Trained Models

This repository allows to:
- Download pre-trained models from Zenodo.
- Build input matrices from our own data.
- Predict gene expression using the pre-trained models.

We provide two types of models: CRE-RF and Binned-CNN. They differ in the feature setup and — as the names suggests — in the type of model used for training. The CRE setup (used for CRE-RF) uses the ENCODE CREs within a 1 MB window around a gene's 5'TSS as features. The Binned setup (for Binned-CNN) splits the 1 MB window 
into consecutive bins of size 100 bp. Both types of approaches produced models for ~28,000 genes each. 

## Usage

### 1. Download Pre-Trained Models from Zenodo

Before running the predictions, you need to download the pre-trained models from Zenodo. The models for each method are provided as a tar file containing ~ 28,000 models, one for each gene.

**Steps:**
- Go to the [Zenodo record](https://zenodo.org/uploads/13992024).
  - If you would like to use the CRE-RF models, download CRE_RF_models.tar.
  - If you would like to use the Binned-CNN method, download CNN_best_models.tar.gz.
- Uncompress the files into one directory on your local machine.

### 2. Software and packages

There are three ways to get the required packages:
1. Create a **conda environment** with the YAML-file we provide. In theory, it should give you all the correct versions (if conda is in your favour), but system dependencies can cause issues.
2. **Manually** install the packages [listed here under dependencies](https://github.com/SchulzLab/EpiExpress/blob/main/Prediction/condaenv.yml).
3. Pull our **Docker** image and work within a Docker container. Docker promises complete reproducibility of the software, but requires the installation of Docker and is less convenient when it comes to setting paths.

#### Conda environment
If you decided to try the installation with conda, clone the repository and set your working directory to Prediction/.
```bash
git clone https://github.com/SchulzLab/ExpressionPredictionModels.git
```
To create the environment run the following:
```bash
conda env create -f condaenv.yml
conda activate epiexpress
```
Depending on the system, the creation command might fail. In this case, you might want to try:
```bash
conda env create --platform osx-64 -f condaenv.yml
```

If you have [mamba](https://mamba.readthedocs.io/en/latest/) installed (same as conda but way faster), replace _conda_ with _mamba_.
When this was successful, you can move to point 3. Generate Input Data.

#### Docker container
As alternative to conda or manual installation, you can download our Docker image, create a container from it and work within that container. It requires a running installation of Docker. To pull the image, run:
```bash
docker image pull dennisheck/epiexpress
```
When starting the container, you will want to mount your local directory in which the downloaded models and your other data is located. Otherwise, the container doesn't have access to your local file system and changes made within the container won't be reflected in your file system (that's the not-so-convenient part mentioned before). 
```
docker run -itd --mount type=bind,source=<absolute local path>,target=</folder name for inside the container> --name epiexpress dennisheck/epiexpress bash
```
If this successfully created the container, you can start an interactive bash window in it with:
```bash
docker start -ai epiexpress
```
To exit the container again, call:
```
exit
```
Please note that when working with the Docker container, all provided paths must match the paths in the container, not those in your local file system. For example, if you have a file _/User/birds/pelican_run.JSON_ and you ran the docker run command with _--mount type=bind,source=/Users/birds/,target=/container_birds_, then the file is found in the container under _container_birds/pelican_run.JSON_. 

### 3. Generate Input Data

To get gene expression predictions on your own data, you do not only need to download the pre-trained models (step above), but also generate the input matrices in the right format. For that, you will need:
- BigWig-files for each sample with H3K27ac ChIP-seq in hg38 that contain the fold-change over the control.
- To download the folder Provided_Input/ with the other required data we provide on [Zenodo](https://zenodo.org/uploads/13992024).
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
Both setups are based on a matrix of samples*genomic regions filled with the H3K27ac ChiP-seq signal. To test the input generation and the subsequent expression prediction, 
you can follow the commands explained here, which will use a small example data set. The entire output of this example
run can also be found in the [ExampleOutput folder](https://github.com/SchulzLab/ExpressionPredictionModels/tree/main/Prediction/ExampleOutput).
For your own data, exchange the example JSON file with your own. If you cloned the directory, move into the Prediction/ folder. Inside the folder you can call:

```bash
python3 GenerateInput/BigWigToInput.py Example_Run.JSON
```

The script BigWigToInput.py writes the output based on the information given by Example_Run.JSON, or your own JSON-file.
Once the input matrices are generated, you can proceed with the next step and get the models' expression prediction. 
If the input matrix for a gene for a feature setup was not written, it is listed in the file _FailedGenes.txt_ with a note saying why.


### 4. Expression Prediction
This step uses the same JSON file as the previous input generation step.
To predict the expression by CRE-RF, run this command:

```bash
Rscript Prediction_CRE_RF/Prediction_CRE_RF.R Example_Run.JSON
```

This will generate a sub-folder CRE_RF_result/ inside the out_folder you declared in the JSON file. For each gene, 
one csv-file will be written with a predicted expression value for each sample found in the bigwigs folder.

In the same manner, you can get the predictions from the Binned-CNN:
```bash
python3 Prediction_Binned_CNN/Prediction_Binned_CNN.py Example_Run.JSON
```

## Issues

If you face issues with running any of the steps, please open an issue here on GitHub and we can work on a solution.

## manuscript codes

This repository also contains another folder named "manuscript_code" which Contains the code and analysis we performed during the development of the project and the manuscript. This includes training models, generating plots, and conducting the analysis presented in the paper. Please see the README file in the folder to get a description of the files. Please note, the scripts are meant as documentation and less as ready-to-use software. 


## License

This project is licensed under the MIT License - see the [LICENSE] file for details.

