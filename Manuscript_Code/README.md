# Manuscript Code  

This folder contains all the scripts and analyses developed during this project, including model training, data analysis, and visualization. These scripts document the full workflow used to train models, perform analyses, and generate plots, some of which are included in the manuscript *"Harnessing machine learning models for epigenome to transcriptome association studies."*  


## Folder Structure  

### 1. Training  
Contains scripts used to train different methods for predicting gene expression from H3K27ac data:

- CRE based methods: Random Forest (RF) and Multi Layer Perceptron (MLP).
- Binned based methods: Convolutional Neural Network (CNN) and Random Forest (RF).
- STITCHIT: A linear regression model with elastic net.

#### Key Scripts:  
- `CRE_RF.R`: Code for training RF models using CRE regions from ENCODE.  
- `CRE_NN1.R`: Code for training one hidden layer MLP models using CRE regions from ENCODE.
- `CRE_NN2.R`: Code for training two hidden layers MLP models using CRE regions from ENCODE.
- `Binned_CNN.py`: Code for training CNN models using binned genomic regions.  
- `Binned_RF.py`: Code for training RF models using binned genomic regions.
- `STITCHIT.R` : Code for training STITCHIT models using its own segmentation approach for defining features.

  
---


### 2. ISP (In Silico Perturbation)  
Scripts for performing in silico perturbation (ISP) analyses, which measure the influence of individual features on gene expression predictions by setting features to zero and observing the resulting changes.  

#### Key Scripts:  
- `ISP_CRE_RF_Crispri.R`: Calculates ISP scores using CRE-RF models.  
- `ISP_Binned_CNN_Crispri.py`: Calculates ISP scores using Binned-CNN models.
- `CRISPRi_Intersection.py` with `CRISPRi_Helpers.py`: Collect the ISP for the K562 sample from the different models and intersects the features with experimentally validated CRISPRi-interactions.
- `eQTL_Recall.py` with `eQTL_Helpers.py`: Collect the ISP score across tissue-mapped samples for all models and calculate the recall among the 100,000 highest scored interactions with GSEA.
 
---

### 3. Plots  
Includes scripts for generating visualizations presented in the manuscript, such as model performance comparisons, validation results, and gene characteristics.  

#### Key Scripts:  
- `Dotplot_Overlap_DisgenetGenes_and_SignificantGenes.R`: Generates plots showing the overlap between DisGeNET's leukemia genes and the significant genes identified by our models in the analysis.
- `Scatterplots_logfoldChange_RF_CNN_significantGenes.R`: Generates two scatterplots: one showing the correlation of fold change between leukemia genes and significant genes identified by CRE-RF, and another showing the correlation of fold change for Binned-CNN
- `CRISPRi_Plots.py`: Plots Precision-Recall curves and the resulting AUPRC for model comparison.
- `eQTL_Plots.py`: Plots the recall of recovering eQTL-gene interactions among the top 100,000 interactions of each model.

---

## Purpose  
These scripts were developed as part of the research process and were used internally to:  
1. Train machine learning models.  
2. Analyze the impact of regulatory features on gene expression.  
3. Validate model predictions.  
4. Generate figures and data for the manuscript.  

The prediction functionality has been made user-friendly and is available in the `Prediction` folder. Users interested in predicting gene expression using the pre-trained models should refer to the main repository README file for guidance.  

## Notes  
- These scripts are specific to the data and workflow used in the project.  
- They are not intended for standalone use outside the context of the project’s analyses.  

