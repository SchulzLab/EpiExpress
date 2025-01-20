import os
import pandas as pd
import numpy as np
import pickle
import zipfile
import subprocess
import gzip

import math
import sys

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from matplotlib.backends.backend_pdf import PdfPages

import plotly.express as px


################################################################
### Set a random seed for the results to be reproducible
import random
seed_value = 12345
random.seed(seed_value)
args= sys.argv
gene_name= args[1]
celltype= args[2]
#gene_name= "ENSG00000154760"

expand_flag= True
expand_len= 10
#gene_name= "ENSG00000060749"# "ENSG00000162927" # "ENSG00000110172" 

#fo = open("/projects/apog/work/IHEC/HKandCTGenes.txt", "r")
#genes = fo.readlines()
path = "/projects/apog/work/IHEC/Binned_Activity/w1MB_100bs_cleaned/"
out_path= "/projects/apog/work/IHEC/" + celltype + "_predictions/CNN/"

if os.path.exists(out_path + "scaled_space/" + gene_name + ".txt"):
    exit("Already exists!")
import keras
import shap
fo = open("/projects/apog/work/CNN_bin/test_cells_cleaned.txt", "r")
test_cells = [x.rstrip("\n") for x in fo.readlines()]

val_loss= pd.read_csv("/projects/apog/work/CNN_bin/scripts/CNN_validation_df_warmseeds_90_percent_training_sigmoid.txt", sep= "\t")
val_loss_gene= val_loss[(val_loss["gene"] == gene_name)].iloc[0]
mtag= val_loss_gene["topology"]
seed_tag= val_loss_gene["seed"]
model_name= "/projects/apog/work/CNN_bin/experiments/models_" + gene_name + "_w1MB_100bs_BinnedActivityindiv_" + mtag + "_part0_" + seed_tag + "_stride_5_0_LR_0.001_bs_32/CNN_" + gene_name + "_w1MB_100bs_BinnedActivityindiv_" + mtag + "_part0_" + seed_tag + "_stride_5_0_LR_0.001_bs_32_w1MB_100bs_" + gene_name + ".h5"

if os.path.exists(model_name):
    model= keras.models.load_model(model_name)
    cor= val_loss_gene["test_cor"] #pd.read_csv(cor_path, sep= "\t").iloc[1, 1]
    err= val_loss_gene["test_err"] #pd.read_csv(err_path, sep= "\t").iloc[1, 1]
    brain_file_name= "/projects/apog/work/IHEC/Binned_Activity/w1MB_100bs_"+ celltype + "/" + gene_name + "_w1MB_100bs_" + celltype + "_BinnedActivity.txt.gz"
    file_name= path + gene_name + "_w1MB_100bs_BinnedActivity.txt.gz"
    data = pd.read_csv(file_name, sep= "\t")
    brain_data = pd.read_csv(brain_file_name, sep= "\t").to_numpy()
    col_names= data.columns[1:-1]
    data= data.to_numpy()
    
    row_names = data[:, 0]
    y = np.log2((1 + data[:, data.shape[1]-1].astype(float)))
    x = np.log2((1 + data[:, 1:(data.shape[1]-1)].astype(float)))
    
    brain_x = np.log2((1 + brain_data[:, 1:brain_data.shape[1]].astype(float)))

    test_idx = [row_names.tolist().index(x) for x in test_cells]
    train_idx = [idx for idx in range(0, len(y)) if idx not in test_idx]

    row_names_train= row_names[train_idx]
    row_names_test= row_names[test_idx]

    x_train = x[train_idx, :]
    y_train = y[train_idx]

    x_test = x[test_idx, :]
    y_test = y[test_idx]
    mn = min(y_train)
    mx = max(y_train)
    y_train = (y_train - mn) / (mx - mn);
    
    predicted= model.predict(brain_x) 
    predicted_bs= (predicted * (mx - mn) + mn) ** 2 - 1
    df= pd.DataFrame(predicted)
    df_bs= pd.DataFrame(predicted_bs)
    
    df.to_csv(out_path + "scaled_space/" + gene_name + ".txt", sep= "\t")
    df_bs.to_csv(out_path + "backscaled_space/" + gene_name + "_bs.txt", sep= "\t")





