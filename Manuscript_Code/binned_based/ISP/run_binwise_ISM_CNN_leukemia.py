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

import keras
import shap
################################################################
### Set a random seed for the results to be reproducible
import random
seed_value = 12345
random.seed(seed_value)
gene_name= "ENSG00000001461"

expand_flag= False
expand_len= 10

#fo = open("/projects/apog/work/IHEC/HKandCTGenes.txt", "r")
#genes = fo.readlines()
path = "/projects/apog/work/IHEC/Binned_Activity/w1MB_100bs_cleaned/"
path_leuk = "/projects/apog/work/IHEC/Binned_Activity/w1MB_100bs_leukemia/"


#val_loss= pd.read_csv("/projects/apog/work/CNN_bin/scripts/CNN_validation_df_warmseeds.txt", sep= "\t")
val_loss= pd.read_csv("/projects/apog/work/CNN_bin/scripts/CNN_validation_df_warmseeds_90_percent_training_sigmoid.txt", sep= "\t")


samples= [x.rstrip("\n") for x in open("/projects/apog/work/models/1MB/Luekemia/adult_luekemia_sample.txt", "r").readlines()]
fn= args.get("input_file")
test_file= args.get("test_file")
application_file= args.get("application_file") # path to leukemia
samples_file= args.get("samples_file") # path to samples in leukemia used for ISP
model_file= args.get("model_file")
seed_value= int(args.get("seed_value"))
seed_file= args.get("seed_file")
output_file= args.get("output_file")
        
for sample_ID in samples:
    if os.path.exists(model_name):
        model= keras.models.load_model(model_file)
        data = pd.read_csv(fn, sep= "\t")
        col_names= data.columns[1:-1]
        data= data.to_numpy()
        test_cells= pd.read_csv(test_file).iloc[:, 0]
        samples= [x.rstrip("\n") for x in open(samples_file, "r").readlines()]

        row_names = data[:, 0]
        y = np.log2((1 + data[:, data.shape[1]-1].astype(float)))
        x = np.log2((1 + data[:, 1:(data.shape[1]-1)].astype(float)))

        test_idx = [row_names.tolist().index(x) for x in test_cells]
        train_idx = [idx for idx in range(0, len(y)) if idx not in test_idx]

        row_names_train= row_names[train_idx]
        row_names_test= row_names[test_idx]

        x_train = x[train_idx, :]
        y_train = y[train_idx]

        x_test = x[test_idx, :]
        y_test = y[test_idx]
        ## Min-Max scale the response to make it compatible with the sigmoid activation function (ENet didn't like the min-max scaling, so I'm moving that chucnk of code here to have it applied only to CNN)
        mn = min(y_train)
        mx = max(y_train)
        y_train = (y_train - mn) / (mx - mn);
        
        ## read the application data in
        file_name= application_file + gene_name + ".txt.gz"
        data = pd.read_csv(file_name, sep= "\t")
        col_names= data.columns[1:-1]
        data= data.to_numpy()
        
        row_names = data[:, 0]
        y = np.log2((1 + data[:, data.shape[1]-1].astype(float)))
        x = np.log2((1 + data[:, 1:(data.shape[1]-1)].astype(float)))
        chrs, starts, ends, predicteds, ISMs= [], [], [], [], []
        file_flag= False
        for i in range(len(col_names)):
            print(i)
            chr, start, end= col_names[i].split("-")
            file_flag= True
            sample_hit= np.where(row_names == sample_ID)[0]
            sample_data= x[sample_hit, :]

            predicted= (model.predict(sample_data, verbose= 0)[0][0] * (mx - mn) + mn) ** 2 - 1
            if predicted < 0:
                predicted= 0
            x_train_ISM= sample_data
            x_train_ISM[:, i]= 0
            ISM= (model.predict(x_train_ISM, verbose= 0)[0][0] * (mx - mn) + mn) ** 2 - 1
            if ISM < 0:
                ISM= 0
            #val_res= pd.concat([val_res, pd.DataFrame((chr, crispri_hit["start"].iloc[reg], crispri_hit["end"].iloc[reg], gene_name, predicted, ISM))], axis= 1, ignore_index= True)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            predicteds.append(predicted)
            ISMs.append(ISM)
        dict= {'chr': chrs, 'start': starts, 'end': ends, 'predicted': predicteds, 'ISM': ISMs}
        df_res = pd.DataFrame(dict)
        ratio= np.log2((df_res["ISM"] + 1) / (df_res["predicted"] + 1))
            #my_list= (chr, crispri_hit["start"].iloc[reg], crispri_hit["end"].iloc[reg], gene_name, predicted, ISM)
            #val_res= pd.concat([val_res, pd.DataFrame((chr, crispri_hit["start"].iloc[reg], crispri_hit["end"].iloc[reg], gene_name, predicted, ISM))], axis= 0, ignore_index= True)
        if file_flag: 
            output_filename= output_folder + '/' + gene_name + '.txt'
            f= open(output_filename, 'w')
            f.write('#pearson_r=' + str(cor) + '\n#backscaled_MSE=' + str(err) + '\n')
            f.close()
            df_res.to_csv(output_filename, mode= 'a', index= False, header= True, sep= "\t")
            subprocess.call("gzip -f " + output_filename, shell=True)

    else:
        print("File not found for gene: " + gene_name)

print(gene_name + " done!")






