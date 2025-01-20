import json
import os
import fnmatch
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


expand_flag= True
expand_len= 10

def main(args):
    input_folder= args.get("input_folder")
    test_file= args.get("test_file")
    application_folder= args.get("application_folder") # path to leukemia
    samples_file= args.get("samples_file") # path to samples in leukemia used for ISP
    model_folder= args.get("model_folder")
    output_folder= args.get("output_folder")
    bed_folder= args.get("bed_folder")
    gene_name= args.get("gene_name")
    
    
    bed_file= bed_folder + gene_name + ".bed"
    if not os.path.exists(bed_file):
        print(bed_file + " not found!")#continue
    eqtl= pd.read_csv(bed_file, sep= "\t", header= None)
    eqtl.columns= ["chr", "start", "end"]

    samples= pd.read_csv(samples_file, header= None).to_numpy()
    for sample_ID_1 in samples:
        sample_ID= sample_ID_1[0]
        model_name= model_folder + "/" + fnmatch.filter(os.listdir(model_folder), "*" + gene_name + "*")[0]
        if os.path.exists(model_name):
            model= keras.models.load_model(model_name)
            file_name= input_folder + gene_name + "_w1MB_100bs_BinnedActivity.txt.gz"
            data = pd.read_csv(file_name, sep= "\t")
            col_names= data.columns[1:-1]
            data= data.to_numpy()
            test_cells= pd.read_csv(test_file).iloc[:, 0]
            
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
            file_name= application_folder + gene_name + "_w1MB_100bs_leukemia_BinnedActivity.txt.gz"
            data = pd.read_csv(file_name, sep= "\t")
            col_names= data.columns[1:-1]
            data= data.to_numpy()
            x = np.log2((1 + data[:, 1:(data.shape[1])].astype(float)))
            row_names_x = data[:, 0]
            
            chrs, starts, ends, predicteds, ISMs= [], [], [], [], []
            file_flag= False
            for reg in range(0, eqtl.shape[0]):
                #print(reg)
                ##Here I need to implement the ISM based on the start and ed regions
                bin_hits= []
                for i in range(len(col_names)):
                    chr, start, end= col_names[i].split("-")
                    if ((int(start) >= eqtl["start"].iloc[reg]) & (int(end) <= eqtl["end"].iloc[reg]) | (int(start) <= eqtl["start"].iloc[reg]) & (int(end) >= eqtl["start"].iloc[reg]) | (int(start) <= eqtl["end"].iloc[reg]) & (int(end) >= eqtl["end"].iloc[reg])):
                        #print(str(i))
                        bin_hits.append(i)            
                if len(bin_hits) > 0:
                    file_flag= True
                    sample_hit= np.where(row_names_x == sample_ID)[0]
                    sample_data= x[sample_hit, :]
                    if expand_flag:
                        offset= expand_len - len(bin_hits)
                        if offset > 0:
                            bin_hit_0= bin_hits[0]
                            bin_hit_last= bin_hits[len(bin_hits)-1]
                            bin_hits.extend(list(np.arange(bin_hit_0 - math.floor(offset/2), bin_hit_0)))
                            bin_hits.extend(list(np.arange(bin_hit_last, bin_hit_last + math.ceil(offset/2))))
                            valid_idx= np.where((np.asarray(bin_hits) >= 0) & (np.asarray(bin_hits) < x.shape[1]))[0]
                            bin_hits= np.array(bin_hits)[valid_idx]
                    predicted= (model.predict(sample_data, verbose= 0)[0][0] * (mx - mn) + mn) ** 2 - 1
                    if predicted < 0:
                        predicted= 0
                    x_train_ISM= sample_data
                    x_train_ISM[:, bin_hits]= 0
                    ISM= (model.predict(x_train_ISM, verbose= 0)[0][0] * (mx - mn) + mn) ** 2 - 1
                    if ISM < 0:
                        ISM= 0
                    chrs.append(chr)
                    starts.append(eqtl["start"].iloc[reg])
                    ends.append(eqtl["end"].iloc[reg])
                    predicteds.append(predicted)
                    ISMs.append(ISM)
                    dict= {'chr': chrs, 'start': starts, 'end': ends, 'predicted': predicteds, 'ISM': ISMs}
                    df_res = pd.DataFrame(dict)
            if file_flag:
                output_filename= output_folder + '/' + gene_name + '.txt'
                if not os.path.exists(output_folder):
                    os.makedirs(output_folder, exist_ok=True)
                df_res.to_csv(output_filename, mode= 'a', index= False, header= True, sep= "\t")
                subprocess.call("gzip -f " + output_filename, shell=True)

        else:
            print("File not found for gene: " + gene_name)
    print(gene_name + " done!")
    
if __name__ == "__main__":
    json_file= sys.argv[1]
    try:
        with open(json_file, 'r') as file:
            args = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)
         
main(args)




