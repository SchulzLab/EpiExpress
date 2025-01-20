import json
import pickle
import gzip
import numpy as np
import math
import pandas as pd
import os
import sys
import statistics
import re
import joblib


###
from sklearn.ensemble import RandomForestRegressor
###
from sklearn import preprocessing
from sklearn.metrics import mean_squared_error
###
from scipy import stats

### Set a random seed for the results to be reproducible
import random

################################################################
################################################################



def main(args): 
    ##### Read the data and begin constructing the model
    fn= args.get("input_file")
    test_file= args.get("test_file")
    seed_value= int(args.get("seed_value"))
    output_file= args.get("output_file")
    
    ## Read in the samples that should be used for test set
    test_cells= pd.read_csv(test_file).iloc[:, 0]

    random.seed(seed_value)


    data = pd.read_csv(fn, sep= "\t").to_numpy()
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

    ########################################################################################################################################

    mn = min(y_train)
    mx = max(y_train)
    y_train = (y_train - mn) / (mx - mn);

    regr_rf = RandomForestRegressor(max_features= "sqrt", random_state= seed_value, oob_score= True)
    regr_rf.fit(x_train, y_train)
    y_pred = regr_rf.predict(x_test)

    y_pred_train = regr_rf.predict(x_train)


    ## training error and cor
    train_cor = np.corrcoef(y_pred_train.squeeze(), y_train)[0, 1].round(2)
    y_pred_train = y_pred_train * (mx - mn) + mn
    y_train = y_train * (mx - mn) + mn
    train_err = mean_squared_error(y_pred_train, y_train).round(3)
    ## test error and cor
    test_cor = np.corrcoef(y_pred.squeeze(),y_test)[0, 1].round(2)
    test_spcor= stats.spearmanr(y_pred.squeeze(),y_test).statistic

    y_pred = y_pred * (mx - mn) + mn
    test_err = mean_squared_error(y_pred.squeeze(), y_test).round(3)

    print("RF test correlation")
    print(test_cor)

    if not os.path.exists(output_file):
        os.makedirs(output_file, exist_ok=True)
    joblib.dump(regr_rf, output_file + '/RF_model.joblib')

    df_cor = pd.DataFrame([train_cor, test_cor, test_spcor])
    df_cor.index= ["RF_train_cor", "RF_test_Pearson_cor", "RF_test_Spearman_cor"]
    df_err = pd.DataFrame([train_err, test_err])
    df_err.index= ["RF_train_err", "RF_test_err"]

    df_cor.to_csv(output_file + '/correlation_summary_RF.csv', sep= '\t')
    df_err.to_csv(output_file+ '/mse_summary_RF.csv', sep= '\t')


if __name__ == "__main__":
    json_file= sys.argv[1]
    try:
        with open(json_file, 'r') as file:
            args = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)
         
main(args)