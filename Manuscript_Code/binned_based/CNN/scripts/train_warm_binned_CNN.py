import os
import sys
import re
import json

import pickle
import gzip
import numpy as np
import math
import pandas as pd
import statistics

from sklearn.metrics import mean_squared_error

import keras

from keras import Input, Model
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv1D, MaxPooling1D, BatchNormalization, Normalization, LeakyReLU

from keras.callbacks import EarlyStopping

from keras import backend as K

### Set a random seed for the results to be reproducible
import random

################################################################
import tensorflow as tf

################################################################
################################################################
################################################################
def main(args): 
    def contracting(lr, stride_1st, pretrained_model):
        
        opt = keras.optimizers.Adam(learning_rate= lr)

        norm_layer = Normalization(axis= 0)
        model = Sequential()
        model.add(Input(shape= (x_train.shape[1], 1)))
        model.add(BatchNormalization(input_shape = (x_train.shape[1], 1)))
        model.add(Conv1D(10, 100, activation= 'sigmoid', input_shape = (x_train.shape[1], 1), strides= stride_1st))
        model.add(Conv1D(1, 100, activation= 'relu', input_shape = (x_train.shape[1], 1), strides= 1))
        model.add(LeakyReLU(alpha= 0.1))

        model.add(Conv1D(5, 50, activation= 'relu', strides= 1))
        model.add(LeakyReLU(alpha= 0.1))
        model.add(Conv1D(5, 10, activation= 'relu', input_shape = (x_train.shape[1], 1), strides= 1))
        model.add(LeakyReLU(alpha= 0.1))
        model.add(Conv1D(1, 10, activation= 'relu', input_shape = (x_train.shape[1], 1), strides= 1))

        model.add(Flatten())
        model.add(Dense(5, activation= 'relu'))
        model.add(Dropout(0.2))
        model.add(Dense(1, activation= 'sigmoid'))
        Dropout(0.1)
        #model.compile(loss="mse", optimizer="adam")
        if x_train.shape[1] == 10000:
            model.set_weights(pretrained_model.get_weights())
        model.compile(loss="mse", optimizer= opt, metrics=['accuracy'])

        model.build(input_shape= (x_train.shape[1], 1))
        model.summary()
        return(model)
    ################################################################
    ################################################################
    def expanding(lr, stride_1st, pretrained_model):
        opt = keras.optimizers.Adam(learning_rate= lr)
        norm_layer = Normalization(axis= 0)
        model = Sequential()
        model.add(Input(shape= (x_train.shape[1], 1)))
        model.add(BatchNormalization(input_shape = (x_train.shape[1], 1)))

        model.add(Conv1D(20, 10, activation= 'sigmoid', strides= stride_1st))
        model.add(Conv1D(10, 10, activation= 'relu', strides= 2))
        model.add(Conv1D(5, 50, activation= 'relu', strides= 2))
        model.add(Conv1D(5, 100, activation= 'relu', strides= 2))

        model.add(Flatten())
        model.add(Dense(5, activation= 'relu'))
        model.add(Dropout(0.2))
        model.add(Dense(1, activation= 'sigmoid'))
        Dropout(0.1)
        if x_train.shape[1] == 10000:
            model.set_weights(pretrained_model.get_weights())
        model.compile(loss="mse", optimizer= opt, metrics=['accuracy'])
        model.build(input_shape= (x_train.shape[1], 1))
        model.summary()
        return(model)
    ################################################################
    ################################################################

    fn= args.get("input_file")
    test_file= args.get("test_file")
    model_type= args.get("model_type")
    seed_value= int(args.get("seed_value"))
    seed_file= args.get("seed_file")
    output_file= args.get("output_file")

    stride_1st= 5
    lr= 0.001
    batch_size= 32

    pattern = '[\w-]+?(?=\.)'
    fn_trim= re.search(pattern, fn).group()
    gene_name= fn_trim.split("_")[0]

    ## Read in the samples that should be used for test set
    test_cells= pd.read_csv(test_file).iloc[:, 0]


    ##### Read the data and begin constructing the model
    random.seed(seed_value)
    tf.random.set_seed(seed_value)

    print("reading warm seed: " + seed_file)
    pretrained_model= keras.models.load_model(seed_file)

    data = pd.read_csv(fn, sep= "\t").to_numpy()
    row_names = data[:, 0]

    y = np.log2((1 + data[:, data.shape[1]-1].astype(float)))
    x = np.log2((1 + data[:, 1:(data.shape[1]-1)].astype(float)))

    ## Split data into training and test sets
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
    epochs= 1000
    callback = EarlyStopping(monitor='loss', patience= 5)
    #    model.summary(print_fn=lambda x: f.write(x + '\n'))
    if model_type == "contracting":
        model= contracting(lr, stride_1st, pretrained_model)
    else:
        model= expanding(lr, stride_1st, pretrained_model)

    print("training " + model_type)
    history= model.fit(x_train, y_train, batch_size= batch_size, epochs= epochs, validation_split= 0.2, callbacks=[callback], verbose= 2) #, verbose= 0)

    y_pred = model.predict(x_test)
    y_pred_train = model.predict(x_train)
    print(model.evaluate(x_train, y_train))

    ## training error and cor
    train_cor = np.corrcoef(y_pred_train.squeeze(), y_train)[0, 1].round(2)
    y_pred_train = y_pred_train * (mx - mn) + mn
    y_train = y_train * (mx - mn) + mn
    train_err = mean_squared_error(y_pred_train, y_train).round(3)
    ## test error and cor
    test_cor = np.corrcoef(y_pred.squeeze(), y_test)[0, 1].round(2)
    ## Backscale the response to make it comparable with ENet
    y_pred = y_pred * (mx - mn) + mn

    test_err = mean_squared_error(y_pred.squeeze(), y_test).round(3)

    print("CNN test correlation")
    print(test_cor)

    if not os.path.exists(output_file):
        os.makedirs(output_file, exist_ok=True)
    model.save(output_file + '/CNN_' + model_type+ '.h5')
    with open(output_file + '/CNN_history_'+ model_type + '.pkl', 'wb') as file_pi:
        pickle.dump(history.history, file_pi)



    df_cor = pd.DataFrame([train_cor, test_cor])
    df_cor.index= ['CNN_train_cor', 'CNN_test_cor']
    df_err = pd.DataFrame([train_err, test_err])
    df_err.index= ['CNN_train_err', 'CNN_test_err']

    df_cor.to_csv(output_file + '/correlation_summary_' + model_type  + '.csv', sep= '\t')
    df_err.to_csv(output_file+ '/mse_summary_' + model_type + '.csv', sep= '\t')
   
if __name__ == "__main__":
    json_file= sys.argv[1]
    try:
        with open(json_file, 'r') as file:
            args = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)
         
main(args)