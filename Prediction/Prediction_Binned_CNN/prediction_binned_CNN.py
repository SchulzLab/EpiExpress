import json
import os
import pandas as pd
import numpy as np
import sys

import keras
################################################################
### Set a random seed for the results to be reproducible
import random
seed_value = 12345
random.seed(seed_value)


def main(args):
    input_file= args.get("input_folder")
    model_file= args.get("model_folder")
    output_folder= args.get("output_folder")
    minmax_file= args.get("minmax_file")
    
    if os.path.exists(model_file):
        model= keras.models.load_model(model_file)
        data = pd.read_csv(input_file, sep= "\t")
        data= data.to_numpy()
        row_names = data[:, 0]
        mn, mx= pd.read_csv(minmax_file).to_numpy()
        x= np.log2((1 + data[:, 1:data.shape[1]].astype(float)))

        ## Generate predictions and backscale the values using the given min and max
        prediction= (model.predict(x, verbose= 0) * (mx - mn) + mn) ** 2 - 1
        
        df_res = pd.DataFrame(prediction, index=row_names, columns=['predicted_expr'])
        output_filename= output_folder + '/prediction.txt'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder, exist_ok=True)
        df_res.to_csv(output_filename)

    else:
        print("File: " + model_file+ " not found!")
    print("Prediction done!")
    
if __name__ == "__main__":
    json_file= sys.argv[1]
    try:
        with open(json_file, 'r') as file:
            args = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)
         
main(args)




