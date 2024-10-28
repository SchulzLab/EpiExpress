
---

### R Script for CRE-RF Predictions (CRE_predict_RF.R)

The script below will load the model, normalize the test data, make predictions, and save the results based on the JSON configuration.

```r
# CRE_predict_RF.R

# Load required packages
library(dplyr)
library(randomForest)
library(jsonlite)

# Function to normalize test data using pre-saved min and max values
normalize_test_data <- function(data, min_values, max_values) {
  normalized_data <- data.frame(
    lapply(names(data), function(col) {
      min_val <- min_values[col]
      max_val <- max_values[col]
      
      if (max_val == min_val) {
        return(rep(0, length(data[[col]])))
      } else {
        (data[[col]] - min_val) / (max_val - min_val)
      }
    })
  )
  colnames(normalized_data) <- colnames(data)
  return(normalized_data)
}

# Load configuration from JSON file
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1]
config <- fromJSON(config_path)

models_path <- config$models_path
input_path <- config$input_path
output_path <- config$output_path
min_max_file <- config$min_max_file

# Load min and max values for normalization
min_max_values <- readRDS(min_max_file)

# Loop through each gene input file
input_files <- list.files(input_path, pattern = "\\.txt.gz$", full.names = TRUE)
for (file in input_files) {
  gene_name <- tools::file_path_sans_ext(basename(file))
  
  # Load input data
  df <- read.table(file, header = TRUE, sep = "\t")
  df <- df[order(df$Sample),]
  rownames(df) <- NULL
  df <- df[, colSums(df != 0) > 0]
  
  # Normalize data using stored min/max values
  df_data <- log2(df[,-1] + 1)
  normalized_data <- normalize_test_data(df_data, min_values = min_max_values[[gene_name]]$min, max_values = min_max_values[[gene_name]]$max)
  
  # Load model
  model_file <- file.path(models_path, paste0(gene_name, ".RDS"))
  if (file.exists(model_file)) {
    rf_model <- readRDS(model_file)
    
    # Predict expression values
    predicted_values <- predict(rf_model, normalized_data)
    
    # Rescale predicted values
    min_origin <- min_max_values[[gene_name]]$target_min
    max_origin <- min_max_values[[gene_name]]$target_max
    rescaled_predictions <- (predicted_values * (max_origin - min_origin)) + min_origin
    
    # Save results
    output_file <- file.path(output_path, paste0(gene_name, "_predictions.txt"))
    write.table(data.frame(Sample = df$Sample, Prediction = rescaled_predictions), output_file, row.names = FALSE, sep = "\t")
    print(paste("Predictions saved for", gene_name))
  } else {
    print(paste("Model for gene", gene_name, "not found."))
  }
}
