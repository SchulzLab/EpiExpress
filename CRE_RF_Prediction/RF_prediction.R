# Load necessary libraries
library(jsonlite)
library(randomForest)
library(dplyr)
library(utils)

# Load the JSON file path from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]  # User-specified JSON file path
config <- fromJSON(config_file)

# Extract paths from JSON
model_dir <- config$model_path               # Directory containing model files and min-max files
input_dir <- config$out_folder               # Directory for input gene files (txt.gz) and output predictions

# Function to scale data based on min-max values for each column
scale_data <- function(data, min_values, max_values) {
  scaled_data <- sweep(data, 2, min_values, FUN = "-")
  scaled_data <- sweep(scaled_data, 2, max_values - min_values, FUN = "/")
  return(scaled_data)
}

# Load min-max values from a single file per gene
load_min_max <- function(gene_name) {
  min_max_file <- file.path(model_dir, paste0(gene_name, "_min_max.txt"))
  
  # Read min-max values for each column 
  min_max_data <- read.table(min_max_file, header = TRUE, sep = "\t")  # Ensure the correct separator
  
  # Check if the required columns are present
  if (!all(c("Feature", "Min", "Max") %in% colnames(min_max_data))) {
    stop("Min-max file must contain 'Feature', 'Min', and 'Max' columns.")
  }
  
  # Extract min and max values for each feature
  min_values <- min_max_data$Min
  max_values <- min_max_data$Max
  
  return(list(min = min_values, max = max_values))
}

# Loop through all gene input files in input_dir with .txt.gz extension
input_files <- list.files(input_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

for (input_file in input_files) {
  # Extract gene name from file name
  gene_name <- sub("\\.txt\\.gz$", "", basename(input_file))
  
  # Load model for the gene
  model_file <- file.path(model_dir, paste0(gene_name, ".RDS"))
  if (!file.exists(model_file)) {
    message(paste("Model not found for gene:", gene_name))
    next
  }
  model <- readRDS(model_file)
  
  # Load input data
  data <- read.table(input_file, header = TRUE)
  
  # Apply log transformation: log2(data + 1)
  log_transformed_data <- log2(data + 1)
  
  # Load min-max values for the gene and scale data
  min_max <- load_min_max(gene_name)
  
  # Check if input data columns match feature names in min-max
  if (!all(colnames(data) %in% min_max$Feature)) {
    stop("Some input data columns do not match features in the min-max file.")
  }
  
  # Scale the data
  scaled_data <- scale_data(log_transformed_data, min_max$min, min_max$max)
  
  # Run predictions using the loaded model
  predictions <- predict(model, newdata = scaled_data)
  
  # Save predictions to output_dir with the same gene name, as .csv
  output_file <- file.path(input_dir, paste0(gene_name, "_predictions.csv"))
  write.csv(predictions, output_file, row.names = FALSE)
  
  message(paste("Predictions saved for gene:", gene_name, "in", output_file))
}

message("Prediction process completed.")
