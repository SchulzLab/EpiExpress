# Load necessary libraries
library(jsonlite)
library(randomForest)
library(dplyr)
library(utils)

print("reading the json file")
# Load the JSON file path from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]  # User-specified JSON file path
config <- fromJSON(config_file)
# Extract paths from JSON
model_dir <- config$model_folder               # Directory containing model files and min-max files
input_dir <- file.path(config$out_folder, "CRE_input") # Directory for input gene files (txt.gz) 
min_max_dir <- file.path(config$provided_input, "MinMax_CRE_RF")  # Directory for min-max files 

# Define the path for the result subfolder
result_dir <- file.path(config$out_folder, "CRE_RF_result")

# Create the result subfolder if it doesn't already exist
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

min_max_normalize_test <- function(data, min_values, max_values) {
  normalized_data <- data.frame(
    lapply(names(data), function(col) {
      min_val <- min_values[col]
      max_val <- max_values[col]
      print(min_val)
      print(max_val)
      if (max_val == min_val) {
        rep(0, length(data[[col]]))
      } else {
        (data[[col]] - min_val) / (max_val - min_val)
      }
    })
  )
  colnames(normalized_data) <- colnames(data)
  return(normalized_data)
}

# Load min-max values from a single file per gene
load_min_max <- function(gene_name) {
  min_max_file <- file.path(min_max_dir, paste0(gene_name, "_min_max.txt.gz"))
  
  # Read min-max values for each column
  min_max_data <- read.table(gzfile(min_max_file), header = TRUE, sep = "\t")  
  # Ensure Min and Max columns are numeric
  min_values <- as.numeric(min_max_data$Min)
  max_values <- as.numeric(min_max_data$Max)
  # Handle non-numeric columns 
  if (any(is.na(min_values)) | any(is.na(max_values))) {
    stop("Error: Min or Max columns contain non-numeric values or NAs.")
  }
  #return(list(min = min_values, max = max_values))
  return(list(min = min_values, max = max_values, data = min_max_data))  # Include full data
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
  # Convert dashes (-) in column names to periods (.)
  colnames(data) <- gsub("-", ".", colnames(data))  
  # Storing the 'Sample' column separately before removal
  sample_column <- data$Sample
  # Removing 'Sample' column if it exists
  if ("Sample" %in% colnames(data)) {
    data <- data[, !colnames(data) %in% "Sample"]  # Exclude 'Sample' column
  }
  
  min_max <- load_min_max(gene_name)
  min_max_data <- as.data.frame(min_max$data)  # Extract the full min-max data frame
  train_min_values <- setNames(min_max_data$Min, min_max_data$Feature)
  train_max_values <- setNames(min_max_data$Max, min_max_data$Feature)

  # Apply log transformation: log2(data + 1)
  log_transformed_data <- log2(data + 1)

  colnames(log_transformed_data) <- gsub("X","", colnames(log_transformed_data))
  # min-max normalization
  scaled_data <- min_max_normalize_test(data = log_transformed_data, min_values = train_min_values, max_values = train_max_values) #remove last element which is expr

  colnames(scaled_data) <- paste0("X", colnames(scaled_data))
  # Run predictions using the loaded model
  predictions <- predict(model, newdata = scaled_data)

  min_exp <- as.numeric((min_max_data$Min[min_max_data$Feature == "Expression"]))
  max_exp <- as.numeric((min_max_data$Max[min_max_data$Feature == "Expression"]))
  back_scaled_predictions <- (predictions * (max_exp - min_exp)) + min_exp
  
  # Reverse log transformation
  original_predictions <- 2^back_scaled_predictions - 1
  
  # Create final output with both columns
  final_output <- data.frame(Sample = sample_column, 
                             Prediction = predictions, 
                             Original_space_Prediction = original_predictions)
  
  # Save predictions to the result directory with the same gene name, as .csv
  output_file <- file.path(result_dir, paste0(gene_name, "_predictions.csv"))
  write.csv(final_output, output_file, row.names = FALSE)
  
  message(paste("Predictions saved for gene:", gene_name, "in", output_file))
}

message("Prediction process completed.")
