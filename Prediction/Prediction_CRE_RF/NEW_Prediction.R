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


# new function to check column names
scale_data <- function(data, min_max_data) {
  # Extract min-max values based on matching column names
  min_values <- setNames(min_max_data$Min, min_max_data$Feature)
  max_values <- setNames(min_max_data$Max, min_max_data$Feature)
  
  # Find common columns
  common_cols <- intersect(colnames(data), names(min_values))
  
  # Warn if some expected columns are missing in data
  missing_cols <- setdiff(names(min_values), colnames(data))
  if (length(missing_cols) > 0) {
    warning("Missing columns in input data: ", paste(missing_cols, collapse = ", "))
  }

  # Warn if there are extra columns in data
  extra_cols <- setdiff(colnames(data), names(min_values))
  if (length(extra_cols) > 0) {
    warning("Extra columns in input data (will not be scaled): ", paste(extra_cols, collapse = ", "))
  }

  # Apply min-max normalization only on matching columns
  scaled_data <- data
  for (col in common_cols) {
    min_val <- min_values[col]
    max_val <- max_values[col]
    range_val <- max_val - min_val

    if (range_val != 0) {
      scaled_data[[col]] <- (data[[col]] - min_val) / range_val
    } else {
      scaled_data[[col]] <- 0  # Avoid division by zero
    }
  }
  
  return(scaled_data)
}


# Load min-max values from a single file per gene
load_min_max <- function(gene_name) {
  min_max_file <- file.path(min_max_dir, paste0(gene_name, "_min_max.txt.gz"))
  
  # Read min-max values for each column
  min_max_data <- read.table(gzfile(min_max_file), header = TRUE, sep = "\t") #reading gz files 
 # min_max_data <- read.table(min_max_file, header = TRUE, sep = "\t")  
  # Ensure Min and Max columns are numeric
  min_values <- as.numeric(min_max_data$Min)
  max_values <- as.numeric(min_max_data$Max)
  # Handle non-numeric columns 
  if (any(is.na(min_values)) | any(is.na(max_values))) {
    stop("Error: Min or Max columns contain non-numeric values or NAs.")
  }
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
  # Convert dashes (-) in column names to periods (.)
  colnames(data) <- gsub("-", ".", colnames(data))  
  # Storing the 'Sample' column separately before removal
  sample_column <- data$Sample
  # Removing 'Sample' column if it exists
  if ("Sample" %in% colnames(data)) {
    data <- data[, !colnames(data) %in% "Sample"]  # Exclude 'Sample' column
  }
  
  min_max <- load_min_max(gene_name)
  
  # Apply log transformation: log2(data + 1)
  log_transformed_data <- log2(data + 1)
  # data scaling
  scaled_data <- scale_data(log_transformed_data, min_max$min, min_max$max)
  
  # Run predictions using the loaded model
  predictions <- predict(model, newdata = scaled_data)
  
  # Add the 'Sample' column back to the predictions data frame
  final_output <- cbind(Sample = sample_column, Prediction = predictions)
  
  # Save predictions to the result directory with the same gene name, as .csv
  output_file <- file.path(result_dir, paste0(gene_name, "_predictions.csv"))
  write.csv(final_output, output_file, row.names = FALSE)
  
  message(paste("Predictions saved for gene:", gene_name, "in", output_file))
}

message("Prediction process completed.")
