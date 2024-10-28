# Load necessary libraries
library(randomForest)
library(dplyr)
library(utils)

# Define your file paths
zenodo_url <- "https://zenodo.org/record/your_record_id/files/models.zip"  # Replace with the actual Zenodo URL
download_path <- "/path/to/download/directory"
unzip_path <- "/path/to/unzip/directory"
input_data_dir <- "/path/to/test_samples/directory"  # This now contains user input data
output_dir <- "/path/to/output/directory"

# Function to download the models zip file from Zenodo
download_models_from_zenodo <- function(zenodo_url, download_path, unzip_path) {
  # Download the zip file
  zip_file <- file.path(download_path, "models.zip")
  download.file(zenodo_url, zip_file, mode = "wb")
  
  # Unzip the file to the specified directory
  unzip(zip_file, exdir = unzip_path)
  message("Models downloaded and unzipped successfully.")
}

# Define normalization functions
normalize_min_max <- function(data) {
  if(max(data) == min(data)) {
    return(0)
  }
  return((data - min(data)) / (max(data) - min(data)))
}

min_max_normalize_test <- function(data, min_values, max_values) {
  normalized_data <- data.frame(
    lapply(names(data), function(col) {
      min_val <- min_values[col]
      max_val <- max_values[col]
      
      if (max_val == min_val) {
        rep(0, length(data[[col]]))
      } else {
        (data[[col]] - min_val) / (max_val - min(min_values[[col]], max_values[[col]]))
      }
    })
  )
  colnames(normalized_data) <- colnames(data)
  return(normalized_data)
}

# Rescaling function to reverse normalization
rescale <- function(x, min_val, max_val) {
  return((x * (max_val - min_val)) + min_val)
}

# Function to load model and predict gene expression
predict_gene_expression <- function(gene_name, input_data_dir, models_dir, output_dir) {
  # Load data for the specific gene
  gene_file_path <- file.path(input_data_dir, paste0(gene_name, ".txt.gz"))
  
  if (file.exists(gene_file_path)) {
    df <- read.table(gene_file_path, header = TRUE, sep = "\t")
    df <- df[order(df$Sample), ]
    rownames(df) <- NULL
    df <- df[, colSums(df != 0) > 0]  # Remove columns with all 0s
    df_samples <- df[, 1]
    df2 <- df[, 2:ncol(df)]
    df2 <- log2(df2 + 1)
    df3 <- cbind(df_samples, df2)
    
    # Prompt user for sample names
    sample_names_input <- readline(prompt = "Enter the sample names separated by commas: ")
    sample_names <- unlist(strsplit(sample_names_input, ","))
    
    # Select user-defined samples
    user_data <- df3[df3$df_samples %in% trimws(sample_names), ]
    
    # Normalize data
    train_min_values <- apply(user_data, 2, min)
    train_max_values <- apply(user_data, 2, max)
    user_data_normalized <- as.data.frame(apply(user_data[, 2:ncol(user_data)], 2, normalize_min_max))
    
    # Load the pre-trained model from the local directory
    model_path <- file.path(models_dir, paste0(gene_name, ".RDS"))
    if (!file.exists(model_path)) {
      stop(paste("Model for gene", gene_name, "not found"))
    }
    Best_model <- readRDS(model_path)
    
    # Predict on user data
    predicted_test <- predict(Best_model, user_data_normalized)
    
    # Save the results
    results <- data.frame(
      Gene = gene_name,
      Predicted_Values = predicted_test
    )
    write.csv(results, file.path(output_dir, paste0(gene_name, "_user_input_prediction_results.csv")), row.names = FALSE)
    
    return(results)
  } else {
    stop(paste("Gene data for", gene_name, "not found"))
  }
}

# Main function to process multiple genes
predict_multiple_genes <- function(gene_list_file, input_data_dir, models_dir, output_dir) {
  gene_list <- read.table(gene_list_file, header = FALSE, col.names = "Gene")
  
  results <- lapply(gene_list$Gene, function(gene_name) {
    predict_gene_expression(gene_name, input_data_dir, models_dir, output_dir)
  })
  
  final_results <- do.call(rbind, results)
  write.csv(final_results, file.path(output_dir, "all_genes_user_input_prediction_results.csv"), row.names = FALSE)
}

# Example of calling the main function
# predict_multiple_genes("gene_list.txt", input_data_dir, models_dir, output_dir)
