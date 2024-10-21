# Load necessary libraries
library(randomForest)
library(dplyr)
library(utils)

# Define your file paths
zenodo_url <- "https://zenodo.org/record/your_record_id/files/models.zip"  # Replace with the actual Zenodo URL (its it README file)
download_path <- "/path/to/download/directory"
unzip_path <- "/path/to/unzip/directory"

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
predict_gene_expression <- function(gene_name, train_samples_file, test_samples_file, input_data_dir, models_dir, output_dir) {
  
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
    
    # Split into train and test samples
    train_samples <- read.csv(train_samples_file)[, 1]
    test_samples <- read.csv(test_samples_file)[, 1]
    train_data <- df3[df3$df_samples %in% train_samples, ]
    test_data <- df3[df3$df_samples %in% test_samples, ]
    
    # Normalize train data
    train_min_values <- apply(train_data, 2, min)
    train_max_values <- apply(train_data, 2, max)
    train_data_normalized <- as.data.frame(apply(train_data[, 2:ncol(train_data)], 2, normalize_min_max))
    
    # Normalize test data using train min/max values
    test_data_normalized <- min_max_normalize_test(test_data[, 2:ncol(test_data)], train_min_values, train_max_values)
    
    # Train data target
    training <- train_data_normalized[, 1:ncol(train_data_normalized)-1]
    training_target <- train_data_normalized[, ncol(train_data_normalized)]
    
    # Test data target
    test1 <- test_data_normalized[, 1:ncol(test_data_normalized)-1]
    test1_target <- test_data_normalized[, ncol(test_data_normalized)]
    
    # Load the pre-trained model from the local directory
    model_path <- file.path(models_dir, paste0(gene_name, ".RDS"))
    if (!file.exists(model_path)) {
      stop(paste("Model for gene", gene_name, "not found"))
    }
    Best_model <- readRDS(model_path)
    
    # Predict on test data
    predicted_test <- predict(Best_model, test1)
    
    # Calculate MSE and correlation
    MSE1 <- mean((test1_target - predicted_test)^2)
    cor1 <- cor(test1_target, predicted_test)
    
    # Rescale the predictions back to the original scale
    min_origin_target <- min(training_target)
    max_origin_target <- max(training_target)
    scaled_predicted_target <- rescale(predicted_test, min_origin_target, max_origin_target)
    scaled_original_target <- rescale(test1_target, min_origin_target, max_origin_target)
    MSE_backscaled <- mean((scaled_original_target - scaled_predicted_target)^2)
    
    # Save the results
    results <- data.frame(
      Gene = gene_name,
      MSE_RF = MSE1,
      Cor_RF = cor1,
      MSE_back_scaled_RF = MSE_backscaled
    )
    write.csv(results, file.path(output_dir, paste0(gene_name, "_prediction_results.csv")), row.names = FALSE)
    
    return(results)
  } else {
    stop(paste("Gene data for", gene_name, "not found"))
  }
}

# Main function to process multiple genes
predict_multiple_genes <- function(gene_list_file, train_samples_file, test_samples_file, input_data_dir, models_dir, output_dir) {
  gene_list <- read.table(gene_list_file, header = FALSE, col.names = "Gene")
  
  results <- lapply(gene_list$Gene, function(gene_name) {
    predict_gene_expression(gene_name, train_samples_file, test_samples_file, input_data_dir, models_dir, output_dir)
  })
  
  final_results <- do.call(rbind, results)
  write.csv(final_results, file.path(output_dir, "all_genes_prediction_results.csv"), row.names = FALSE)
}
