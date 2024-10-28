library(jsonlite)

# Load config
config <- fromJSON("path/to/config.json")

# Extract paths from config
model_path <- config$model_path
input_dir <- config$input_data
output_path <- config$output_path
min_max_file <- config$min_max_file

# Load min/max values
min_max_data <- read.csv(min_max_file, row.names = 1)

# Function to normalize and predict gene expression
predict_expression_rf <- function(gene_file) {
  gene_name <- gsub(".txt.gz$", "", basename(gene_file))
  model_file <- file.path(model_path, paste0(gene_name, ".RDS"))
  
  if (file.exists(model_file)) {
    rf_model <- readRDS(model_file)
    gene_data <- read.csv(gzfile(gene_file))
    
    # Apply min/max normalization if values are available for the gene
    if (gene_name %in% rownames(min_max_data)) {
      min_val <- min_max_data[gene_name, "min"]
      max_val <- min_max_data[gene_name, "max"]
      gene_data <- (gene_data - min_val) / (max_val - min_val)
    } else {
      message(paste("Min/max values for gene", gene_name, "not found. Skipping normalization."))
    }
    
    predictions <- predict(rf_model, gene_data)
    write.csv(predictions, file.path(output_path, paste0(gene_name, "_predictions.csv")))
  } else {
    message(paste("Model for gene", gene_name, "not found. Skipping."))
  }
}

# Run predictions for all gene files in the input directory
gene_files <- list.files(input_dir, pattern = "\\.txt.gz$", full.names = TRUE)
lapply(gene_files, predict_expression_rf)

