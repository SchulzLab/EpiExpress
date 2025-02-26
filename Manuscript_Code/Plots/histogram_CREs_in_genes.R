library(parallel)

# Define file path and get list of files
file_path <- "/projects/apog/work/input/IHEC_Activity_1MB_hg38/"
files <- list.files(file_path, full.names = TRUE)

# Set number of cores (leave one free for stability)
num_cores <- min(42, detectCores() - 1)

# Function to process each file
process_files <- function(file) {
  gene_name <- sub("\\.txt\\.gz$", "", basename(file))
  
  # Try to read file, return NA if fails
  gene <- tryCatch({
    read.table(gzfile(file), header = TRUE)
  }, error = function(e) return(NA))
  
  # Check if data is valid before processing
  if (!is.data.frame(gene)) return(setNames(NA, gene_name))
  
  # Count columns minus one
  tmp <- length(colnames(gene)) - 1  
  return(setNames(tmp, gene_name)) 
}

# Apply parallel processing
results <- mclapply(files, process_files, mc.cores = num_cores)

# Convert list to named numeric vector, removing NAs
col_vec <- unlist(results)
col_vec <- col_vec[!is.na(col_vec)]  # Remove NA values

# Read the kept genes list (assuming gene names are in the second column)
kept_genes <- read.table("/projects/apog/work/CNN_bin/miscellaneous/kept_genes.txt",
                         header = FALSE, skip = 1)[, 2]

# Filter valid genes
kept_genes <- kept_genes[kept_genes %in% names(col_vec)]
col_vec_kept <- col_vec[kept_genes]

# Save results
write.table(col_vec_kept, "/projects/apog/work/models/1MB/results/number_of_CREs_kept_genes.txt",
            col.names = FALSE, quote = FALSE)
write.table(col_vec, "/projects/apog/work/models/1MB/results/number_of_CREs_All_genes.txt",
            col.names = FALSE, quote = FALSE)

# Generate histogram only if col_vec is non-empty
if (length(col_vec) > 0) {
  pdf(file = "/projects/apog/work/models/1MB/results/hist_CREs_kep_genes.pdf", width = 5, height = 5)
  hist(col_vec_kept, breaks = 30, main = "Histogram of Number of CREs per Gene",
       xlab = "Number of CREs", col = "blue", border = "black")
  dev.off()
} else {
  message("⚠️ Warning: No valid data for histogram!")
}
