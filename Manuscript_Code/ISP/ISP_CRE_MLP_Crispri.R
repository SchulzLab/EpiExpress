print("ISP for best MLP(NN1 or NN2), all features set to 0, backscaled version")

library(GenomicRanges)
library(dplyr)
library(magrittr)
library(plyr) 
library(boot)
library(tfruns)
library(MLmetrics)
library(parallel)
library(snow)
library(reticulate)
library(tensorflow)
library(keras)
library(doParallel)
library(foreach)


#reading the data for each gene
result_dataframe = data.frame()
file_path = ("/projects/apog/work/models/1MB/genes_CRISPR.txt")
file_contents <- readLines(file_path)
file_contents <- file_contents
MLP2 <- readRDS("/projects/apog/work/models/1MB/merged_data.RDS")

num_cores <- 36  # Change this to the desired number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Define a function for parallel processing
process_gene <- function(p) {
  
  set.seed(200)
  tryCatch({
    gene_name <- p
    print(paste("Processing gene:", gene_name))
    
    
    
    if (gene_name %in% MLP2$Name) {
      # Find the corresponding row in MLP2
      gene_row <- MLP2[MLP2$Name == gene_name, ]
      
      # Check the Source column value
      source <- gene_row$MLP
    } else {
      # If gene_name is not found in MLP2, consider source as NN2
      source <- "NN2"
    }
    
    # Load the file based on the source
    if (source == "NN1") {
      
      if(file.exists(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"))){
        print("Data Exist")
        df = read.table(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
        
        if(file.exists(paste0("/projects/apog/work/models/1MB/NN1/",gene_name, ".hdf5")))
        {
          print("Model Exist")
          
          my_model <- load_model_hdf5(paste0("/projects/apog/work/models/1MB/NN1/",gene_name, ".hdf5"), custom_objects = NULL, compile = TRUE)
          
          df <- df[order(df$Sample),] #ordering the sample names
          rownames(df) = NULL
          df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
          df_samples = df[,1] #GET THE SAMPLE name, first col
          df2 = df[, 2:length(df)] 
          df2 = log2(df2 + 1) #whole data after log2, 
          df3 = cbind(df_samples, df2) #df3 is my dataset after log2 + 1
          
          train_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_train.csv"
          
          
          train_samples2 <- read.csv(train_samples)
          train_samples3 <- train_samples2[,1]
          train_data <- df3[df3$df_samples %in% train_samples3, ]
          sample_final <- train_data[,1]
          
          #rename the train and test data
          train_data = train_data[, 2:length(train_data)] #whole
          train_traget_just_log_no_norm = train_data[,ncol(train_data)] #target original
          train_data2 = train_data[,1: (ncol(train_data)) - 1] #whole without target
          
          train_data_normalized <- as.data.frame(apply(train_data2, 2, normalize_min_max))
          rownames(train_data_normalized) <- NULL
          rownames(train_data_normalized) <- sample_final
          
          geneA <- train_data_normalized[rownames(train_data_normalized) == "IHECRE00001887",]
          
          
          
          my_matrix <- matrix(geneA, nrow = 1, ncol = length(geneA))
          colnames(my_matrix) <- names(geneA)
          geneA <- my_matrix
          # removing last column which is expression
          # geneA <- as.matrix(geneA)
          
          
          
          column_names <- colnames(geneA)
          original_geneA1 <- geneA
          original_target <- train_traget_just_log_no_norm #train target column
          
          
          
          a <- predict(my_model, original_geneA1)
          print(a)
          print(class(original_target))
          print(length(original_target))
          # a <- inverse_normalize_min_max(a, original_target)
          a <- a * (max(original_target) - min(original_target)) + min(original_target)
          a <- (2^a) - 1
          # a <- reverse_log2(a)
          print(a)
          
          
          result_list <- list()
          result_list <- lapply(column_names, function(col_name) {
            # print(col_name)
            
            column_parts <- strsplit(col_name, "\\.")[[1]]
            chr_numeric <- gsub("X", "", column_parts[1])
            if (is.na(chr_numeric)) {
              chr_numeric <- "X"
            } else {
              print("its ok")
            }
            
            geneA[, col_name] <- 0
            
            b <- predict(my_model, geneA)
            
            # Reverse Normalization for ISM
            b <- inverse_normalize_min_max(b, original_target)
            
            # Reverse Log2 Transformation for ISM
            b <- reverse_log2(b)
            result_row <- data.frame(
              chr = chr_numeric,
              start = as.numeric(column_parts[2]),
              end = as.numeric(column_parts[3]),
              EnsemblID = gene_name,
              predicted = a,
              ISM = b,
              stringsAsFactors = FALSE
            )
            
            # print("made the dataframe for current column")
            return(result_row)
          })
          
          # Combine the list of result rows into a data frame
          result_df <- do.call(rbind, result_list)
          
          write.table(result_df, file = paste0("/projects/apog/work/IHEC/ValidateInteractions/ENCODE_MLP_ism_all/IHECRE00001887/",gene_name,"_NN1_all_reverse.txt"),
                      sep = "\t", row.names = FALSE,  quote = FALSE)
          
          
          
          
        }else{
          print("the model dosnt exist")
          
          
        }
      } else{
        print("NO data")
        
      }
      #NN2
    } else{
      if(file.exists(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"))){
        print("Data Exist")
        df = read.table(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
        
        if(file.exists(paste0("/projects/apog/work/models/1MB/NN2/",gene_name, ".hdf5")))
        {
          print("Model Exist")
          
          my_model <- load_model_hdf5(paste0("/projects/apog/work/models/1MB/NN2/",gene_name, ".hdf5"), custom_objects = NULL, compile = TRUE)
          
          df <- df[order(df$Sample),] #ordering the sample names
          rownames(df) = NULL
          df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
          df_samples = df[,1] #GET THE SAMPLE name, first col
          df2 = df[, 2:length(df)] 
          df2 = log2(df2 + 1) #whole data after log2, 
          df3 = cbind(df_samples, df2) #df3 is my dataset after log2 + 1
          
          train_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_train.csv"
          
          
          train_samples2 <- read.csv(train_samples)
          train_samples3 <- train_samples2[,1]
          train_data <- df3[df3$df_samples %in% train_samples3, ]
          sample_final <- train_data[,1]
          
          #rename the train and test data
          train_data = train_data[, 2:length(train_data)] #whole
          train_traget_just_log_no_norm = train_data[,ncol(train_data)] #target original
          train_data2 = train_data[,1: (ncol(train_data)) - 1] #whole without target
          
          train_data_normalized <- as.data.frame(apply(train_data2, 2, normalize_min_max))
          rownames(train_data_normalized) <- NULL
          rownames(train_data_normalized) <- sample_final
          
          geneA <- train_data_normalized[rownames(train_data_normalized) == "IHECRE00001887",]
          
          
          
          my_matrix <- matrix(geneA, nrow = 1, ncol = length(geneA))
          colnames(my_matrix) <- names(geneA)
          geneA <- my_matrix
          # removing last column which is expression
          
          column_names <- colnames(geneA)
          original_geneA1 <- geneA
          original_target <- train_traget_just_log_no_norm #train target column
          
          a <- predict(my_model, original_geneA1)
          print(a)
          print(class(original_target))
          print(length(original_target))
          # a <- inverse_normalize_min_max(a, original_target)
          a <- a * (max(original_target) - min(original_target)) + min(original_target)
          a <- (2^a) - 1
          # a <- reverse_log2(a)
          print(a)
          
          
          result_list <- list()
          result_list <- lapply(column_names, function(col_name) {
            # print(col_name)
            
            column_parts <- strsplit(col_name, "\\.")[[1]]
            chr_numeric <- gsub("X", "", column_parts[1])
            if (is.na(chr_numeric)) {
              chr_numeric <- "X"
            } else {
              print("its ok")
            }
            
            geneA[, col_name] <- 0
            
            b <- predict(my_model, geneA)
            
            # Reverse Normalization for ISM
            b <- inverse_normalize_min_max(b, original_target)
            
            # Reverse Log2 Transformation for ISM
            b <- reverse_log2(b)
            result_row <- data.frame(
              chr = chr_numeric,
              start = as.numeric(column_parts[2]),
              end = as.numeric(column_parts[3]),
              EnsemblID = gene_name,
              predicted = a,
              ISM = b,
              stringsAsFactors = FALSE
            )
            
            # print("made the dataframe for current column")
            return(result_row)
          })
          
          # Combine the list of result rows into a data frame
          result_df <- do.call(rbind, result_list)
          
          write.table(result_df, file = paste0("/projects/apog/work/IHEC/ValidateInteractions/ENCODE_MLP_ism_all/IHECRE00001887/",gene_name,"_NN2_all_reverse.txt"),
                      sep = "\t", row.names = FALSE,  quote = FALSE)
          
          
          
          
        }else{
          print("the model dosnt exist")
          # no_model[gene_name] <- gene_name
          
        }
      } else{
        print("NO data")
        
      }
    }
    
    
    
    
  },  error=function(e){
    print(paste("Error occurred at gene", gene_name)); message(e)
    
    
  })
  
}
# Parallelize the loop over genes
foreach(p = file_contents, .packages=c("GenomicRanges", "dplyr", "magrittr", "plyr", "boot", "tfruns", "MLmetrics", "parallel", "snow", "reticulate", "tensorflow", "keras", "doParallel")) %dopar% {
  process_gene(p)
}

# Stop the parallel cluster
stopCluster(cl)

