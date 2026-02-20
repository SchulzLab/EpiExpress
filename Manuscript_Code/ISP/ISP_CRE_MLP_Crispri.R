print("last run MLP all 0 scaled")
# Load the GenomicRanges library
library(GenomicRanges)
library(dplyr)
library(magrittr)
library(plyr) 
library(boot)
library(tfruns)
print("after tfruns")
library(MLmetrics)
print("after MLmetrics")
library(parallel)
library(snow)
print("after snow")
#library(tidyverse)
library(reticulate)
library(tensorflow)
library(keras)
library(doParallel)
library(foreach)


#reverse normalization(in this code we do it just for train because K562 is in train so we dont need to consider it for test)
inverse_normalize_min_max <- function(data1, original_data) {
  return (data1 * (max(original_data) - min(original_data)) + min(original_data))
}

#reverse log2 and +1
reverse_log2 <- function(data) {
  original_data <- (2^data) - 1
  return(original_data)
}

#Normalization function
normalize_min_max <- function(data) {
  if(max(data) == min(data)){
    return(0)
  }
  
  return ((data - min(data)) / (max(data) - min(data)))
}


# Iterate over gene names in file_contents


#reading the data for each gene
result_dataframe = data.frame()
file_path = ("/projects/apog/work/models/1MB/genes_CRISPR.txt")
file_contents <- readLines(file_path)
print(length(file_contents))
print("_____________________________________")

MLP2 <- readRDS("/path/to/data//Best_MLP_seeds.RDS")
# Set the number of cores to use
num_cores <- 48 # Change this to the desired number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Define a function for parallel processing
process_gene <- function(p) {
  
  set.seed(200)
  tryCatch({
    gene_name <- p
    print(paste("Processing gene:", gene_name))
    
    
    if (gene_name %in% MLP2$Name) {
      
      gene_row <- MLP2[MLP2$Name == gene_name, ]
      MLP <- gene_row$MLP
      
      # Get Pearson, MSE, and scaled MSE values
      pearson_r <- gene_row$Cor_test
      MSE <- gene_row$MSE_test
      scaled_MSE <- gene_row$MSE_S_test
      
      
      if (MLP == "NN1") {
        
        if(file.exists(paste0("/path/to/data//", gene_name, ".txt.gz"))){
          print("Data Exist")
          df = read.table(paste0("/path/to/data//", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
          
          if(file.exists(paste0("/path/to/data//",gene_name, ".hdf5")))
          {
            print("Model Exist")
            
            my_model <- load_model_hdf5(paste0("/path/to/data//",gene_name, ".hdf5"), custom_objects = NULL, compile = TRUE)
            
            df <- df[order(df$Sample),] #ordering the sample names
            rownames(df) = NULL
            df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
            df_samples = df[,1] #GET THE SAMPLE name, first col
            df2 = df[, 2:length(df)] 
            df2 = log2(df2 + 1) #whole data after log2, 
            df3 = cbind(df_samples, df2) #df3 is my dataset after log2 + 1
            
            train_samples <- "/path/to/data//partition0_train.csv"
            
            
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
            
      
            geneA <- as.matrix(geneA)
            
            
            
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
            # print(a)
            
            
            result_list <- list()
            result_list <- lapply(column_names, function(col_name) {
              # print(col_name)
              
              column_parts <- strsplit(col_name, "\\.")[[1]]
              chr_numeric <- gsub("X", "", column_parts[1])
              if (is.na(chr_numeric) || chr_numeric == "") {
                chr_numeric <- "X"
              }else {
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
            
            # write.table(result_df, file = paste0("/path/to/data//",gene_name,"_MLP_seed_checked.txt"),
            #             sep = "\t", row.names = FALSE,  quote = FALSE)
            
            
            
            header_line <- paste0("# pearson_r = ", pearson_r, " .  MSE = ", MSE, " .  scaled_MSE = ", scaled_MSE)
            output_file <- paste0("/path/to/data//", gene_name, "_MLP_all_seeds.txt")
            write(header_line, file = output_file)
            write.table(result_df, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
            
            
            
            
          }else{
            print("the model dosnt exist")
            # no_model[gene_name] <- gene_name
            
          }
        } else{
          print("NO data")
          
        }
        #NN2
      } else{
        if(file.exists(paste0("/path/to/data//", gene_name, ".txt.gz"))){
          print("Data Exist")
          df = read.table(paste0("/path/to/data//", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
          
          if(file.exists(paste0("/path/to/data//",gene_name, ".hdf5")))
            
          {
            print("Model Exist")
            
            my_model <- load_model_hdf5(paste0("/path/to/data//",gene_name, ".hdf5"), custom_objects = NULL, compile = TRUE)
            
            df <- df[order(df$Sample),] #ordering the sample names
            rownames(df) = NULL
            df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
            df_samples = df[,1] #GET THE SAMPLE name, first col
            df2 = df[, 2:length(df)] 
            df2 = log2(df2 + 1) #whole data after log2, 
            df3 = cbind(df_samples, df2) #df3 is my dataset after log2 + 1
            
            train_samples <- "/path/to/data//partition0_train.csv"
            
            
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
            
            
            
            # my_matrix <- matrix(geneA, nrow = 1, ncol = length(geneA))
            # colnames(my_matrix) <- names(geneA)
            # geneA <- my_matrix
            # removing last column which is expression
            geneA <- as.matrix(geneA)
            
            
            
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
              if (is.na(chr_numeric) || chr_numeric == "") {
                chr_numeric <- "X"
              }else {
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
            
            # write.table(result_df, file = paste0("/path/to/data//",gene_name,"_MLP_seed_checked.txt"),
            #             sep = "\t", row.names = FALSE,  quote = FALSE)
            # 
            
            
            header_line <- paste0("# pearson_r = ", pearson_r, " .  MSE = ", MSE, " .  scaled_MSE = ", scaled_MSE)
            output_file <- paste0("/path/to/data//", gene_name, "_MLP_all_seeds.txt")
            write(header_line, file = output_file)
            write.table(result_df, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
            
            
            
          }else{
            print("the model dosnt exist")
            # no_model[gene_name] <- gene_name
            
          }
        } else{
          print("NO data")
          
        }
      }
      
    } else {
      # If gene_name is not found in MLP2, consider source as NN2
      print("we dont have this gene in our statistics")
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

print("end one")
