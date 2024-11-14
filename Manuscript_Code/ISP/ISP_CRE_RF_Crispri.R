#ISP for RF (setting all features to 0 and predicting the expression, ISP is log ratio of old to new predicted value)
#loading libraries
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
library(caret)
library(randomForest)
library(doParallel)
library(foreach)
no_model <- c()
#reading the data for each gene
result_dataframe = data.frame()
file_path = ("/path/to/CRISPRi/genes/genes_CRISPR.txt")
file_contents <- readLines(file_path)
file_contents <- file_contents
print("_____________________________________")
# Set the number of cores to use
num_cores <- 8  # Change this to the desired number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

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

print("start of the function")
process_gene <- function(p) {
  
  set.seed(200)
  tryCatch({
    gene_name <- p
    print(paste("Processing gene:", gene_name))
    if(file.exists(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"))){
      print("Data Exist")
      df = read.table(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
      if(file.exists(paste0("/projects/apog/work/models/1MB/RF/",gene_name, ".RDS"))){
        print("Model Exist")
        
        my_model <- readRDS(paste0("/projects/apog/work/models/1MB/RF/RF_last/",gene_name, ".RDS"))
        df <- df[order(df$Sample),] #ordering the sample names
        rownames(df) = NULL
        df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
        df_samples = df[,1] 
        df2 = df[, 2:length(df)] 
        df2 = log2(df2 + 1) #whole data after log2, 
        df3 = cbind(df_samples, df2) #df3 is my dataset after log2 + 1
        train_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_train.csv"
        
        
        train_samples2 <- read.csv(train_samples)
        train_samples3 <- train_samples2[,1]
        train_data <- df3[df3$df_samples %in% train_samples3, ]
        
        #rename the train and test data
        train_data = train_data[, 2:length(train_data)] #whole
        train_traget_just_log_no_norm = train_data[,ncol(train_data)] #target original
        train_data2 = train_data[,1: (ncol(train_data2)) - 1] #whole without target
        
        train_data_normalized <- as.data.frame(apply(train_data2, 2, normalize_min_max))
        rownames(train_data_normalized) <- NULL
        rownames(train_data_normalized) <- sample_final
        
        geneA <- train_data_normalized[rownames(train_data_normalized) == "IHECRE00001887",]
        
        
        
        my_matrix <- matrix(geneA, nrow = 1, ncol = length(geneA))
        colnames(my_matrix) <- names(geneA)
        geneA <- my_matrix
    
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
          
          return(result_row)
        })
        
        # Combine the list of result rows into a data frame
        result_df <- do.call(rbind, result_list)
        
        write.table(result_df, file = paste0("/projects/apog/work/IHEC/ValidateInteractions/ENCODE_RF_ism_all_reverse/IHECRE00001887/",gene_name,"_RF_all_reverse.txt"),
                    sep = "\t", row.names = FALSE,  quote = FALSE)
        
        
        
        
      }else{
        print("the model dosnt exist")
        no_model[gene_name] <- gene_name
        
      }
    } else{
      print("NO data")
      
    }
    
    
    
    
  },  error=function(e){
    print(paste("Error occurred at gene", gene_name)); message(e)
    
    
  })
  
}
# Parallelize the loop over genes
foreach(p = file_contents, .packages=c("dplyr", "magrittr", "plyr", "boot", "tfruns", "MLmetrics", "parallel", "snow", "reticulate", "randomForest", "keras", "doParallel")) %dopar% {
  process_gene(p)
}

# Stop the parallel cluster
stopCluster(cl)

print("done")