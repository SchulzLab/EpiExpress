print("hi RF all to 0 CRISPERI proper reverse version for luekemia")

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
#sig genes for RF
file_path = ("/projects/apog/work/models/1MB/results/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
file_contents <- readLines(file_path)
print(length(file_contents))
print("_____________________________________")
# Set the number of cores to use
num_cores <- 16  # Change this to the desired number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#train normalization
normalize_min_max <- function(data) {
  if(max(data) == min(data)){
    return(0)
  }
  
  return ((data - min(data)) / (max(data) - min(data)))
}


#test normalization
min_max_normalize_test <- function(data, min_values, max_values) {
  normalized_data <- data.frame(
    lapply(names(data), function(col) {
      min_val <- min_values[col]
      max_val <- max_values[col]
      
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


#back scaling
inverse_normalize_min_max <- function(data1, original_data) {
  return (data1 * (max(original_data) - min(original_data)) + min(original_data))
}

#reverse log2 and +1
reverse_log2 <- function(data) {
  original_data <- (2^data) - 1
  return(original_data)
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
      if(file.exists(paste0("/projects/apog/work/models/1MB/RF/RF_last/",gene_name, ".RDS"))){
        print("Model Exist")
        my_model <- readRDS(paste0("/projects/apog/work/models/1MB/RF/RF_last/",gene_name, ".RDS"))
        
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
        train_min_values <- apply(train_data, 2, min)
        train_max_values <- apply(train_data, 2, max)
        train_traget_just_log_no_norm = train_data[,ncol(train_data)] #target original
        min_train_target_for_rescaling = min(train_traget_just_log_no_norm)
        max_train_target_for_rescaling = max(train_traget_just_log_no_norm)
        train_data2 = train_data[,1: (ncol(train_data)) - 1] #whole without target
        
        train_data_normalized <- as.data.frame(apply(train_data2, 2, normalize_min_max))
        rownames(train_data_normalized) <- NULL
        rownames(train_data_normalized) <- sample_final
        
        #test data: here is leukemia input 
        df_test = read.table(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38_leukemia/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
        rownames(df_test) = NULL
        df_test = df_test[, colSums(df_test != 0) > 0] #if all values of a column is 0 Ill remove that column
        df_samples = df_test[,1] #GET THE SAMPLE name, first col
        df2_test = df_test[, 2:length(df_test)] 
        df2 = log2(df2_test + 1) #whole data after log2, 
        df3_test = cbind(df_samples, df2_test) #df3 is my dataset after log2 + 1
        # new version which outleires are removed (we pick desired sample, helthy and disease)
        luekemia_samples <- readLines("/projects/apog/work/models/1MB/results/samples_leukemia_adult_removed_outliers.txt")
        print(length(luekemia_samples))
        test_data <- df3_test[df3_test$df_samples %in% luekemia_samples, ]
        dim(test_data)
        
        #rename the train and test data
        test_data2 = test_data[, 2:length(test_data)] #whole
        test_data_normalized <- min_max_normalize_test(data = test_data2, min_values = train_min_values, max_values = train_max_values)
        rownames(test_data_normalized) <- NULL
        rownames(test_data_normalized) <- luekemia_samples
        
        
        geneA <- test_data_normalized
        geneA <- as.matrix(geneA)
        dim(geneA)
        result_list <- list()
        
        result_S <- matrix(nrow = nrow(geneA), ncol = ncol(geneA))
        colnames(result_S) <- colnames(geneA)
        rownames(result_S) <- rownames(geneA)
        
        column_names <- colnames(geneA)
        original_geneA1 <- geneA
        original_target <- train_traget_just_log_no_norm #train target column
        #prediction of original data
        a <- predict(my_model, original_geneA1)
        print(a)
        print(class(original_target))
        print(length(original_target))
        #scaling part
        a <- a * (max_train_target_for_rescaling - min_train_target_for_rescaling) + min_train_target_for_rescaling
        a <- (2^a) - 1
        
      
        
        row_counter <- 1
        
        result_list <- list()
        # Iterate over each column in geneA
        for (col_index in seq_along(column_names)) {
          col_name <- column_names[col_index]
          column_parts <- strsplit(col_name, "\\.")[[1]]
          chr_numeric <- gsub("X", "", column_parts[1])
          
          # Make a copy of geneA to modify
          modified_geneA <- geneA
          
          # Set the column(feature) to zero for the feature being evaluated
          modified_geneA[, col_name] <- 0
          
          # Predict with the modified data
          b <- predict(my_model, modified_geneA)
          b <- b * (max_train_target_for_rescaling - min_train_target_for_rescaling) + min_train_target_for_rescaling
          b <- (2^b) - 1  # Reverse Log2 Transformation for ISM
          
          # Calculate S for each sample
          S <- (log2((b + 1) / (a + 1)))
          S_txt <- S
          S_df <- -1 * (S)
          # Place the S values in the correct row of result_S
          result_S[, col_index] <- S_df
          
          # Create a data frame for the current column and add it to the result_list
          result_row <- data.frame(
            chr = chr_numeric,
            start = as.numeric(column_parts[2]),
            end = as.numeric(column_parts[3]),
            EnsemblID = gene_name,
            predicted = a,
            ISM = b,
            S = S_txt,
            stringsAsFactors = FALSE
          )
          
          result_list[[col_index]] <- result_row
        }
        
        # Combine the list of result rows into a data frame
        result_df <- do.call(rbind, result_list)
        
        
        write.table(result_df, file = paste0("provide/a/path/for/saving/the/result/",gene_name,"_RF_all_reverse.txt"),
                    sep = "\t", row.names = FALSE,  quote = FALSE)
        saveRDS(result_S, file = paste0("provide/a/path/for/saving/the/result/",gene_name,".RDS"))
        
        
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

foreach(p = file_contents, .packages = c("dplyr", "magrittr", "plyr", "boot", "tfruns", "MLmetrics", "parallel", "snow", "reticulate", "randomForest", "doParallel", "foreach")) %dopar% {
  process_gene(p)
}


# Stop the parallel cluster
stopCluster(cl)

print("end of code")