
#RandomForest method

#required packages
library(plyr)
library(boot)
library(tfruns)
library(MLmetrics)
library(parallel)
library(snow)
library(tidyverse)
library(reticulate)
library(doParallel)
library(foreach)
library(randomForest)
library(dplyr)

set.seed(0)
# Paralizatio: Set the number of cores
num_cores <- 40 # assign number of available cores here
cl <- makeCluster(num_cores)
registerDoParallel(cl)

df = NULL
acc_nested_list = list()

file_path = ("path/to/kept_genes.txt") #set path to list of gens 
file_contents <- read.table(file_path, header = FALSE, col.names = c("Index", "Gene"))$Gene
file_contents = file_contents


results <- foreach(p = file_contents, .combine = bind_rows) %dopar% {
  set.seed(0)
  tryCatch({
    library(dplyr)  
    #loading related test data
    gene_name <- p
    
    #loading input file
    if(file.exists(paste0("/path/to/input/files/", gene_name, ".txt.gz"))){
      print("TRUE")
      df = read.table(paste0("/path/to/input/files/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
    }
    
    df <- df[order(df$Sample),] #ordering the sample names
    rownames(df) = NULL
    df = df[, colSums(df != 0) > 0] #if all values of a column is 0, remove that column
    df_samples = df[,1] #GET THE SAMPLE name, first col
    df2 = df[, 2:length(df)] 
    df2 = log2(df2 + 1) #whole data after log2, 
    df3 = cbind(df_samples, df2) 
    
    test_samples <- "/path/to/partition0_test.csv"
    train_samples <- "path/to/partition0_train.csv"
    
    test_sample2 <- read.csv(test_samples)
    test_sample3 <- test_sample2[,1]
    test_data <- df3[df3$df_samples %in% test_sample3, ]
    
    train_samples2 <- read.csv(train_samples)
    train_samples3 <- train_samples2[,1]
    train_data <- df3[df3$df_samples %in% train_samples3, ]
    
    #rename the train and test data
    train_data = train_data[, 2:length(train_data)]
    test_data = test_data[, 2:length(test_data)]
    train_traget_just_log_no_norm = train_data[,ncol(train_data)]
    
    
    train_min_values <- apply(train_data, 2, min)
    train_max_values <- apply(train_data, 2, max)
    #Normalization function
    normalize_min_max <- function(data) {
      if(max(data) == min(data)){
        return(0)
      }
      
      return ((data - min(data)) / (max(data) - min(data)))
    }
    
    train_data_normalized <- as.data.frame(apply(train_data, 2, normalize_min_max))
    
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
    
    # Applying the normalization function to the test data
    normalized_test <- min_max_normalize_test(data = test_data, min_values = train_min_values, max_values = train_max_values)
    
    test_data_normalized = normalized_test
    
    #after normalization(split the data and targets)
    training = train_data_normalized[,1:ncol(train_data_normalized)-1]
    training_target = train_data_normalized[,ncol(train_data_normalized)]
    test1 = test_data_normalized[,1:ncol(test_data_normalized)-1]
    test1_target = test_data_normalized[,ncol(test_data_normalized)]
    
    #rescaling
    rescale <- function(x, min_val, max_val) {
      return((x * (max_val - min_val)) + min_val)
    }
    
    temp = NULL
    a2 = floor(sqrt(ncol(training)))
    mtry1 = a2
    ntree1 = 501 

    Best_model <- randomForest::randomForest(x = training, y = training_target, ntree = ntree1, mtry = a2, importance = TRUE)
    
    name1 = paste0("path/to/directory/for/saving/models/",gene_name,".RDS") #saving the best model
    saveRDS(Best_model, name1)
    
    predicted_test = Best_model%>% predict(test1)
    predicted_test = as.numeric(predicted_test)
    
    MSE1 = mean((test1_target - predicted_test)^2)
    cor1 = cor(test1_target, predicted_test)
    
    min_origin_target =  min(train_traget_just_log_no_norm)
    max_origin_target = max(train_traget_just_log_no_norm)
    
    scaled_predicted_target <- rescale(predicted_test, min_val = min_origin_target, max_val = max_origin_target)
    scaled_original_target <- rescale(test1_target, min_val = min_origin_target, max_val = max_origin_target)
    
    MSE_backscaled = mean((scaled_original_target - scaled_predicted_target)^2)

    accuracy1 = data.frame(
      MSE_RF = MSE1,
      Cor_RF = cor1, 
      MSEs_back_scaled_RF = MSE_backscaled,
      Name = gene_name
    )
    
    
    # Run for-loop over lists
    acc_nested_list[[gene_name]] <- (accuracy1)
    
  }, error = function(e) {
    print(paste("Error occurred at gene", gene_name)); message(e)
    return(NULL)
  })
}

# Combine results from all cores
final_results <- bind_rows(results)

# Save the final results
name2 <- "path/to/directory/for/saving/statistics/RF_1MB_last.RDS"
saveRDS(final_results, name2)

# Stop the parallel cluster
stopCluster(cl)

print("end of the code")

