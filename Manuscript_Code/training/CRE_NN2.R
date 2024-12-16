#libraries
library(plyr) 
library(boot)
library(keras)
library(tfruns)
library(MLmetrics)
library(parallel)
library(snow)
library(reticulate)
library(tensorflow)
library(keras)
set.seed(123)



#V2 to 5601 on Bohr and tmux is: MLP2


file_path <- "/projects/apog/work/CNN_bin/miscellaneous/kept_genes.txt"

# Read the file, skipping the first row (header)
file_contents <- read.table(
  file_path,
  header = FALSE,       # Since the header is irregular, treat the file as headerless
  skip = 1,             # Skip the first line
  col.names = c("Index", "Gene"), # Define column names manually
  stringsAsFactors = FALSE # Keep strings as characters
)$Gene

print(length(file_contents))

# Generate reproducible seeds
set.seed(123)
seeds <- sample(1:10000, 30)
print(seeds)

df = NULL
acc_nested_list = list()

validation_split = 0.3

for(p in file_contents) { 
  
  tryCatch({
    
    gene_name <- p
    
    start_time <- Sys.time() # Record start time
    
    if(file.exists(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"))){
      print("TRUE")
      df = read.table(paste0("/projects/apog/work/input/IHEC_Activity_1MB_hg38/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
    } else{
      print("NOTTRUE")
    }
    
    
    
    print(gene_name)
    cl <- makeCluster(52)
    
    
    df <- df[order(df$Sample),] #ordering the sample names
    rownames(df) = NULL
    df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
    df_samples = df[,1] 
    df2 = df[, 2:length(df)] 
    df2 = log2(df2 + 1) #whole data after log2, 
    df3 = cbind(df_samples, df2)
    
    test_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_test.csv"
    train_samples <- "/projects/apog/work/CNN_bin/miscellaneous/partition0_train.csv"
    
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
    print(class(training))
    print(class(training_target))
    #rescaling
    rescale <- function(x, min_val, max_val) {
      return((x * (max_val - min_val)) + min_val)
    }
    training = as.matrix(training)
    temp = NULL
    
    dropouts = c(0.2, 0.4)
    batchsizes = 32 
    Epochs= 1200
    # learning_rates=c(0.001, 0.01)
    layers = c(1,2)
    dropouts2 = c(0.2, 0.4)
    

    best_loss = Inf
    best_loss2 = Inf
    
    
    #------------------------------------------------------
    
    # Main loop with parallelism for seeds
    print("Two hidden layers mode")
    for (j in 1:length(dropouts)) {
      for (jj in 1:length(dropouts2)) {
        # Export variables to parallel cluster
        clusterExport(cl = cl, list = c("dropouts", "dropouts2", "seeds", "j", "jj", "training", "training_target"))
        # print("Before second parLapply")
        
        # Parallel processing over seeds
        a2 <- parLapply(cl, seeds, function(seed) {
          require(keras)
          tensorflow::tf$random$set_seed(seed) # Set the seed for TensorFlow randomness
          
          # Build the model
          model <- keras_model_sequential() %>%
            layer_dense(100, activation = "relu", input_shape = c(dim(training)[2])) %>%
            layer_dropout(rate = dropouts[j]) %>%
            layer_dense(50, activation = "relu") %>%
            layer_dropout(rate = dropouts2[jj]) %>%
            layer_dense(units = 1, activation = "linear")
          
          # Compile the model
          model %>% compile(
            loss = 'mse',
            optimizer = optimizer_adam(),
            metrics = c('mse')
          )
          
          # Callbacks
          print_dot_callback <- callback_lambda(
            on_epoch_end = function(epoch, logs) {
              if (epoch %% 100 == 0) cat("\n")
              cat(".")
            })
          early_stop <- callback_early_stopping(monitor = "val_loss", mode = 'min', patience = 30)
          
          # Fit the model
          model1 <- model %>% fit(
            training,
            training_target,
            epochs = 1200,
            batch_size = 32,
            shuffled = FALSE,
            validation_split = 0.3,
            verbose = 0,
            callbacks = list(early_stop, print_dot_callback)
          )
          
          # Calculate the loss
          temp_loss <- mean(model1$metrics$val_loss)
          
          return(list(
            model = model1,
            loss = temp_loss,
            seed = seed # Include seed in the result
          ))
        })
        
        # Evaluate and save the best model
        for (w in 1:length(a2)) {
          if (a2[[w]]$loss < best_loss2) {
            best_loss2 <- a2[[w]]$loss
            b_hyper_dropouts <- dropouts[j]
            b_hyper_dropouts2 <- dropouts2[jj]
            best_seed <- a2[[w]]$seed
            best_model2 <- a2[[w]]$model
          }
        }
      } # End of dropouts2 loop
    } # End of dropouts loop
    

    
    validation_err = best_loss2
    
    #train again again
    tensorflow::tf$random$set_seed(best_seed)
    Best_model <- keras_model_sequential()
    Best_model %>%
      layer_dense(100, activation = "relu", input_shape = c(dim(training)[2])) %>%
      layer_dropout(rate = b_hyper_dropouts) %>%
      layer_dense(50, activation = "relu") %>%
      layer_dropout(rate = b_hyper_dropouts2) %>%
      layer_dense(units=1, activation ="linear")
    
    Best_model %>% compile(
      loss = 'mse',
      optimizer = optimizer_adam(), #default
      metrics = c('mse'))
    
    # Define callbacks (NEW)
    early_stop <- callback_early_stopping(monitor = "val_loss", mode = 'min', patience = 30)
    print_dot_callback <- callback_lambda(
      on_epoch_end = function(epoch, logs) {
        if (epoch %% 100 == 0) cat("\n")
        cat(".")
      })
    
    
    Best_ModelFited <- Best_model %>% fit(
      training, 
      training_target,
      epochs = 1200, 
      batch_size = 32,
      shuffled = F,
      validation_split = 0.3,
      verbose = 0,
      callbacks = list(early_stop, print_dot_callback) # Include callbacks (new)
    )
    

    
    # Save best model
    name1 = paste0("/projects/apog/work/models/1MB/new_MLP_seed/NN2/models/", gene_name, ".hdf5")
    save_model_hdf5(Best_model, name1)
    print("Best model is saved")
    
    min_origin_target =  min(train_traget_just_log_no_norm)
    max_origin_target = max(train_traget_just_log_no_norm)
    
    #test
    test1 = as.matrix(test1)
    predicted_test = Best_model%>% predict(test1)
    predicted_test = as.numeric(predicted_test)
    print("test_statistics")
    MSE1 = mean((test1_target - predicted_test)^2)
    cor1 = cor(test1_target, predicted_test, method = "pearson")
    scaled_predicted_target <- rescale(predicted_test, min_val = min_origin_target, max_val = max_origin_target)
    scaled_original_target <- rescale(test1_target, min_val = min_origin_target, max_val = max_origin_target)
    MSE_backscaled = mean((scaled_original_target - scaled_predicted_target)^2)
    
    
    #train
    predicted_train = Best_model%>% predict(training)
    predicted_train = as.numeric(predicted_train)
    print("train_statistics")
    MSE_train = mean((training_target - predicted_train)^2)
    Cor_train = cor(training_target, predicted_train, method = "pearson")
    scaled_predicted_target_train <- rescale(predicted_train, min_val = min_origin_target, max_val = max_origin_target)
    scaled_original_target_train <- rescale(training_target, min_val = min_origin_target, max_val = max_origin_target)
    MSE_backscaled_train = mean((scaled_original_target_train - scaled_predicted_target_train)^2)
    
    
    end_time <- Sys.time() # Record end time
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    
    accuracy1 = data.frame(
      Name = gene_name,
      MSE_test_NN2 = MSE1,
      Cor_test_NN2 = cor1, 
      MSE_S_test_NN2 = MSE_backscaled,
      MSE_train_NN2 = MSE_train,
      Cor_train_NN2 = Cor_train,
      MSE_S_train_NN2 = MSE_backscaled_train,
      val_error = validation_err,
      runtime_secs = runtime # Include runtime in seconds
    )
    # print("I saved the stat in a dataframe")
    
    acc_nested_list[[p]] <- (accuracy1)
    
    #saving the statistics
    saveRDS(acc_nested_list, "/projects/apog/work/models/1MB/new_MLP_seed/NN2/stat/NN2_seed_statistics_V2.RDS") #local
    
    
    print(gene_name)
  }, error = function(e) {cat("ERROR:", conditionMessage(e), "\n")})
  
  stopCluster(cl) 
}
