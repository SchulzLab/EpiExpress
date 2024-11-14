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
set.seed(200)
tensorflow::tf$random$set_seed(200)

acc_nested_list = list()
pr_nested_list = list()
acc_nested_list2 = list()


file_path = ("path/to/gene/list/gene_list.txt")
file_contents <- readLines(file_path)

df = NULL
acc_nested_list = list()

validation_split = 0.3

for(p in file_contents) { 
  set.seed(200)
  tryCatch({
    
    gene_name <- p
    
    if(file.exists(paste0("path/to/input/files/", gene_name, ".txt.gz"))){
      print("TRUE")
      df = read.table(paste0("path/to/input/files/", gene_name, ".txt.gz"), header = TRUE, sep = "\t")
    } else{
      print("NOTTRUE")

    }
    
    mainlist1 = list()
    print(gene_name)
    cl <-makeCluster(50)
    clusterEvalQ(cl, { #keep reproducibility
      tensorflow::tf$random$set_seed(200)
    })
    
    df <- df[order(df$Sample),] #ordering the sample names
    rownames(df) = NULL
    df = df[, colSums(df != 0) > 0] #if all values of a column is 0 Ill remove that column
    df_samples = df[,1] 
    df2 = df[, 2:length(df)] 
    df2 = log2(df2 + 1) #whole data after log2, 
    df3 = cbind(df_samples, df2)
    
    test_samples <- "path/to/test/samples/partition0_test.csv"
    train_samples <- "path/to/train/samples/partition0_train.csv"
    
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
    learning_rates=c(0.001, 0.01)
    layers = c(1,2)
    dropouts2 = c(0.2, 0.4)
    
    best_loss = Inf
    best_loss2 = Inf
    
    
    #------------------------------------------------------
    
    # else{
    print(" two hidden layers mode")
    # for (m in 1: length(batchsizes)) {
    for (j in 1: length(dropouts)) {
      for (jj in 1: length(dropouts2)) {
        
        clusterExport(cl=cl, list=c("dropouts","learning_rates", "dropouts2", "j","jj","training", "training_target"))
        print("before second parLapply")
        a2 <-  parLapply(cl, 1:length(learning_rates),
                         fun2 <- function(f){require(keras)
                           tensorflow::tf$random$set_seed(200) #to controll randomness
                           
                           model <- keras_model_sequential()
                           model %>%
                             layer_dense(100, activation = "relu"
                                         # , kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01)
                                         ,input_shape = c(dim(training)[2])
                             ) %>%
                             
                             layer_dropout(rate = dropouts[j]) %>%
                             layer_dense(50, activation = "relu"
                                         # , kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01)
                             ) %>%
                             
                             layer_dropout(rate = dropouts2[jj]) %>%
                             
                             layer_dense(units=1, activation ="linear")
                           
                           
                           #####c) Compiling the DNN model
                           model %>% compile(
                             loss = 'mse',
                             optimizer = optimizer_adam(learning_rates[f]),
                             metrics = c('mse'))
                           model
                           # }
                           print_dot_callback <- callback_lambda(
                             on_epoch_end = function(epoch, logs) {
                               if (epoch %% 100 == 0) cat("\n")
                               cat(".")})
                           
                           early_stop <- callback_early_stopping(monitor = "val_loss", mode='min',patience =30)
                           # model<-build_model()
                           # model %>% summary()
                           
                           
                           ###########d) Fitting the DNN model#################
                           
                           model1<-model %>% fit(
                             training,
                             training_target,
                             epochs = 1200,
                             batch_size = 32,
                             shuffled=F,
                             validation_split = 0.3,
                             verbose=0,
                             callbacks = list(early_stop, print_dot_callback)
                           )
                           
                           
                           temp_loss = mean(model1$metrics$val_loss)
                           
                           
                           return(list(model=model1, loss = temp_loss
                                       #  ,o_dropouts = dropouts[j],
                                       # o_batchsizes = batchsizes[m], o_learning_rates = learning_rates[k],
                           ))
                         })
        
        
        for(w in 1:length(a2)){
          
          if(a2[[w]]$loss < best_loss2) {
            best_loss2 = a2[[w]]$loss
            
            b_hyper_dropouts = dropouts[j]
            # b_hyper_batchsizes = batchsizes[m]
            b_hyper_learning_rates = learning_rates[w]
            b_hyper_dropouts2 = dropouts2[jj]
            best_model2 = a2[[w]]$model
            
          }
          
        }
        
        
        
        
        
        
        # }#f
        # }#i2
        # }#k2
      }#jj
    }#j2
    # }#m2
    #second loss
    
    
    
    # }
    
    # }#b
    
    
    validation_err = best_loss2
    #Best model with best paramethers
    Best_model<-keras_model_sequential()
    Best_model %>%
      layer_dense(100, activation = "relu", 
                  # kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01), 
                  input_shape = c(dim(training)[2])
      ) %>%
      layer_dropout(rate = b_hyper_dropouts) %>%
      layer_dense(50, activation = "relu",  
                  # kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01)
      ) %>%
      layer_dropout(rate = b_hyper_dropouts2) %>%
      
      layer_dense(units=1, activation ="linear")
    
    
    #####c) Compiling the DNN model
    Best_model %>% compile(
      loss = 'mse',
      optimizer = optimizer_adam(b_hyper_learning_rates),
      metrics = c('mse'))
    Best_model
    # }
    print("---------------------------------------")
    Best_model %>% summary()
    print("---------------------------------------")
    
    ###########d) Fitting the DNN model#################
    Best_ModelFited<-Best_model %>% fit(
      training, 
      training_target,
      epochs = 1200, 
      batch_size = 32,
      shuffled=F,
      validation_split = 0.3,
      verbose=0,
      # callbacks = list(early_stop, print_dot_callback)
    )
    Best_ModelFited
    Best_ModelFited%>% summary()
    print("---------------------------------------")
    
    # }
    
    name1 = paste0("path/to/ottput/file/for/saving/best/model/",gene_name,".hdf5")
    save_model_hdf5(Best_model, name1)
    print("the best model is saved")
    
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
    
    
    
    
    accuracy1 = data.frame(
      Name = gene_name,
      MSE_test_NN2 = MSE1,
      Cor_test_NN2 = cor1, 
      MSE_S_test_NN2 = MSE_backscaled,
      MSE_train_NN2 = MSE_train,
      Cor_train_NN2 = Cor_train,
      MSE_S_train_NN2 = MSE_backscaled_train,
      val_error = validation_err 
    )
    print("I saved the stat in a dataframe")

    acc_nested_list[[p]] <- (accuracy1)

    #saving the statistics
    saveRDS(acc_nested_list, "path/to/output/for/saving/statistics/NN2_statistics.RDS") #local

    
    
    
    # }
  },  error=function(e){
    print(paste("Error occurred at gene", p)); message(e)
    
    
  }
  
  
  )
  
  
  stopCluster(cl) 
}#p



print("done")

