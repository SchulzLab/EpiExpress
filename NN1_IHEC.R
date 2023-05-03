print("One hidden layer Neural Network(MLP) with a search grid which is parallized on a CPU")
#It saves best model on the train dataset then test it on the test dataset.
#We find the best model based on search grid

library(tensorflow)
# install_keras()R
library(keras)
library(plyr) 
library(boot)
print("after library")

library(tfruns)
print("after tfruns")
library(MLmetrics)
print("after MLmetrics")
library(parallel)
library(snow)
print("after snow")
#library(tidyverse)
library(reticulate)

print("2") 
# checkvector = c()
acc_nested_list = list()
pr_nested_list = list()

df = NULL

acc_nested_list = list()

fileList <- list.files(path="/projects/apog/work/input/IHEC_Activity", pattern=".txt.gz") #21582 genes


for(p in 1:length(fileList)) { 
  
  tryCatch({
    
    # mainlist = list()
    mainlist1 = list()
    
    print(p)
    gene_name <- strsplit(fileList[p], "\\.")[[1]][1]
    
    mainlist1[[p]] = read.table(paste0("/projects/apog/work/input/IHEC_Activity/", fileList[p]), header = TRUE, sep = "\t"
                                # , row.names = 1
    )
    
    print("Data is ready")
    
    
    cl <-makeCluster(32)
    
    # checkvector[p] = mainlist1[[p]]
    df = mainlist1[[p]]
    
    V1 = var(df["Expression"]) * ((ncol(df)-1)/ncol(df))
    print("VAriance is:")
    print(V1)
    if (V1 >= 10){
      print("inside the loop")
      
      
      print("//////////////////////////////")
      df <- df[order(df$Sample),]
      df = df[, colSums(df != 0) > 0]
      print(dim(df))
      
      row.names = 1
      df2 = df[, 2:length(df)]
      
      #Normalization
      normalize <- function(x) {
        if(max(x) == min(x)){
          return(0)
        }
        return ((x - min(x)) / (max(x) - min(x)))
      }
      maxmindf <- as.data.frame(lapply(df2, normalize))
      
      # attach(maxmindf)
      df_norm<-as.matrix(maxmindf)
      
      set.seed(200)
      # Determine sample size and split data
      ind <- sample(2, nrow(df_norm), replace=TRUE, prob=c(0.80, 0.20))
      training <- df_norm[ind==1, 1:ncol(df_norm)-1]
      test1 <- df_norm[ind==2, 1:ncol(df_norm)-1]
      training_target <- df_norm[ind==1, ncol(df_norm)]
      test1_target <- df_norm[ind==2, ncol(df_norm)]
      
      temp = NULL
      # nodes1 = 150
      dropouts = c(0.2, 0.4)
      batchsizes = 32 #its related to samples which are 28 in our training dataset
      Epochs= 1200
      learning_rates=c(0.001, 0.01)
      # layers = c(1,2)
      dropouts2 = c(0.2, 0.4)
      
      best_loss = Inf
      best_loss2 = Inf
      
      # cl <- makeCluster(mc <- getOption("cl.cores", 2))
      #--------------------------------------------------------
      
      
      
      
      # for (m in 1: length(batchsizes)) {#m
      
      for (j in 1: length(dropouts)) {#j
        print("j loop")
        clusterExport(cl, list=c("dropouts", "learning_rates","j","training", "training_target"))
        print("before first parLapply")
        a <-  parLapply(cl, 1:length(learning_rates),
                        fun1 <- function(i){require(keras)
                          model <- keras_model_sequential()
                          model %>%
                            layer_dense(100, activation = "relu",
                                        # kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01),
                                        input_shape = c(dim(training)[2])
                            ) %>%
                            
                            
                            layer_dropout(rate = dropouts[j]) %>%
                            # regularizer_l1_l2(l1 = 0.01, l2 = 0.01) %>
                            layer_dense(units=1, activation ="linear")
                          
                          #####c) Compiling the DNN model
                          model %>% compile(
                            loss = 'mse',
                            optimizer = optimizer_adam(learning_rates[i]),
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
        
        
        
        for(y in 1:length(a)){
          print("hello")
          if(a[[y]]$loss < best_loss) {
            best_loss = a[[y]]$loss
            # a_hyper_nodes1 = nodes1[y]
            # print(hyper_nodes1)
            a_hyper_dropouts = dropouts[j]
            # a_hyper_batchsizes = batchsizes[m]
            a_hyper_learning_rates = learning_rates[y]
            print(a_hyper_learning_rates)
            best_model = a[[y]]$model
            
          }
          
        }
        # stop()
        
        
        
        # }#k1
      }#j1
      # }#m1
      
      
      
      # } # if(layers)
      
      
      
      
      #------------------------------------------
      #Best model with best paramethers
      Best_model<-keras_model_sequential()
      Best_model %>%
        
        layer_dense(units = 100, activation = "relu",  
                    # kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.01), 
                    input_shape = c(dim(training)[2])) %>%
        layer_dropout(rate = a_hyper_dropouts) %>%
        
        layer_dense(units=1, activation ="linear")
      
      
      #####c) Compiling the DNN model
      Best_model %>% compile(
        loss = 'mse',
        optimizer = optimizer_adam(a_hyper_learning_rates),
        metrics = c('mse'))
      Best_model
      # }
      
      Best_model %>% summary()
      
      
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
      
      name1 = paste0("/projects/apog/work/models/IHEC_Activity_NN1_/",gene_name,".hdf5")
      save_model_hdf5(Best_model, name1)
      
      #test
      Yhat = Best_model%>% predict(test1)
      Yhat_t = Best_model%>% predict(training)
      y_p = Yhat
      y_p_t = Yhat_t
      y_p_tst = as.numeric(y_p)
      y_p_trn = as.numeric(y_p_t)
      
      
      MSE1 = mean((test1_target - y_p_tst)^2)
      MAPE1 = MAPE(test1_target,y_p_tst)
      cor1 = cor(test1_target, y_p_tst)
      
      dimm = ncol(mainlist1[[p]])
      # accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
      accuracy1 = data.frame(
        features = dimm,
        MSEs= MSE1,
        MAPEs = MAPE1,
        Corelations = cor1, 
        Var = V1,
        Name = gene_name,
        GeneInd= p
      )
      
      
      # Run for-loop over lists
      acc_nested_list[[p]] <- (accuracy1)
      name2 = paste0("/projects/apog/work/models/IHEC_Activity_NN1_/IHEC_Activity_NN1_002",".RDS")
      # name2 = paste0("~/Desktop/PhD/ML_project/ml_scripts/IHEC_test/IHEC_Activity_NN1_001",".RDS")
      print(name2)
      saveRDS(acc_nested_list, name2) #local
      
      
    }else{
      print("The variance of Expression Column for this gene is less than 2 so I put it out")
      print(p)
      
      dimm = ncol(mainlist1[[p]])
      # accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
      accuracy1 = data.frame(
        features = dimm,
        MSEs= NA,
        MAPEs = NA,
        Corelations = NA, 
        Var = V1,
        Name = gene_name,
        GeneInd= p
      )
      
      
      # Run for-loop over lists
      acc_nested_list[[p]] <- (accuracy1)
      name3 = paste0("/projects/apog/work/models/IHEC_Activity_NN1_/IHEC_Activity_NN1_002",".RDS")
      # name3 = paste0("~/Desktop/PhD/ML_project/ml_scripts/IHEC_test/IHEC1",".RDS")
      # print(name2)
      saveRDS(acc_nested_list, name3) #local
      
      
    }
    
    
  },  error=function(e){
    print(paste("Error occurred at gene", p)); message(e)
    
    
  }
  
  
  )
  
  
  stopCluster(cl) 
}#p

