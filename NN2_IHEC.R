print("two hidden layers NN(MLP) with a search grid to find the best model which is parallized on a CPU")

library(tensorflow)
# install_keras()R
library(keras)

library(plyr) 
library(boot)
# install.packages("keras")
# library(keras)
print("after library")

print("after two new package")

library(tfruns)
print("after tfruns")
library(MLmetrics)
print("after MLmetrics")
library(parallel)
library(snow)
print("after snow")
#library(tidyverse)
library(reticulate)

# checkvector = c()
acc_nested_list = list()
pr_nested_list = list()

df = NULL

acc_nested_list = list()




fileList <- list.files(path="/projects/apog/work/input/IHEC_Activity", pattern=".txt.gz")
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
    
    
    cl<-makeCluster(36)
    
    
    if(file.exists(paste0("/projects/apog/work/models/IHEC_Activity_NN2/", gene_name, ".hdf5"))){
      print("it exist, so let's go to the next gene")
      
      
    }else{
      print("file not exist")
      df = mainlist1[[p]]
      
      V1 = var(df["Expression"]) * ((ncol(df)-1)/ncol(df))
      print("VAriance is:")
      V1 = V1[1,1]
      print(V1)
      if (V1 >= 10){
        print("inside the loop")
        
        print("inside the loop")
        print("dim of df befor removing:")
        print(dim(df))
        print("//////////////////////////////")
        df <- df[order(df$Sample),]
        df = df[, colSums(df != 0) > 0]
        
        print("dim of df after removing:")
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
        
        dropouts = c(0.2, 0.4)
        batchsizes = 32 #its related to samples which are 28 in our training dataset
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
                # b_hyper_nodes2 = nodes2[w]
                # print(hyper_nodes1)
                # b_hyper_nodes1 = nodes1[i]
                # print(hyper_nodes1)
                b_hyper_dropouts = dropouts[j]
                # b_hyper_batchsizes = batchsizes[m]
                b_hyper_learning_rates = learning_rates[w]
                b_hyper_dropouts2 = dropouts2[jj]
                best_model2 = a2[[w]]$model
                # name1 = paste0("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/","geneee___",1,".hdf5")
                # # save_model_hdf5(Best_model, "/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/name1.hdf5")
                # save_model_hdf5(a[[y]]$model, name1)
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
        
        
        
        
        
        print("the best mse and Hyperparameters are:")
        print(best_loss)
        print("--------")
        print(best_loss2)
        
        print("The model with two hidden layers gives better result")
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
        
        
        # }
        
        name1 = paste0("/projects/apog/work/models/IHEC_Activity_NN2/",gene_name,".hdf5")
        # save_model_hdf5(Best_model, "/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/name1.hdf5")
        save_model_hdf5(Best_model, name1)
        
        
        
        
        
        #test
        Yhat = Best_model%>% predict(test1)
        Yhat_t = Best_model%>% predict(training)
        y_p = Yhat
        y_p_tst = as.numeric(y_p)
        y_p_t = Yhat_t
        y_p_trn = as.numeric(y_p_t)
        
        # plot(test1_target,y_p_tst)+
        #   abline(0,1)ssh
        
        
        MSE1 = mean((test1_target - y_p_tst)^2)
        MAPE1 = MAPE(test1_target,y_p_tst)
        cor1 = cor(test1_target, y_p_tst)
        # cor2 = cor(training_target, y_p_trn)
        
        
        
        # accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
        dimm = ncol(mainlist1[[p]])
        # accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
        accuracy1 = data.frame(
          features = dimm,
          MSEs= MSE1,
          MAPEs = MAPE1,
          Corelations = cor1, 
          Var = V1,
          Name = gene_name,
          GeneInd= p)
        
        # Run for-loop over lists
        acc_nested_list[[p]] <- (accuracy1)
        name2 = paste0("/projects/apog/work/models/IHEC_Activity_NN2/IHEC_Activity_NN2_002_002",".RDS")
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
        name3 = paste0("//projects/apog/work/models/IHEC_Activity_NN2/IHEC_Activity_NN2_002_002",".RDS")
        # name3 = paste0("~/Desktop/PhD/ML_project/ml_scripts/IHEC_test/IHEC1",".RDS")
        # print(name2)
        saveRDS(acc_nested_list, name3) #local
        
        
      }
    }
    
    # }
  },  error=function(e){
    print(paste("Error occurred at gene", p)); message(e)
    
    
  }
  
  
  )
  
  
  stopCluster(cl) 
}#p