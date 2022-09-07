library(plyr) 
library(boot)
# library(neuralnet)
library(keras)
# library(corrplot)
library(tensorflow)
library(kerasR) 
library(tidyverse)
library(tfruns)
library(MLmetrics)
library(parallel)
# library(foreach)
# library(doParallel)
library(snow)
cl<-makeCluster(4)

# setwd("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1")

set.seed(200)

acc_nested_list = list()
pr_nested_list = list()
alpha1 = 0.6
alpha2 = 0.5
r = 0

# accuracy1 = data.frame(
#   feature1 = c(),
#   MSE11 = c(),
#   MAPE11 = c(),
#   Core11 = c())

fileList <- list.files(path="/projects/apog/work/input/Model_sh/EpiRegio_Mats_chr1", pattern=".txt")

mainlist = list()
for (v in 1: 20) {

  mainlist[[v]] = read.delim(paste0("/projects/apog/work/input/Model_sh/EpiRegio_Mats_chr1/", fileList[v]), header = TRUE, sep = "\t")

}
 for(p in 1:length(mainlist)) {   


#first data
df = mainlist[[p]]
df2 = df[, 2:length(df)]

#Normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
maxmindf <- as.data.frame(lapply(df2, normalize))

attach(maxmindf)
df_norm<-as.matrix(maxmindf)


# Determine sample size and split data
ind <- sample(2, nrow(df_norm), replace=TRUE, prob=c(0.80, 0.20))
training <- df_norm[ind==1, 1:ncol(df_norm)-1]
test1 <- df_norm[ind==2, 1:ncol(df_norm)-1]
training_target <- df_norm[ind==1, ncol(df_norm)]
test1_target <- df_norm[ind==2, ncol(df_norm)]



#number of nodes in the first hidden layer
u1_1 = ceiling((1/2) * (ncol(training)+1))
u2_1 = ceiling(1* (ncol(training)+1))
u3_1 = ceiling((2/3) * (ncol(training)+1))
u4_1 = ceiling(2*(ncol(training)))

#number of nodes in the second hidden layer
u1_2 = ceiling((2/3) * (u1_1)+1)
u2_2 = ceiling((2/3) * (u2_1)+1)
u3_2 = ceiling((2/3) * (u3_1)+1)
u4_2 = ceiling((2/3) * (u4_1)+1)
# u5_2 = ceiling((1/2) * (u1_1)+1)
# u6_2 = ceiling((1/2) * (u2_1)+1)
# u7_2 = ceiling((1/2) * (u3_1)+1)
# u8_2 = ceiling((1/2) * (u4_1)+1)

temp = NULL


nodes1 = c(u1_1,u2_1, u3_1, u4_1)
nodes2 = c(u1_2,u2_2, u3_2, u4_2)
           # , u5_2, u6_2, u7_2, u8_2)
# activations = c("relu")
dropouts = c(0,0.05)
batchsizes = c(4, 8, 16) #its related to samples which are 28 in our training dataset
Epochs= 1000
learning_rates=c(0.001, 0.01)
layers = c(1,2)
dropouts2 = c(0,0.05)

best_loss = Inf
best_loss2 = Inf


#------------------------------------------------------




# cl <- makeCluster(mc <- getOption("cl.cores", 2))
#--------------------------------------------------------
for (b in 1: length(layers)) {
  
  if(layers[b] == 1){
    print("enter to the one hidden layer mode")
    for (m in 1: length(batchsizes)) {
      print("m loop")
      for (j in 1: length(dropouts)) {
        print("j loop")
        for (k in 1:length(learning_rates)) {
          print("k loop")
          clusterExport(cl=cl, list=c("batchsizes", "dropouts", "learning_rates", "nodes1","m","i", "j","k","training", "training_target"))
          print("before first parLapply")
          a <-  parLapply(cl, 1:length(nodes1), 
                          fun1 <- function(i){require(keras)
                            model <- keras_model_sequential()
                            model %>%
                              layer_dense(nodes1[i], activation = "relu", input_shape = c(dim(training)[2])
                              ) %>%
                              
                              layer_dropout(rate = dropouts[j]) %>%
                              
                              layer_dense(units=1, activation ="linear")
                            
                            #####c) Compiling the DNN model
                            model %>% compile(
                              loss = 'mse',
                              optimizer = optimizer_adam(learning_rates[k]),
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
                              epochs = 1000, 
                              batch_size = batchsizes[m],
                              shuffled=F,
                              validation_split = 0.2,
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
            
            if(a[[y]]$loss < best_loss) {
              best_loss = a[[y]]$loss
              a_hyper_nodes1 = nodes1[y]
              # print(hyper_nodes1)
              a_hyper_dropouts = dropouts[j]
              a_hyper_batchsizes = batchsizes[m]
              a_hyper_learning_rates = learning_rates[k]
              best_model = a[[y]]$model
              # name1 = paste0("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/","geneee___",1,".hdf5")
              # # save_model_hdf5(Best_model, "/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/name1.hdf5")
              # save_model_hdf5(a[[y]]$model, name1)
            }
            
          }
          # stop()
          
          
          
        }#k1
      }#j1
    }#m1
    
    

  } # if(layers)    
  
  
  else{
    print("enter to two hidden layers mode")
    for (m in 1: length(batchsizes)) {
      for (j in 1: length(dropouts)) {
        for (jj in 1: length(dropouts2)) {
          for (k in 1:length(learning_rates)) {
            
            for (i in 1: length(nodes1) ) {
              # for (f in 1: length(nodes2) ) {
              #   print("here is one f")
                
                clusterExport(cl=cl, list=c("nodes1", "nodes2","dropouts", "dropouts2", "learning_rates","batchsizes","m","i", "j","k","jj","training", "training_target"))
                print("before second parLapply")
                a2 <-  parLapply(cl, 1:length(nodes2), 
                                fun2 <- function(f){require(keras)
                
                
                model <- keras_model_sequential()
                model %>%
                  layer_dense(nodes1[i], activation = "relu", input_shape = c(dim(training)[2])
                  ) %>%
                  layer_dropout(rate = dropouts[j]) %>%
                  layer_dense(nodes2[f], activation = "relu")   %>%
                  
                  layer_dropout(rate = dropouts2[jj]) %>%
                  
                  layer_dense(units=1, activation ="linear")
                
                
                #####c) Compiling the DNN model
                model %>% compile(
                  loss = 'mse',
                  optimizer = optimizer_adam(learning_rates[k]),
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
                  epochs = 1000, 
                  batch_size = batchsizes[m],
                  shuffled=F,
                  validation_split = 0.2,
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
                    b_hyper_nodes2 = nodes2[w]
                    # print(hyper_nodes1)
                    b_hyper_nodes1 = nodes1[i]
                    # print(hyper_nodes1)
                    b_hyper_dropouts = dropouts[j]
                    b_hyper_batchsizes = batchsizes[m]
                    b_hyper_learning_rates = learning_rates[k]
                    b_hyper_dropouts2 = dropouts2[jj]
                    best_model2 = a2[[w]]$model
                    # name1 = paste0("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/","geneee___",1,".hdf5")
                    # # save_model_hdf5(Best_model, "/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/name1.hdf5")
                    # save_model_hdf5(a[[y]]$model, name1)
                  }
                  
                }
                
                

                
                
                
              # }#f
            }#i2
          }#k2
        }#jj
      }#j2
    }#m2
    #second loss
    
    
    
  }
  
}#b





print("the best mse and Hyperparameters are:")
print(best_loss)
print("--------")
print(best_loss2)
if(best_loss <= best_loss2){
  print("The model with  one hidden layer gives better result")
  #------------------------------------------
  #Best model with best paramethers
  Best_model<-keras_model_sequential()
  Best_model %>%
    
    layer_dense(units = a_hyper_nodes1, activation = "relu", input_shape = c(dim(training)[2])) %>%
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
    epochs = 1000, 
    batch_size = a_hyper_batchsizes,
    shuffled=F,
    validation_split = 0.2,
    verbose=0,
    # callbacks = list(early_stop, print_dot_callback)
  )
  Best_ModelFited
  Best_ModelFited%>% summary()
  
}else{
  print("The model with two hidden layers gives better result")
  #Best model with best paramethers
  Best_model<-keras_model_sequential()
  Best_model %>%
    layer_dense(b_hyper_nodes1, activation = "relu", input_shape = c(dim(training)[2])
    ) %>%
    layer_dropout(rate = b_hyper_dropouts) %>%
    layer_dense(b_hyper_nodes2, activation = "relu") %>%
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
    epochs = 1000, 
    batch_size = b_hyper_batchsizes,
    shuffled=F,
    validation_split = 0.2,
    verbose=0,
    # callbacks = list(early_stop, print_dot_callback)
  )
  Best_ModelFited
  Best_ModelFited%>% summary()
  
  
}

name1 = paste0("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/","gene___",p,".hdf5")
# save_model_hdf5(Best_model, "/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/models/name1.hdf5")
save_model_hdf5(Best_model, name1)





#test
Yhat = Best_model%>% predict(test1)
y_p = Yhat
y_p_tst = as.numeric(y_p)
#y_tst=y[tst_set]
plot(test1_target,y_p_tst)+
  abline(0,1)

# a = test1_target
# a  = cbind(a, y_p_tst)
# a
MSE1 = mean((test1_target - y_p_tst)^2)
MAPE1 = MAPE(test1_target,y_p_tst)
cor1 = cor(test1_target, y_p_tst)
# print("MAPE is")
# print(MAPE1)
# print("Correlation of predicted data is:")
# print(cor1)
# print("The MSE is:")
# print(MSE1)

dimm = ncol(mainlist[[p]])
accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
# accuracy1 = data.frame(
#   feature1 = dimm,
#   MSE11 = MSE1,
#   MAPE11 = MAPE1,
#   Core11 = cor1)
# accuracy1$feature1 = dimm
# accuracy1@MSE1 = MSE1

# Run for-loop over lists
acc_nested_list[[p]] <- (accuracy1)
name2 = paste0("/Users/shamim/Desktop/PhD/ML project/EpiRegio_Mats_chr1/acc_metrics",".RDS")
saveRDS(acc_nested_list, name2) #local
acc_metrics <- readRDS("~/Desktop/PhD/ML project/EpiRegio_Mats_chr1/acc_metrics.RDS")

corelations =  data.frame(do.call(rbind.data.frame, acc_metrics))
for (u in 1: length(corelations)) {
  if(alpha2 <= corelations[,4]){
    r = r+1
  }
}
print("The precentage of genes which passed the thereshould are")
print(r / p)
}#p
stopCluster(cl) 
