### Import libraries

library(ggplot2)
library(rsample) 
install.packages("randomForest")
library(randomForest) # basic implementation
library(caret)        # an aggregator package for performing many machine learning models
library(Metrics)
library(MLmetrics)

set.seed(200)
fileList <- list.files(path="/projects/apog/work/input/MPAL_Mats/", pattern=".txt")

mainlist = list()
for (v in 1: 3000) {
  print(v)
  mainlist[[v]] = read.delim(paste0("/projects/apog/work/input/MPAL_Mats//", fileList[v]), header = TRUE, sep = "\t")
  
}

name1 = paste0("/projects/apog/work/input/MPAL_Mats/main_data",".RDS")
saveRDS(mainlist, name1)
mainlist <- readRDS("/projects/apog/work/input/MPAL_Mats/main_data.RDS")
print("after reading the dataset")
acc_nested_list = list()
for(j in 1522:length(mainlist)){
  print(j)
  dataset1 = mainlist[[j]]
  if(ncol(dataset1) >= 339){
    # print(dim(mainlist[[j]]))}
    # Splitting data into Train and Test
    # dataset = dataset1[, 2:length(dataset1)]
    # ind <- sample(2, nrow(dataset), replace=TRUE, prob=c(0.70, 0.30))
    # train_data <- dataset[ind==1, 1:ncol(dataset)]
    # test_data <- dataset[ind==2, 1:ncol(dataset)]
    
    
    X = dataset1[,2: ncol(dataset1) - 1]
    Y = dataset1[ ,ncol(dataset1)]
    index <- createDataPartition(Y, p=0.80, list=FALSE)
    X_train <- X[index, ]
    X_test <- X[-index, ]
    y_train <- Y[index]
    y_test<-Y[-index]
    
    
    # Train the model 
    regr <- randomForest::randomForest(x = X_train, y = y_train
                         , mtry = (ncol(X_train)/3)
                         , ntree = 501)
    
    # Make prediction
    predictions <- predict(regr, X_test)
    df1 = cbind(predictions, y_test)
    Cor1 = cor(predictions, y_test)
    print(Cor1)
    print(j)
    # print(paste0('MAE: ' , mae(y_test,predictions) ))
    # print(paste0('MSE: ' ,caret::postResample(predictions , y_test)['RMSE']^2 ))
    # print(paste0('R2: ' ,caret::postResample(predictions , y_test)['Rsquared'] )) #R-squared scores
    MSE1 = mean((predictions -  y_test)^2)
    MAPE1 = MAPE(predictions, y_test)
    dimm = ncol(mainlist[[j]])
    # accuracy1 = c(dimm ,MSE1, MAPE1, cor1) #for each gene/model we store these 3 variables
    accuracy1 = data.frame(
      features = dimm,
      MSEs= MSE1,
      MAPEs = MAPE1,
      Corelations = Cor1,
      Gene1 = j)
    
    acc_nested_list[[j]] <- (accuracy1)
    name2 = paste0("/projects/apog/work/input/MPAL_Mats/RF_acc7",".RDS")
    saveRDS(acc_nested_list, name2) #local
    # acc1 <- readRDS("/Users/shamim/Desktop/PhD/ML project/MPAL_Mats/acc_metrics1.RDS")
  }
}

