#SHAP values for the Random Forest Method on IHEC dataset

library(randomForest)
library(randomForestExplainer)
library(keras)
library(fastshap)
library(readr)
library(tensorflow)
library(keras)

library(plyr) 
library(boot)
library(tfruns)
library(MLmetrics)
library(parallel)
library(snow)
library(tidyverse)
library(reticulate)

# # #preparing the validation file
# # 
# # mytext <- readLines("Validated_Genes.txt")
# # 
# # # Remove everything after the dot on each line
# # mytext_clean <- gsub("\\..*", "", mytext)
# # 
# # # Save the cleaned text as a new file in txt format
# # writeLines(mytext_clean, "mytextfile_clean.txt")

lst =list()
myCon = file(description = "/projects/apog/work/models/SHAP_RF/Validated_Genes_clean.txt", open="r", blocking = TRUE)
myCon <- readLines(myCon)
fileList <- list.files(path="/projects/apog/work/input/IHEC_Activity/data1", pattern=".txt.gz")

for(gene_name in myCon)
{
  print(gene_name)
  
     my_model <- readRDS(paste0("/projects/apog/work/models/test_RF_grid/RF_NEW_Format/", gene_name,".RDS"))
    my_data = read.table(paste0("/projects/apog/work/input/IHEC_Activity/data1/", gene_name, ".txt.gz"), header = TRUE, sep = "\t" )
    df = my_data
    df <- df[order(df$Sample),]
    df = df[, colSums(df != 0) > 0]
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
    #SHAP part
    print("data is ready")
    train_data = cbind(training,training_target)
    print("Hello")
    exp_forest <- fastshap::explain(my_model,X = subset(train_data, select = -training_target),pred_wrapper= predict,nsim = 1)
    # autoplot(exp_forest, type = "contribution")  # explain first row of X
    # autoplot(exp_forest)  # Shapley-based importance plot
    #save the found SHAP matrix with gene's name
    
    lst[[gene_name]] = exp_forest
    name1 = paste0("/projects/apog/work/models/SHAP_RF/SHAP_List1",".RDS")
    saveRDS(lst,name1)
  }
  print("END of the Validation File!")