library(glmnet)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)
library(doSNOW)


eqtl.dir <- "~/ISM_Mini_eQTL_hg38"
models.dir <- "~/regression_segmentation_partition0_all_genes" 
feature.mat.dir <- "~/feature_counts_partition0_test_all_genes" 
segmentation.dir <- "~/segmentation_partition0" 
output.dir <- "~/STITCHIT_ism_all_eqtl_backscaled"
performance.file <- "~/performance_test_partition0_all_genes.tsv"

eqtl.files <- list.files(eqtl.dir, full.names = TRUE)

#load model performance  
	#load per gene eQTL interactions (bed format): chr | start | end | ihec_sample_id
        # extract sample ids: input ism
performance.df <- read.table(performance.file, header=TRUE, sep="\t")

num_cores <- 30

for (eqtl.file in eqtl.files){
  
  gene.id <- strsplit(basename(eqtl.file), "_")[[1]][1]
  print(gene.id)
  
  eqtl.df <- read.table(eqtl.file, header = TRUE, sep = "\t" )
  eqtl.samples <- unique(eqtl.df[,4])
  print(eqtl.samples)
  
  #load input model
  model.file <- file.path(models.dir, paste0("Elasticnet_Regression_Model_Segmentation_", gene.id, "_Pearson_10.RData"))
  if (!file.exists(model.file)){
    print(paste(gene.id, ": model does not exist"))
    next
  }
  load(model.file)
  print(paste("loaded model for gene id:", gene.id))
  
  # features without intercept
  features <- rownames(coef(elasticnet_model[[1]], s="lambda.min"))[-1]
  
  #predict gene expression: for each predicted interaction set corresponding feature to 0

  #load segmentation file (training data)
  segmentation.file <- file.path(segmentation.dir, paste0("Segmentation_",gene.id,"_Pearson_10.txt"))
  training.counts <- read.table(segmentation.file, header = TRUE, sep = "\t" )
  rownames(training.counts) <- training.counts$X
  training.counts <- training.counts[,-1]
  
  #scaling parameters STITCHIT
  scaled.training.counts <- scale(log2(training.counts+1), center=TRUE, scale=TRUE)
  training.means <- attr(scaled.training.counts,"scaled:center")
  training.sds <- attr(scaled.training.counts,"scaled:scale")
  scaled.training.counts <- data.frame(scaled.training.counts) 

  #load test data
  test.counts.file <- file.path(feature.mat.dir, paste0(gene.id, "_stitchit_feature_counts.txt.gz"))
  test.counts <- read.table(gzfile(test.counts.file), header = TRUE, sep = "\t")
  rownames(test.counts) <- test.counts$Sample_id
  test.counts <- test.counts[,-1]
  
  #scale test data
  log.test.counts <- log2(test.counts+1)
  scaled.test.counts <- data.frame(scale(log.test.counts, center = training.means, scale = training.sds))

  #split expression
  x.train <- as.matrix(scaled.training.counts[,-which(names(scaled.training.counts) == "Expression")]) 
  x.test <- as.matrix(scaled.test.counts[,-which(names(scaled.test.counts) == "Expression")]) 
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  sample.df.list <- foreach(sample  = eqtl.samples, .packages=c("glmnet","GenomicRanges", "tidyr", "stringr"), .final = function(x) setNames(x, eqtl.samples)) %dopar% {
	
    sample.dir <- file.path(output.dir, sample)
    if(!dir.exists(sample.dir)){
      dir.create(sample.dir)
      print(paste(sample.dir, "was created!"))
    }
    
    #check if sample is in test or train feature matrix
    if (sample%in%row.names(test.counts)){   
      x.vector <- x.test[which(rownames(scaled.test.counts) == sample),]
    } else {
      x.vector <- x.train[which(rownames(scaled.training.counts) == sample),]
    }
    predict.fit <- predict(elasticnet_model$model, x.vector, s="lambda.min")
	
    #ISP
    ism.df <- data.frame(feature.coords = character(), ism = numeric(), prediction = numeric(), stringsAsFactors = FALSE)
    i=1
    for (feature in features) {
      p.vector <- x.vector
      p.feature <- features[i]
      p.vector[i] <- 0.0
      ism.fit <- predict(elasticnet_model$model, p.vector, s="lambda.min")
      temp.df <- data.frame(feature.coords = p.feature, ism = ism.fit, prediction = predict.fit, stringsAsFactors = FALSE)
      ism.df <- rbind(ism.df,temp.df)
      i=i+1    
    } 
    
    final.df <- separate(ism.df, feature.coords, into = c("chrom", "start", "end"), sep = "\\.")
    colnames(final.df) <- c("chrom", "start", "end", "ism", "prediction")
	
    #backscaling
    expression.mean <- training.means[length(training.means)]
    expression.sd <- training.sds[length(training.sds)]
	    
    #TODO: backscale gene expression into original expression without log transformation
    log.ism <- as.numeric(final.df$ism) * expression.sd + expression.mean
    log.prediction <- as.numeric(final.df$prediction) * expression.sd + expression.mean

    backscaled.ism <- 2^(log.ism)-1
    backscaled.prediction <- 2^(log.prediction)-1
		    
    #remove pseudo-count effect
    backscaled.ism[ backscaled.ism == -1] <- 0
    backscaled.prediction[ backscaled.prediction == -1] <- 0

    backscaled.ism.df <- data.frame(final.df$chrom, final.df$start, final.df$end, backscaled.ism, backscaled.prediction)
    colnames(backscaled.ism.df) <- c("chrom","start","end","ism_backscaled","prediction_backscaled")
   
  }

write.output <- function(x){
	 if(!is.null(sample.df.list[[x]])){ 
	    sample.dir <- file.path(output.dir, x)
	    output.file <- paste0(sample.dir,"/", gene.id ,".txt")
	    #write performance header
	    gene.performance <- performance.df[which(performance.df$GeneID == gene.id),]
	    header <- paste("#pearson_r=", gene.performance$Pearson, "\t", "mse_backscaled=", gene.performance$backscaled_MSE, sep = "")
	    writeLines(header, con = output.file)
	    write.table(sample.df.list[[x]], output.file, sep = "\t", append=TRUE,  row.names = FALSE, col.names = TRUE, quote =FALSE)
	    system(paste("gzip", output.file))
	 }
  }

  lapply(names(sample.df.list), write.output)
  stopCluster(cl)
  
}

