#HAWAS Analysis
#load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("airway")
library(airway)
library(tidyverse)
########################################################### shared genes between kept genes, Binned-CNN and CRE-RF #############################################################################################################################
#28180 kept genes
kept_genes <- read.table("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/kept_genes.txt", header = FALSE, skip = 1)[, 2]
print(length(kept_genes))

#genes from RF
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled.RDS")
genes_from_RF <- rownames(health_Luekemia_scaled)
print(length(genes_from_RF))

#genes from CNN
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/CNN/CNN_Luekemia_backscaled_healthy.RDS")
genes_from_CNN <- rownames(health_Luekemia_scaled)
print(length(genes_from_CNN))


#finding common genes between kept_genes, Rf and CNN
shared_genes1 <- intersect(genes_from_RF, genes_from_CNN)
shared_genes2 <- intersect(shared_genes1, kept_genes)
print(length(shared_genes2)) #27,385 genes are shared

new_kept_genes <- shared_genes2
################################################################################################################################################################################################################################################

######################################################## Deseq on Real data ########################################################################################################################################################################################
#this part is on real RNA expression, we do Deseq2 on 22 samples from adult because we removed 2 far outliers also we removed one unknown sample from sex

kept_genes <- new_kept_genes #new kept genes
#metadata file for IHEC data
metadata_Luekemia <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/IHEC_metadata_harmonization.v1.1_BCellLeukemiaSamples.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
#picking desired colnames from the metadata:
col = c("epirr_id_without_version", "EpiClass_pred_Sex", "EpiClass_pred_Life_stage", "harmonized_sample_disease_high", "RNA_available", "project", "harmonized_tissue_type")
metadata_Luekemia = metadata_Luekemia[,col]
colnames(metadata_Luekemia) = c("sample", "sex", "age", "health_status", "RNA_available", "project", "harmonized_tissue_type")
#filtering samples: adult, those that we have real RNA for them, removing outliers:
metadata_Luekemia$health_status <- gsub('Healthy/None', 'Healthy', metadata_Luekemia$health_status)
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$age == "adult",]
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$RNA_available == "True",]
dim(metadata_Luekemia)
samples = metadata_Luekemia$sample
outliers = c("IHECRE00000785",  "IHECRE00000234","IHECRE00000298", "IHECRE00000732")
samples = setdiff(samples, outliers) #removing outliers

#picking healthy samples:
healthy_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Healthy",]
healthy_sample <- healthy_df$samples
#picking disease samples:
disease_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Cancer",]
disease_sample <- disease_df$samples

#loading real data
real_RNA_data <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/genes_expected_count_DESeq2_H3K27acFormatted.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
rownames(real_RNA_data) = real_RNA_data[,1]
real_RNA_data = real_RNA_data[,-1]
#picking new kept genes
real_RNA_data_filtered <- real_RNA_data[kept_genes,]
real_RNA_data_filtered <- t(real_RNA_data_filtered)
#picking defined samples:
count_data <- real_RNA_data_filtered[rownames(real_RNA_data_filtered) %in% samples,]
dim(count_data)

# read in sample info
metadata_Luekemia$harmonized_tissue_type[metadata_Luekemia$harmonized_tissue_type == ""] <- "unknown"
colData <- metadata_Luekemia[, c("sample", "health_status", "sex", "project", "harmonized_tissue_type")]
#preparing the data for deseq analysis
rownames(colData) <- colData$sample
count_data <- t (count_data)
#making sure the raw names in colData matches to column in countsData
all(colnames(count_data) %in% rownames(colData))
#same order
all(colnames(count_data) == rownames(colData))
colData <- colData[match(colnames(count_data), rownames(colData)),]
all(colnames(count_data) == rownames(colData))

# to make integer (for doing deseq2)
# Corrected custom rounding function based on the .5 threshold
custom_round <- function(x) {
  # Convert negative values to 0
  x[x < 0] <- 0
  # Extract the integer part and the fractional part
  integer_part <- floor(x)
  fractional_part <- x - integer_part
  # Apply rounding based on a comparison with 0.5
  x_rounded <- ifelse(fractional_part >= 0.5, ceiling(x), floor(x))
  return(x_rounded)
}

# Apply the corrected custom rounding function to the count data
count_data2 <- custom_round(count_data)

# Step2: construct a DeSeqDataset object
dds <- DESeqDataSetFromMatrix(countData = count_data2,
                              colData = colData,
                              design = ~ sex + health_status 
)

#pre-filtering: removig rows with low gene counts (this is recomended not required)
#for example, keeping genes which have higher than 10 counts in all cell lines / conditions
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set a factor level (conditions)
relevel(dds$health_status, ref = "Healthy")
# Estimate size factors for normalization
sizeFactors(dds) = 1
#we can see all size factors are 1
sizeFactors(dds)
# Estimate dispersion
dds <- estimateDispersions(dds)
# Fit the GLM (generalized linear model)
dds <- nbinomWaldTest(dds)
# Get results of differential expression analysis
res2 <- results(dds)

#Exploring the results:
results(dds, alpha = 0.05) #cut off is 0.05
sum(res2$padj < 0.05, na.rm=TRUE )
a2 = as.data.frame(res2)
sig_genes_deseq_adult_df <- a2[ a2$padj < 0.05,]
sig_genes_deseq_adult_df <- na.omit(sig_genes_deseq_adult_df)
sig_genes_deseq_adult_df <- sig_genes_deseq_adult_df %>%  filter(abs(log2FoldChange)> 0.585) #based on litretures we filter logfold change
sig_genes_deseq_adult_real <- rownames(sig_genes_deseq_adult_df)
#saving the result:
writeLines(sig_genes_deseq_adult_real, "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_real.txt")
# MA and dispersion plots
plotMA(res2)
################################################################################################################################################################################################################################################

###################################################################### Deseq on RF Predicted values #########################################################################################################################################################################
#siome preprocessing and filtering for samples has been don ein the above part (real data part)
# preparing disease RF (picking samples and genes)
disease_Luekemia_scaled <- readRDS("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/disease_Luekemia_scaled.RDS")
disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
disease_Luekemia_scaled <- disease_Luekemia_scaled [rownames(disease_Luekemia_scaled) %in% samples,]
disease_Luekemia_scaled <- disease_Luekemia_scaled[,kept_genes]

#preparing healthy RF (picking samples and genes)
health_Luekemia_scaled <- readRDS("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled.RDS")
health_Luekemia_scaled <- t(health_Luekemia_scaled)
health_Luekemia_scaled <- health_Luekemia_scaled [rownames(health_Luekemia_scaled) %in% samples,]
health_Luekemia_scaled <-  health_Luekemia_scaled [, kept_genes]

#making count data for the deseq analysis
count_data <- rbind(disease_Luekemia_scaled,health_Luekemia_scaled)

# read in sample info
metadata_Luekemia$harmonized_tissue_type[metadata_Luekemia$harmonized_tissue_type == ""] <- "unknown"
colData <- metadata_Luekemia[, c("sample", "health_status", "sex", "project", "harmonized_tissue_type")]
rownames(colData) <- colData$sample
count_data <- t (count_data)
#making sure the raw names in colData matches to column in countsData
all(colnames(count_data) %in% rownames(colData))
#same order
all(colnames(count_data) == rownames(colData))
colData <- colData[match(colnames(count_data), rownames(colData)),]
all(colnames(count_data) == rownames(colData))

# to make integer
# Corrected custom rounding function based on the .5 threshold
custom_round <- function(x) {
  # Convert negative values to 0
  x[x < 0] <- 0
  # Extract the integer part and the fractional part
  integer_part <- floor(x)
  fractional_part <- x - integer_part
  # Apply rounding based on a comparison with 0.5
  x_rounded <- ifelse(fractional_part >= 0.5, ceiling(x), floor(x))
  return(x_rounded)
}

# Apply the corrected custom rounding function to the count data
count_data2 <- custom_round(count_data)
# Step2: construct a DeSeqDataset object
dds <- DESeqDataSetFromMatrix(countData = count_data2,
                              colData = colData,
                              design = ~ sex + health_status 
)

#pre-filtering: removig rows with low gene counts (this is recomended not required)
#for example, keeping genes which have higher than 10 counts in all cell lines / conditions
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# set a factor level (conditions)
relevel(dds$health_status, ref = "Healthy")
sizeFactors(dds) = 1
#we can see all size factors are 1
sizeFactors(dds)
# Estimate dispersion
dds <- estimateDispersions(dds, )
# Fit the GLM (generalized linear model)
dds <- nbinomWaldTest(dds)
# Get results of differential expression analysis
res2 <- results(dds)
# Explore the result
results(dds, alpha = 0.05)
sum(res2$padj < 0.05, na.rm=TRUE )
a2 = as.data.frame(res2)
sig_genes_deseq_adult_df <- a2[ a2$padj < 0.05,]
sig_genes_deseq_adult_df <- na.omit(sig_genes_deseq_adult_df)
sig_genes_deseq_adult_df <- sig_genes_deseq_adult_df %>%  filter(abs(log2FoldChange)> 0.585) #based on litretures we filter it 

sig_genes_deseq_adult_RF <- rownames(sig_genes_deseq_adult_df)
writeLines(sig_genes_deseq_adult_RF, "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_RF.txt")
# MA and dispersion plots
plotMA(res2)
################################################################################################################################################################################################################################################

########################################################################### Deseq2 on CNN predicted values ###################################################################################################################################################################

#loading the data: #disease CNN
disease_Luekemia_scaled <- readRDS("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/disease_Luekemia_scaled_correct_scale_CNN.RDS") #scale_corrected_data
#disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
disease_Luekemia_scaled <- disease_Luekemia_scaled [rownames(disease_Luekemia_scaled) %in% samples,]
disease_Luekemia_scaled <- disease_Luekemia_scaled[,kept_genes]
#loading the data: #Healthy CNN
health_Luekemia_scaled <- readRDS("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled_correct_scale_CNN.RDS") #scale_corrected_data
#health_Luekemia_scaled <- t(health_Luekemia_scaled)
health_Luekemia_scaled <- health_Luekemia_scaled [rownames(health_Luekemia_scaled) %in% samples,]
health_Luekemia_scaled <- health_Luekemia_scaled[,kept_genes]
#count matrix for deseq analysis
count_data <- rbind(disease_Luekemia_scaled,health_Luekemia_scaled)

# read in sample info
metadata_Luekemia$harmonized_tissue_type[metadata_Luekemia$harmonized_tissue_type == ""] <- "unknown"
colData <- metadata_Luekemia[, c("sample", "health_status", "sex", "project", "harmonized_tissue_type")]
rownames(colData) <- colData$sample
count_data <- t (count_data)
#making sure the raw names in colData matches to column in countsData
all(colnames(count_data) %in% rownames(colData))
#same order
all(colnames(count_data) == rownames(colData))
colData <- colData[match(colnames(count_data), rownames(colData)),]
all(colnames(count_data) == rownames(colData))

# to make integer
# Corrected custom rounding function based on the .5 threshold
custom_round <- function(x) {
  # Convert negative values to 0
  x[x < 0] <- 0
  # Extract the integer part and the fractional part
  integer_part <- floor(x)
  fractional_part <- x - integer_part
  # Apply rounding based on a comparison with 0.5
  x_rounded <- ifelse(fractional_part >= 0.5, ceiling(x), floor(x))
  return(x_rounded)
}

# Apply the corrected custom rounding function to the count data
count_data2 <- custom_round(count_data)

# Step2: construct a DeSeqDataset object
dds <- DESeqDataSetFromMatrix(countData = count_data2,
                              colData = colData,
                              design = ~ sex + health_status 
)

#pre-filtering: removig rows with low gene counts (this is recomended not required)
#for example, keeping genes which have higher than 10 counts in all cell lines / conditions
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set a factor level (conditions)
relevel(dds$health_status, ref = "Healthy")
# Estimate size factors for normalization
sizeFactors(dds) = 1
sizeFactors(dds)
# Estimate dispersion
dds <- estimateDispersions(dds)
# Fit the GLM (generalized linear model)
dds <- nbinomWaldTest(dds)
# Get results of differential expression analysis
res2 <- results(dds)
results(dds, alpha = 0.05)
sum(res2$padj < 0.05, na.rm=TRUE )
a2 = as.data.frame(res2)
sig_genes_deseq_adult_df <- a2[ a2$padj < 0.05,]
sig_genes_deseq_adult_df <- na.omit(sig_genes_deseq_adult_df)
sig_genes_deseq_adult_df <- sig_genes_deseq_adult_df %>%  filter(abs(log2FoldChange)> 0.585)
sig_genes_deseq_adult_CNN <- rownames(sig_genes_deseq_adult_df)

writeLines(sig_genes_deseq_adult_CNN, "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_CNN.txt")
# MA and dispersion plots
plotMA(res2)
################################################################################################################################################################################################################################################

################################################### scatterplot for log fold change for real vs RF (for sig genes from real data) ########################################################################
#kept_genes are already loaded (first part of the script)
#RF part
kept_genes <- new_kept_genes
metadata_Luekemia <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/IHEC_metadata_harmonization.v1.1_BCellLeukemiaSamples.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
col = c("epirr_id_without_version", "EpiClass_pred_Sex", "EpiClass_pred_Life_stage", "harmonized_sample_disease_high", "RNA_available", "project", "harmonized_tissue_type")
metadata_Luekemia = metadata_Luekemia[,col]
colnames(metadata_Luekemia) = c("sample", "sex", "age", "health_status", "RNA_available", "project", "harmonized_tissue_type")
metadata_Luekemia$health_status <- gsub('Healthy/None', 'Healthy', metadata_Luekemia$health_status)
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$age == "adult",]
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$RNA_available == "True",]
dim(metadata_Luekemia)
samples = metadata_Luekemia$sample
outliers = c("IHECRE00000785",  "IHECRE00000234","IHECRE00000298", "IHECRE00000732")
samples = setdiff(samples, outliers)
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$sample %in% samples,]

# samples <- setdiff(samples,"IHECRE00000298")
healthy_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Healthy",]
healthy_sample <- healthy_df$sample
disease_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Cancer",]
disease_sample <- disease_df$sample

#loading singnificant genes from real data
significant_gene_names_real= readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_real.txt")
length(significant_gene_names_real)

#load healthy part
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled.RDS")
health_Luekemia_scaled <- t(health_Luekemia_scaled)
RF_health <- health_Luekemia_scaled[rownames(health_Luekemia_scaled) %in% healthy_sample,]
#filter real sig genes from healthy part
RF_health <- RF_health[,significant_gene_names_real]
dim(RF_health)

#load RF disease
disease_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/disease_Luekemia_scaled.RDS")
disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
RF_disease <- disease_Luekemia_scaled[rownames(disease_Luekemia_scaled) %in% disease_sample,]
dim(RF_disease)
#filter real sig genes from disease part
RF_disease <- RF_disease[,significant_gene_names_real]
dim(RF_disease)

# #loading real data itself (already is loaded)
# real_RNA_data <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/genes_expected_count_DESeq2_H3K27acFormatted.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
# dim(real_RNA_data)
# rownames(real_RNA_data) = real_RNA_data[,1]
# real_RNA_data = real_RNA_data[,-1]
# dim(real_RNA_data)
#picking filtered genes:
real_RNA_data_filtered <- real_RNA_data[significant_gene_names_real,]
#picking healthy samples
real_RNA_data_healthy <- real_RNA_data_filtered[,rownames(RF_health)]
real_RNA_data_healthy <- t(real_RNA_data_healthy)
health_lable <- paste0("Healthy", 1:10)
rownames(real_RNA_data_healthy) <- health_lable
#picking disease samples
real_RNA_data_disease <- real_RNA_data_filtered[, rownames(RF_disease)]
real_RNA_data_disease <- t(real_RNA_data_disease)
Disease_lable <- paste0("D", 1:12)
rownames(real_RNA_data_disease) <- Disease_lable
dim(real_RNA_data_disease)

library(ggplot2)
# Calculate the log2 fold change for RF data
RF_disease_means <- colMeans(RF_disease)
RF_disease_means[RF_disease_means == -1] <- 0
RF_health_means <- colMeans(RF_health)
RF_health_means[RF_health_means == -1] <- 0
log2_FC_RF <- log2((RF_health_means +1) / (RF_disease_means + 1))

# Calculate the log2 fold change for real RNA data
real_RNA_disease_means <- colMeans(real_RNA_data_disease)
real_RNA_health_means <- colMeans(real_RNA_data_healthy)
log2_FC_real <- log2((real_RNA_health_means +1) / (real_RNA_disease_means + 1))

# Create a data frame for plotting
genes <- colnames(RF_disease)
plot_data <- data.frame(
  Gene = genes,
  log2_FC_RF = log2_FC_RF,
  log2_FC_real = log2_FC_real
)
correlation <- cor(log2_FC_RF, log2_FC_real)
print(correlation)


# #--------------------------- plot part with grids in background ----------------------
# # Plot the scatter plot with matching axes limits
# ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_RF, label = Gene)) +
#   geom_point() +
#   # Set both x and y axis limits to -10 to 10
#   xlim(-10, 10) +
#   ylim(-10, 10) +
#   labs(
#     title = "Scatter plot of log2 Fold Change for ~15,000 sig gene from Real",
#     x = "log2 Fold Change (Real Expression )",
#     y = "log2 Fold Change (RF Expression)"
#   ) +
#   theme_minimal() 
# # + annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
# #       hjust = 1.1, vjust = -1.1, size = 4, color = "black")




#  ---------------------- no grid version of the plot  ----------------------

# Plot the scatter plot with axis limits and without background grid lines
ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_RF, label = Gene)) +
  geom_point() +
  # Set both x and y axis limits to -10 to 10
  xlim(-10, 10) +
  ylim(-10, 10) +
  labs(
    title = "Scatter plot of log2 Fold Change for ~15,000 sig gene from Real",
    x = "log2 Fold Change (Real Expression )",
    y = "log2 Fold Change (RF Expression)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),      # Remove border around the plot
    axis.line = element_line(color = "black")  # Add only x and y axis lines
  ) 
# + annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
#           hjust = 1.1, vjust = -1.1, size = 4, color = "black")
################################################################################################################################################################################################################################################

############################################## Scatter plot for real data vs CNN (on significant genes from real data) ##########################################################################################################################
kept_genes <- new_kept_genes
metadata_Luekemia <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/IHEC_metadata_harmonization.v1.1_BCellLeukemiaSamples.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
col = c("epirr_id_without_version", "EpiClass_pred_Sex", "EpiClass_pred_Life_stage", "harmonized_sample_disease_high", "RNA_available", "project", "harmonized_tissue_type")
metadata_Luekemia = metadata_Luekemia[,col]
colnames(metadata_Luekemia) = c("sample", "sex", "age", "health_status", "RNA_available", "project", "harmonized_tissue_type")
metadata_Luekemia$health_status <- gsub('Healthy/None', 'Healthy', metadata_Luekemia$health_status)
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$age == "adult",]
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$RNA_available == "True",]
dim(metadata_Luekemia)

outliers = c("IHECRE00000785",  "IHECRE00000234","IHECRE00000298", "IHECRE00000732")
samples = setdiff(samples, outliers)
metadata_Luekemia <- metadata_Luekemia[metadata_Luekemia$sample %in% samples,]

# samples <- setdiff(samples,"IHECRE00000298")
healthy_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Healthy",]
healthy_sample <- healthy_df$sample
disease_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Cancer",]
disease_sample <- disease_df$sample

#loading significant genes from real data
significant_gene_names_real= readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_real.txt")
length(significant_gene_names_real)

#load CNN healthy
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled_correct_scale_CNN.RDS")
#health_Luekemia_scaled <- t(health_Luekemia_scaled)
CNN_health <- health_Luekemia_scaled[rownames(health_Luekemia_scaled) %in% healthy_sample,]
CNN_health <- CNN_health[, kept_genes]
dim(CNN_health)
#filter sig genes:
CNN_health <- CNN_health[,significant_gene_names_real]
dim(CNN_health)


#load CNN disease
disease_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/disease_Luekemia_scaled_correct_scale_CNN.RDS")
#disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
CNN_disease <- disease_Luekemia_scaled[rownames(disease_Luekemia_scaled) %in% disease_sample,]
CNN_disease <- CNN_disease[, kept_genes]
#filter sig genes:
CNN_disease <- CNN_disease[,significant_gene_names_real]

# #loading real data itself (its already loaded thats why I cm it here)
# real_RNA_data <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/genes_expected_count_DESeq2_H3K27acFormatted.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
# rownames(real_RNA_data) = real_RNA_data[,1]
# real_RNA_data = real_RNA_data[,-1]

# filtered the sig genes:
real_RNA_data_filtered <- real_RNA_data[significant_gene_names_real,]
#picking healthy samples
real_RNA_data_healthy <- real_RNA_data_filtered[,rownames(RF_health)]
real_RNA_data_healthy <- t(real_RNA_data_healthy)
health_lable <- paste0("Healthy", 1:10)
rownames(real_RNA_data_healthy) <- health_lable
#picking disease samples
real_RNA_data_disease <- real_RNA_data_filtered[, rownames(RF_disease)]
real_RNA_data_disease <- t(real_RNA_data_disease)
Disease_lable <- paste0("D", 1:12)
rownames(real_RNA_data_disease) <- Disease_lable


# Calculate the log2 fold change for CNN data
CNN_disease_means <- colMeans(CNN_disease)
CNN_disease_means[CNN_disease_means == -1] <- 0
CNN_health_means <- colMeans(CNN_health)
CNN_health_means[CNN_health_means == -1] <- 0
log2_FC_CNN <- log2((CNN_health_means +1) / (CNN_disease_means + 1))

# Calculate the log2 fold change for real RNA data
real_RNA_disease_means <- colMeans(real_RNA_data_disease)
real_RNA_health_means <- colMeans(real_RNA_data_healthy)
log2_FC_real <- log2((real_RNA_health_means +1) / (real_RNA_disease_means + 1))


# Create a data frame for plotting
genes <- colnames(CNN_disease)
plot_data <- data.frame(
  Gene = genes,
  log2_FC_CNN = log2_FC_CNN,
  log2_FC_real = log2_FC_real
)
correlation <- cor(log2_FC_CNN, log2_FC_real)
print(correlation) #before removing outliers smaller than -10, greater than +10

library(ggplot2)
# #-------------------------- Plot the scatter plot ---------------------
# ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_CNN, label = Gene)) +
#   geom_point() +
#   # geom_abline(slope = 1, intercept = 0, color = "red")+
#   labs(
#     title = "Scatter plot of log2 Fold Change for ~15,000 real sig genes",
#     x = "log2 Fold Change (Real Expression )",
#     y = "log2 Fold Change (CNN Expression)"
#   ) +
#   theme_minimal() +
#   annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
#            hjust = 1.1, vjust = -1.1, size = 4, color = "black" )
# 
# 
# library(ggplot2)
# # Plot the scatter plot with axis limits
# ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_CNN, label = Gene)) +
#   geom_point() +
#   # Limit both x and y axes between -10 and 10
#   xlim(-10, 10) +
#   ylim(-10, 10) +
#   labs(
#     title = "Scatter plot of log2 Fold Change for ~15,000 real sig genes",
#     x = "log2 Fold Change (Real Expression )",
#     y = "log2 Fold Change (CNN Expression)"
#   ) +
#   theme_minimal() 
# # + annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
# #           hjust = 1.1, vjust = -1.1, size = 4, color = "black")


#----------------- plotting with no grid background --------------- 
library(ggplot2)
# Plot the scatter plot with axis limits and no grid lines
ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_CNN, label = Gene)) +
  geom_point() +
  # Limit both x and y axes between -10 and 10
  xlim(-10, 10) +
  ylim(-10, 10) +
  labs(
    title = "Scatter plot of log2 Fold Change for ~15,000 real sig genes",
    x = "log2 Fold Change (Real Expression )",
    y = "log2 Fold Change (CNN Expression)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),      # Remove the plot border
    axis.line = element_line(color = "black")  # Keep only x and y axis lines
  ) 
# + annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
#           hjust = 1.1, vjust = -1.1, size = 4, color = "black")

# Filter the data to include only points within the desired range
filtered_data <- plot_data[log2_FC_real >= -10 & log2_FC_real <= 10 & 
                             log2_FC_CNN >= -10 & log2_FC_CNN <= 10, ]

# Calculate the correlation for the filtered data
correlation_filtered <- cor(filtered_data$log2_FC_real, filtered_data$log2_FC_CNN)
print(correlation_filtered)
