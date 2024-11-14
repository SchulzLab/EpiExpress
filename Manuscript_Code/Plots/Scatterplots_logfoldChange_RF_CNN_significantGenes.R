################################################### scatterplot for log fold change for real vs RF ########################################################################

#libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("airway")
library(airway)
library(tidyverse)
library(ggplot2)

# shared genes between kept, CNN and RF
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled.RDS")
genes_from_RF <- rownames(health_Luekemia_scaled)
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/CNN/CNN_Luekemia_backscaled_healthy.RDS")
genes_from_CNN <- rownames(health_Luekemia_scaled)
kept_genes <- read.table("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/kept_genes.txt", header = FALSE, skip = 1)[, 2]


shared_genes1 <- intersect(genes_from_RF, genes_from_CNN)
shared_genes2 <- intersect(shared_genes1, kept_genes)
print(length(shared_genes2))
new_kept_genes <- shared_genes2

# kept_genes <- readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/shared_genes_CNN_RF_kept.txt")

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




health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/health_Luekemia_scaled.RDS")
health_Luekemia_scaled <- t(health_Luekemia_scaled)
RF_health <- health_Luekemia_scaled[rownames(health_Luekemia_scaled) %in% healthy_sample,]
dim(RF_health)


significant_gene_names_RF= readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
# significant_gene_names_RF= readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")

length(significant_gene_names_RF)
#filter genes:
RF_health <- RF_health[,significant_gene_names_RF]
dim(RF_health)

#load RF disease
disease_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/disease_Luekemia_scaled.RDS")
disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
RF_disease <- disease_Luekemia_scaled[rownames(disease_Luekemia_scaled) %in% disease_sample,]
dim(RF_disease)

#filter genes:
RF_disease <- RF_disease[,significant_gene_names_RF]
dim(RF_disease)

real_RNA_data <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/genes_expected_count_DESeq2_H3K27acFormatted.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
dim(real_RNA_data)
rownames(real_RNA_data) = real_RNA_data[,1]
real_RNA_data = real_RNA_data[,-1]
dim(real_RNA_data)
#picking filtered genes:
real_RNA_data_filtered <- real_RNA_data[significant_gene_names_RF,]
#picking healthy Luekemia
real_RNA_data_healthy <- real_RNA_data_filtered[,rownames(RF_health)]
real_RNA_data_healthy <- t(real_RNA_data_healthy)
health_lable <- paste0("Healthy", 1:10)
rownames(real_RNA_data_healthy) <- health_lable
#picking disease Luekemia
real_RNA_data_disease <- real_RNA_data_filtered[, rownames(RF_disease)]
real_RNA_data_disease <- t(real_RNA_data_disease)
Disease_lable <- paste0("D", 1:12)
rownames(real_RNA_data_disease) <- Disease_lable
dim(real_RNA_data_disease)

library(ggplot2)
# Calculate the log2 fold change for CNN data
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



# Plot the scatter plot with matching axes limits
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
  annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")




#no grid version

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
  ) +
  annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")



#####CNN part
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

#load CNN healthy
health_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/CNN/CNN_Luekemia_backscaled_healthy.RDS")
health_Luekemia_scaled <- t(health_Luekemia_scaled)
CNN_health <- health_Luekemia_scaled[rownames(health_Luekemia_scaled) %in% healthy_sample,]
CNN_health <- CNN_health[, kept_genes]
dim(CNN_health)

significant_gene_names_CNN = readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_CNN.txt")
# significant_gene_names_CNN = readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")

#filter genes:
CNN_health <- CNN_health[,significant_gene_names_CNN]
dim(CNN_health)

#load CNN disease
disease_Luekemia_scaled <- readRDS("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/CNN/CNN_Luekemia_backscaled_disease.RDS")
disease_Luekemia_scaled <- t(disease_Luekemia_scaled)
CNN_disease <- disease_Luekemia_scaled[rownames(disease_Luekemia_scaled) %in% disease_sample,]
CNN_disease <- CNN_disease[, kept_genes]


#filter genes:
CNN_disease <- CNN_disease[,significant_gene_names_CNN]


real_RNA_data <- read.csv("~/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/genes_expected_count_DESeq2_H3K27acFormatted.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
rownames(real_RNA_data) = real_RNA_data[,1]
real_RNA_data = real_RNA_data[,-1]

#picking filtered genes:
real_RNA_data_filtered <- real_RNA_data[significant_gene_names_CNN,]
#picking healthy Luekemia
real_RNA_data_healthy <- real_RNA_data_filtered[,rownames(RF_health)]
real_RNA_data_healthy <- t(real_RNA_data_healthy)
health_lable <- paste0("Healthy", 1:10)
rownames(real_RNA_data_healthy) <- health_lable
#picking disease Luekemia
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
print(correlation)
library(ggplot2)
# Plot the scatter plot
ggplot(plot_data, aes(x = log2_FC_real, y = log2_FC_CNN, label = Gene)) +
  geom_point() +
  # geom_abline(slope = 1, intercept = 0, color = "red")+
  labs(
    title = "Scatter plot of log2 Fold Change for ~15,000 real sig genes",
    x = "log2 Fold Change (Real Expression )",
    y = "log2 Fold Change (CNN Expression)"
  ) +
  theme_minimal() +
  annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black" )


library(ggplot2)
# Plot the scatter plot with axis limits
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
  annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")


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
  ) +
  annotate("text", x = Inf, y = -Inf, label = paste0("Cor:", round(correlation, 3)),
           hjust = 1.1, vjust = -1.1, size = 4, color = "black")

