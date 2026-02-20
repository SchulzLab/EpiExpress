
#-----------code with ordering disease and picking 5 category with higest enrichment (average among three methods)
# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readxl)
library(jsonlite)

# Define the file path
file_path <- "path/to/data/CLL_genes_Harmonizome.txt"  # 

# Read the file as a string
json_text <- paste(readLines(file_path), collapse = "")

# Parse JSON content
json_data <- fromJSON(json_text)

# Extract gene associations
gene_data <- json_data$associations
gene_names <- gene_data$gene
CLL_gene_names <- gene_names$symbol
CLL_gene_names <- data.frame(Disease = "CLL", Gene = CLL_gene_names)

# Read the Excel file
df <- read_excel("path/to/data/All_Luekemia_gene_names.xlsx")
df <- df[-1,]
AML_datasets1 <- df[,1:2]
colnames(AML_datasets1) <- c("Disease", "Gene")
AML_datasets1 <- as.data.frame(AML_datasets1)

# Combine CLL and AML datasets
AML_datasets <- rbind(CLL_gene_names, AML_datasets1)

# Read gene lists
CNN_genes <- readLines("path/to/data//Symbol_sig_new_CNN_genes.txt")
RF_genes <- readLines("path/to/data//SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
real_genes <- readLines("path/to/data//SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")
all_genes <- readLines("path/to/data/SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_ALL.txt")

# Calculate the unique sets for CNN and RF
common_CNN_RF <- intersect(RF_genes, CNN_genes)
unique_RF <- setdiff(RF_genes, union(real_genes, CNN_genes))
unique_CNN <- setdiff(CNN_genes, union(real_genes, RF_genes))

# Initialize enrichment_data to store results
enrichment_data <- data.frame(Disease = character(), 
                              # unique_RF_adj = numeric(),  
                              # unique_CNN_adj = numeric(),  
                              real_adj = numeric(),
                              RF_adj = numeric(),
                              CNN_adj = numeric(),
                              stringsAsFactors = FALSE)

# Define a function to calculate enrichment using Fisher's test
calc_enrichment <- function(gene_list, disease_genes, all_genes) {
  NotSigGenes <- setdiff(all_genes, gene_list)
  in_Disease_and_SigGenes <- length(intersect(disease_genes, gene_list))
  in_Disease_and_NotSigGenes <- length(intersect(disease_genes, NotSigGenes))
  in_NotDisease_and_SigGenes <- length(intersect(setdiff(all_genes, disease_genes), gene_list))
  in_NotDisease_and_NotSigGenes <- length(intersect(setdiff(all_genes, disease_genes), NotSigGenes))
  
  fisher_result <- fisher.test(
    matrix(c(in_Disease_and_SigGenes, in_Disease_and_NotSigGenes,
             in_NotDisease_and_SigGenes, in_NotDisease_and_NotSigGenes),
           nrow = 2, byrow = TRUE)
  )
  
  return(p.adjust(fisher_result$p.value, method = "BH"))
}

# Count genes per disease
disease_gene_counts <- AML_datasets %>%
  dplyr::filter(Disease != "Childhood Acute Myeloid Leukemia") %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarize(gene_count = dplyr::n()) %>%
  dplyr::arrange(desc(gene_count))

# # Select the top 8 diseases with the highest number of genes
# top_diseases <- disease_gene_counts %>%
#   slice(1:8) %>%
#   pull(Disease)

top_diseases <- c("CLL" ,"Acute Myeloid Leukemia (AML-M2)", "Cytogenetically normal acute myeloid leukemia", 
                  "Treatment related acute myeloid leukaemia", "Adult Acute Myeloblastic Leukemia")
# Filter AML_datasets to only include the top 8 diseases
AML_top_datasets <- AML_datasets %>% filter(Disease %in% top_diseases)

# Loop through each disease to calculate enrichment
for (disease in top_diseases) {
  disease_genes <- AML_top_datasets %>% filter(Disease == disease) %>% pull(Gene)
  
  # Calculate adjusted p-values for each gene list
  real_adj_pval <- calc_enrichment(real_genes, disease_genes, all_genes)
  # unique_RF_adj_pval <- calc_enrichment(unique_RF, disease_genes, all_genes)
  # unique_CNN_adj_pval <- calc_enrichment(unique_CNN, disease_genes, all_genes)
  RF_p_adj_value <- calc_enrichment(RF_genes, disease_genes, all_genes)
  CNN_p_adj_value <- calc_enrichment(CNN_genes, disease_genes, all_genes)
  
  # Store the -log10(adjusted p-values) in the enrichment_data dataframe
  enrichment_data <- rbind(enrichment_data, data.frame(
    Disease = disease,  
    # unique_RF_adj = -log10(unique_RF_adj_pval + 1e-10),
    # unique_CNN_adj = -log10(unique_CNN_adj_pval + 1e-10),
    real_adj = -log10(real_adj_pval + 1e-10),
    RF_adj = -log10(RF_p_adj_value + 1e-10),
    CNN_adj = -log10(CNN_p_adj_value + 1e-10)
  ))
}

# Reshape the enrichment_data for plotting
enrichment_data_melt <- enrichment_data %>%
  pivot_longer(cols = c(
   # "unique_RF_adj", "unique_CNN_adj",
    "real_adj", "RF_adj", "CNN_adj"),  
               names_to = "List", values_to = "Enrichment")

# Specify the desired order of the List factor
enrichment_data_melt$List <- factor(enrichment_data_melt$List,  
                                    levels = c(
                                      # "unique_RF_adj", "unique_CNN_adj", 
                                      "real_adj", "RF_adj", "CNN_adj"))

# Filter out rows with near-zero enrichment values
enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1e-10)

A = as.data.frame(enrichment_data_melt)
# #max
# A$Disease <- factor(A$Disease, 
#                     levels = A %>%
#                       group_by(Disease) %>%
#                       mutate(MaxEnrichment = max(Enrichment, na.rm = TRUE)) %>%
#                       distinct(Disease, MaxEnrichment) %>%
#                       arrange(desc(MaxEnrichment)) %>%
#                       dplyr::pull(Disease) %>%
#                       rev())  # Reverse the order

#average
A$Disease <- factor(A$Disease, 
                    levels = A %>%
                      group_by(Disease) %>%
                      mutate(AvgEnrichment = mean(Enrichment, na.rm = TRUE)) %>%  # Use mean instead of max
                      distinct(Disease, AvgEnrichment) %>%
                      arrange(desc(AvgEnrichment)) %>%
                      dplyr::pull(Disease) %>%
                      rev())  # Reverse the order



# Now plot again
ggplot(A, aes(x = List, y = Disease, size = Enrichment, color = Enrichment)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +  # Adjust dot sizes
  scale_color_gradient(low = "lightgreen", high = "darkgreen") +  # Gradient color
  theme_minimal() +
  labs(title = "Enrichment of Gene Lists (adj)",
       x = "Gene List",
       y = "Disease Dataset",
       size = "-log10(Adjusted adj_p-value)",
       color = "-log10(Adjusted adj_p-value)") +
  theme(axis.text.y = element_text(size = 8),  # Adjust text size
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))  # Rotate x-axis labels





# #----same without sorting disease categories based on enrichment
# #libraries
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(data.table)
# library(ggplot2)
# library(readxl)
# library(jsonlite)
# 
# 
# # Define the file path
# file_path <- "path/to/data/CLL_genes_Harmonizome.txt" 
# 
# # Read the file as a string
# json_text <- paste(readLines(file_path), collapse = "")
# 
# # Parse JSON content
# json_data <- fromJSON(json_text)
# 
# # Extract gene associations
# gene_data <- json_data$associations
# gene_names <- gene_data$gene
# CLL_gene_names <- gene_names$symbol
# CLL_gene_names <- data.frame(Disease = "CLL", Gene = CLL_gene_names)
# 
# df <- read_excel("path/to/data/All_Luekemia_gene_names.xlsx")
# df <- df[-1,]
# AML_datasets1 <- df[,1:2]
# colnames(AML_datasets1) <- c("Disease", "Gene")
# AML_datasets1 <- as.data.frame(AML_datasets1)
# 
# AML_datasets = rbind(CLL_gene_names, AML_datasets1)
# 
# CNN_genes <- readLines("path/to/data/New_running_after_correctin_CNN/Symbol_sig_new_CNN_genes.txt")
# length(CNN_genes)
# RF_genes <- readLines("path/to/data/SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
# length(RF_genes)
# real_genes <- readLines("path/to/data/SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")
# length(real_genes)
# all_genes <-  readLines("path/to/data/Luekemia_Apog_Aplication/SYMBOL_sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_ALL.txt")
# 
# # Calculate the unique sets for CNN and RF
# common_CNN_RF <- intersect(RF_genes, CNN_genes)
# unique_RF <- setdiff(RF_genes, union(real_genes, CNN_genes))
# unique_CNN <- setdiff(CNN_genes, union(real_genes, RF_genes))
# 
# 
# # Initialize enrichment_data to store results
# enrichment_data <- data.frame(Disease = character(), 
#                               unique_RF_adj = numeric(),  
#                               unique_CNN_adj = numeric(),  
#                               real_adj = numeric(),
#                               RF_adj = numeric(),
#                               CNN_adj = numeric(),
#                               stringsAsFactors = FALSE)
# 
# # Define a function to calculate enrichment using Fisher's test
# calc_enrichment <- function(gene_list, disease_genes, all_genes) {
#   NotSigGenes <- setdiff(all_genes, gene_list)
#   in_Disease_and_SigGenes <- length(intersect(disease_genes, gene_list))
#   in_Disease_and_NotSigGenes <- length(intersect(disease_genes, NotSigGenes))
#   in_NotDisease_and_SigGenes <- length(intersect(setdiff(all_genes, disease_genes), gene_list))
#   in_NotDisease_and_NotSigGenes <- length(intersect(setdiff(all_genes, disease_genes), NotSigGenes))
#   
#   fisher_result <- fisher.test(
#     matrix(c(in_Disease_and_SigGenes, in_Disease_and_NotSigGenes,
#              in_NotDisease_and_SigGenes, in_NotDisease_and_NotSigGenes),
#            nrow = 2, byrow = TRUE)
#   )
#   
#   return(p.adjust(fisher_result$p.value, method = "BH"))
#   
# }
# 
# 
# disease_gene_counts <- AML_datasets %>%
#   dplyr::filter(Disease != "Childhood Acute Myeloid Leukemia") %>%
#   dplyr::group_by(Disease) %>%
#   dplyr::summarize(gene_count = dplyr::n()) %>%
#   dplyr::arrange(desc(gene_count))
# 
# 
# # # Select the top 5 diseases with the highest number of genes
# # top_diseases <- disease_gene_counts %>% top_n(5, gene_count) %>% pull(Disease)
# 
# top_diseases <- disease_gene_counts %>%
#   slice(1:8) %>%
#   pull(Disease)
# 
# # Filter AML_datasets to only include the top 5 diseases
# AML_top_datasets <- AML_datasets %>% filter(Disease %in% top_diseases)
# 
# 
# 
# # Loop through each disease to calculate enrichment
# for (disease in top_diseases) {
#   disease_genes <- AML_top_datasets %>% filter(Disease == disease) %>% pull(Gene)
#   
#   # Calculate adjusted p-values for each gene list
#   real_adj_pval <- calc_enrichment(real_genes, disease_genes, all_genes)
#   # unique_RF_adj_pval <- calc_enrichment(unique_RF, disease_genes, all_genes)
#   # unique_CNN_adj_pval <- calc_enrichment(unique_CNN, disease_genes, all_genes)
#   RF_p_adj_value <- calc_enrichment(RF_genes, disease_genes, all_genes)
#   CNN_p_adj_value <- calc_enrichment(CNN_genes, disease_genes, all_genes)
#   
#   
#   
#   # Store the -log10(adjusted p-values) in the enrichment_data dataframe
#   enrichment_data <- rbind(enrichment_data, data.frame(
#     Disease = disease,  
#     # unique_RF_adj = -log10(unique_RF_adj_pval + 1e-10),
#     # unique_CNN_adj = -log10(unique_CNN_adj_pval + 1e-10),
#     real_adj = -log10(real_adj_pval + 1e-10),
#     RF_adj = -log10(RF_p_adj_value + 1e-10),
#     CNN_adj = -log10(CNN_p_adj_value + 1e-10)
#   ))
# }
# 
# # Reshape the enrichment_data for plotting
# enrichment_data_melt <- enrichment_data %>%
#   pivot_longer(cols = c(
#     # "unique_RF_adj", "unique_CNN_adj",
#     "real_adj", "RF_adj", "CNN_adj"),  
#     names_to = "List", values_to = "Enrichment")
# 
# # Specify the desired order of the List factor
# enrichment_data_melt$List <- factor(enrichment_data_melt$List,  
#                                     levels = c(
#                                       # "unique_RF_adj", "unique_CNN_adj", 
#                                       "real_adj", "RF_adj", "CNN_adj"))
# 
# 
# # Filter out rows with near-zero enrichment values
# enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1e-10)
# 
# 
# 
# # Create the dot plot with gradient color
# ggplot(enrichment_data_melt, aes(x = List, y = Disease, size = Enrichment, color = Enrichment)) +
#   geom_point() +
#   scale_size_continuous(range = c(1, 10)) +  # Adjust the range for dot sizes
#   scale_color_gradient(low = "lightgreen", high = "darkgreen") +  # Apply gradient shades of green
#   theme_minimal() +
#   labs(title = "Enrichment of Gene Lists (adj)",
#        x = "Gene List",
#        y = "Disease Dataset",
#        size = "-log10(Adjusted adj_p-value)",  # Adjust label to reflect that these are adjusted
#        color = "-log10(Adjusted adj_p-value)") +  # Include the color legend
#   theme(axis.text.y = element_text(size = 8),  # Adjust text size for better readability
#         axis.text.x = element_text(size = 10, angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
# 
# 
# 
