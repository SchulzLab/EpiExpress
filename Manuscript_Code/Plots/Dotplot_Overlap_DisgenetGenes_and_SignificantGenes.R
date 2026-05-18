############# Dotplot for comaring overlap between Disgenet leukemia genes and genes predicted by our models ##################
#we downloaded Leukemia related datasets from Disgenet websit
#libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readxl)

df <- read_excel("path/to/data/All_Luekemia_gene_names.xlsx")
df <- df[-1,]
AML_datasets <- df[,1:2]
colnames(AML_datasets) <- c("Disease", "Gene")
AML_datasets <- as.data.frame(AML_datasets)

CNN_genes <- readLines("path/to/data/SYMBOL_Sig_CNN_genes.txt")
length(CNN_genes)
RF_genes <- readLines("path/to/data/SYMBOL_Sig_RF_genes.txt")
length(RF_genes)
real_genes <- readLines("path/to/data/SYMBOL_Sig_real_genes.txt")
length(real_genes)
all_genes <-  readLines("path/to/data/SYMBOL_sig_All_genes.txt")

# Calculate the unique sets for CNN and RF
common_CNN_RF <- intersect(RF_genes, CNN_genes)
unique_RF <- setdiff(RF_genes, union(real_genes, CNN_genes))
unique_CNN <- setdiff(CNN_genes, union(real_genes, RF_genes))



# Initialize enrichment_data to store results
enrichment_data <- data.frame(Disease = character(), 
                              unique_RF_adj = numeric(),  
                              unique_CNN_adj = numeric(),  
                              real_adj = numeric(),
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
  
  return(fisher_result$p.value)
  
}

# Count the number of genes for each disease and exclude "Childhood Acute Myeloid Leukemia"
disease_gene_counts <- AML_datasets %>%
  filter(Disease != "Childhood Acute Myeloid Leukemia") %>%
  group_by(Disease) %>%
  summarize(gene_count = n()) %>%
  arrange(desc(gene_count))

# Select the top 5 diseases with the highest number of genes
top_diseases <- disease_gene_counts %>% top_n(5, gene_count) %>% pull(Disease)

# Filter AML_datasets to only include the top 5 diseases
AML_top_datasets <- AML_datasets %>% filter(Disease %in% top_diseases)

# Loop through each disease to calculate enrichment
for (disease in top_diseases) {
  disease_genes <- AML_top_datasets %>% filter(Disease == disease) %>% pull(Gene)
  
  #calculate p-value
  real_pval <- calc_enrichment(real_genes, disease_genes, all_genes)
  unique_RF_pval <- calc_enrichment(unique_RF, disease_genes, all_genes)
  unique_CNN_pval <- calc_enrichment(unique_CNN, disease_genes, all_genes)
  
 
  enrichment_data <- rbind(enrichment_data, data.frame(
    Disease = disease,  
    unique_RF = unique_RF_pval,
    unique_CNN= unique_CNN_pval,
    real= real_pval 
  ))
}
# Reshape, adjust globally across all tests, and then apply -log10
  enrichment_data_melt <- enrichment_data %>%
    pivot_longer(cols = c("unique_RF", "unique_CNN", "real"),  
                 names_to = "List", values_to = "Raw_P") %>%
    mutate(P_adj = p.adjust(Raw_P, method = "BH"),             # BH adjustment
           Enrichment = -log10(P_adj + 1e-10),                 # math conversion
           List = case_match(List, 
                             "unique_RF"  ~ "unique_RF_adj", 
                             "unique_CNN" ~ "unique_CNN_adj", 
                             "real"       ~ "real_adj"))
  
  # Specify the factor order using the newly mapped names
  enrichment_data_melt$List <- factor(enrichment_data_melt$List,       
                                      levels = c("unique_RF_adj", "unique_CNN_adj", "real_adj"))


# get near-zero enrichment values to be able to see all diasese
 enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1e-10)
# get enrichment values greater than 1.3 which is adj p-value 0.05
#enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1.3)

# Plot
ggplot(enrichment_data_melt, aes(x = List, y = Disease, size = Enrichment, color = Enrichment)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +  # Adjust the range for dot sizes
  scale_color_gradient(low = "lightgreen", high = "darkgreen") +  # Apply gradient shades of green
  theme_minimal() +
  labs(title = "Enrichment of Gene Lists (adj)",
       x = "Gene List",
       y = "Disease Dataset",
       size = "-log10(Adjusted adj_p-value)",  # Adjust label to reflect that these are adjusted
       color = "-log10(Adjusted adj_p-value)") +  # Include the color legend
  theme(axis.text.y = element_text(size = 8),  # Adjust text size for better readability
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

