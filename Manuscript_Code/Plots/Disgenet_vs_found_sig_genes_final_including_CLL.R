
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
CNN_genes <- readLines("path/to/data/Sig_CNN_genes.txt")
RF_genes <- readLines("path/to/data/Sig_RF_genes.txt")
real_genes <- readLines("path/to/data/Sig_Real_genes.txt")
all_genes <- readLines("path/to/data/sig_all.txt")

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
  
  return(fisher_result$p.value)
}

# Count genes per disease
disease_gene_counts <- AML_datasets %>%
  dplyr::filter(Disease != "Childhood Acute Myeloid Leukemia") %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarize(gene_count = dplyr::n()) %>%
  dplyr::arrange(desc(gene_count))


top_diseases <- c("CLL" ,"Acute Myeloid Leukemia (AML-M2)", "Cytogenetically normal acute myeloid leukemia", 
                  "Treatment related acute myeloid leukaemia", "Adult Acute Myeloblastic Leukemia")
# Filter AML_datasets to only include the top 8 diseases
AML_top_datasets <- AML_datasets %>% filter(Disease %in% top_diseases)

# Loop through each disease to calculate enrichment
for (disease in top_diseases) {
  disease_genes <- AML_top_datasets %>% filter(Disease == disease) %>% pull(Gene)
  
  # Calculate raw p-values for each gene list
  real_pval <- calc_enrichment(real_genes, disease_genes, all_genes)
  RF_pval <- calc_enrichment(RF_genes, disease_genes, all_genes)
  CNN_pval <- calc_enrichment(CNN_genes, disease_genes, all_genes)
  
  #clean
  enrichment_data <- rbind(enrichment_data, data.frame(
    Disease = disease,  
    real_p = real_pval,
    RF_p = RF_pval,
    CNN_p = CNN_pval
  ))
  
}

# Reshape the enrichment_data for plotting
enrichment_data_melt <- enrichment_data %>%
  pivot_longer(cols = c(
    "real_p", "RF_p", "CNN_p"),  
               names_to = "List", values_to = "Raw_P")%>%
mutate(List = case_match(List, 
                           "real_p" ~ "real_adj", 
                           "RF_p"   ~ "RF_adj", 
                           "CNN_p"  ~ "CNN_adj"))

#global Adj P-value:
enrichment_data_melt <- enrichment_data_melt %>%
  mutate(P_adj = p.adjust(Raw_P, method = "BH"),
         Enrichment = -log10(P_adj + 1e-10))

# Specify the desired order of the List factor
enrichment_data_melt$List <- factor(enrichment_data_melt$List,  
                                    levels = c(
                                      "real_adj", "RF_adj", "CNN_adj"))



#order
enrichment_data_melt$List <- factor(enrichment_data_melt$List, 
                                    levels = c("real_adj", "RF_adj", "CNN_adj"))

# get near-zero enrichment values to be able to see all diasese
enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1e-10)
# get enrichment values greater than 1.3 which is adj p-value 0.05
#enrichment_data_melt <- enrichment_data_melt %>% filter(Enrichment > 1.3)
A = as.data.frame(enrichment_data_melt)

#average
A$Disease <- factor(A$Disease, 
                    levels = A %>%
                      group_by(Disease) %>%
                      mutate(AvgEnrichment = mean(Enrichment, na.rm = TRUE)) %>%  # Use mean instead of max
                      distinct(Disease, AvgEnrichment) %>%
                      arrange(desc(AvgEnrichment)) %>%
                      dplyr::pull(Disease) %>%
                      rev())  # Reverse the order



# Plot
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

