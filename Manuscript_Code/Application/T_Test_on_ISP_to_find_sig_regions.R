print("T-test on all regions from all genes")

directory <- "/path/to/data/" #input file(genes with ISP)

# Get a list of files in the directory
files <- list.files(directory,pattern = "\\.RDS$", full.names = TRUE)
print(length(files))


#finding healthy and disease samples for adults 
metadata_Luekemia <- read.csv("/path/to/data//IHEC_metadata_harmonization.v1.1_BCellLeukemiaSamples.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE)
print(dim(metadata_Luekemia))
# 


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

metadata_Luekemia <- metadata_Luekemia[ metadata_Luekemia$sample %in% samples,]
healthy_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Healthy",]
healthy_sample <- healthy_df$sample
disease_df <- metadata_Luekemia[metadata_Luekemia$health_status == "Cancer",]
disease_sample <- disease_df$sample

print(length(healthy_sample))
print(length(disease_sample))

# Iterate over each file
for (gene_name in files) {
  print(gene_name)
  SHAP_Sig_ENSG <- readRDS(gene_name)
  # Extract base name of the gene file
  gene_name <- basename(gene_name)
  gene_name <- sub("\\.RDS$", "", gene_name)  # Remove file extension
  
  # 2)  T.test ---------------------------------------------------------------------------
  
  df_shap = as.data.frame(SHAP_Sig_ENSG)
  
  # Split SHAP values based on labels
  shap_healthy <- df_shap[healthy_sample, ]
  shap_disease <- df_shap[disease_sample, ]
  A = rbind(shap_healthy, shap_disease)
  # Split SHAP values based on labels
  shap_healthy <- A[1:10, ]
  shap_disease <- A[11:22, ]
  print(dim(shap_healthy))
  print(dim(shap_disease))
  # Perform t-test for each feature
  t_test_results <- sapply(names(A), function(feature_name) {
    t_test_result <- t.test(shap_healthy[[feature_name]], shap_disease[[feature_name]], paired = FALSE)
    return(t_test_result$p.value)
  })
  
  T_test_results <- data.frame(names = names(t_test_results), values = as.numeric(t_test_results))

  #save all result regardless of p-value amount
  write.table(T_test_results, 
              file = paste0("provide/a/path/to/save/the/result/T_Test_all_",gene_name, ".txt"),
              sep = "\t", row.names = FALSE,  quote = FALSE)
  
}
print("end of the code")
