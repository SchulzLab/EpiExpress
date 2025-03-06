#converting ID Gene to Symbol by gtf file

Gene_ID <- readLines("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/sig_genes_CNN_corrected_scaled.txt")
# Replace 'annotation.gtf.gz' with the path to your GTF file
annotation <- "/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/gencode.v38.annotation.gtf"

# Read the GTF file line by line
lines <- readLines(gzfile(annotation))

# Initialize an empty mapping
gene_name_map <- vector("character", length = length(Gene_ID))
names(gene_name_map) <- Gene_ID

# Iterate through each line
for (line in lines) {
  # Skip comment lines and lines that are not about genes
  if (!startsWith(line, "#") && strsplit(line, "\t")[[1]][3] == "gene") {
    # Extract gene ID and name
    gene_id <- gsub('.*gene_id "([^"]+)".*', '\\1', line)
    gene_name <- gsub('.*gene_name "([^"]+)".*', '\\1', line)
    
    # Extract only the gene ID before the dot (if present)
    gene_id <- gsub("\\..*", "", gene_id)
    
    # Assign gene name to the mapping if gene ID exists in Gene_ID
    if (gene_id %in% Gene_ID) {
      gene_name_map[gene_id] <- gene_name
    }
  }
}



# Print the mapping
print(gene_name_map)

gene_symbol_df <- data.frame(Gene_ID = names(gene_name_map), Symbol = gene_name_map, stringsAsFactors = FALSE)

sum(is.na(gene_symbol_df$Symbol))
length(gene_symbol_df$Symbol)

writeLines(gene_symbol_df$Symbol,"/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/New_running_after_correctin_CNN/Symbol_sig_new_CNN_genes.txt")


