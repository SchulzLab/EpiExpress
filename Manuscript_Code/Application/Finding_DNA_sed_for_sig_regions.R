# Load necessary libraries
library(data.table)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

# Read the data from the file (including the Gene column) #enter either CNN or RF
# data <- fread("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_unique_regions_after_FDR_smaller_0.05_Luekemia_RF_W2_with_gene_names_ranked.txt")
data <- fread("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/Luekemia_Apog_Aplication/sig_unique_regions_after_FDR_smaller_0.05_Luekemia_RF_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt")

# Filter the data for p-values smaller than 0.05
data <- data %>%
  filter(adj_p_value < 0.05)

# Function to transform the 'names' values to match genome coordinates
transform_values <- function(value) {
  value <- gsub("^X(\\d+)", "chr\\1", value) # Replace 'X' followed by digits with 'chr'
  value <- gsub("^X\\.", "chrX.", value)     # Replace 'X.' with 'chrX.'
  value <- gsub("^Y\\.", "chrY.", value)     # Replace 'Y.' with 'chrY.'
  return(value)
}

# Apply the function to the 'names' column
data$names <- sapply(data$names, transform_values)
# data = a
# Split the 'names' column to extract chromosome, start, and end positions
regions <- strsplit(data$names, "[.]")
regions <- lapply(regions, function(x) {
  c(chromosome = x[1], start = as.numeric(x[2]), end = as.numeric(x[3]))
})

# Convert to data frame
regions_df <- as.data.frame(do.call(rbind, regions), stringsAsFactors = FALSE)

# Convert to GRanges object
gr <- GRanges(seqnames = regions_df$chromosome,
              ranges = IRanges(start = as.numeric(regions_df$start), end = as.numeric(regions_df$end)))

# Write to BED file
export(gr, "regions.bed")

#or 

library(rtracklayer)
export.bed(gr,con='regions1.bed')

# Import the regions from the BED file (to ensure the correct format)
regions <- import("regions1.bed")

# Load the genome sequence
genome <- readDNAStringSet("/Users/shamim/Desktop/PhD/ML_project/ml_scripts/IHEC_Project/All_IHEC_result/1MB/hg38.fa")

# Extract sequences based on the regions
sequences <- getSeq(genome, regions)
sig_regions_seq <- sequences

# Add Gene names as headers to the significant region sequences
names(sig_regions_seq) <- data$Gene

# Save significant region sequences to a FASTA file with Gene names as identifiers
writeXStringSet(sig_regions_seq, filepath="provid/path/to/save/DNA_seq_for_sig_regions.fasta")
