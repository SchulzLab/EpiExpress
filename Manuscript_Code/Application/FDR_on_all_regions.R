

###################################### I pick all regions without filtering and give them to the FDR  ##################

# Load required libraries
library(dplyr)
library(readr)
library(stats)


directory <- "path/to/T_Test_all_with_gene_name2"
# Get the list of all .txt files in the directory
file_list <- list.files(path = directory, pattern = "*.txt", full.names = TRUE)
print(length(file_list))


# Define a function to read files with specified column types
read_file_with_types <- function(file) {
  read_tsv(file, col_types = cols(
    names = col_character(),
    values = col_double(),
    Gene = col_character() 
  ), show_col_types = FALSE)
}

# Read and concatenate all files
all_data <- file_list %>%
  lapply(read_file_with_types) %>%
  bind_rows()

# Count the number of lines before picking unique regions
num_lines_before <- nrow(all_data)
print(paste("Number of lines before picking unique regions:", num_lines_before))



# Perform FDR (Benjamini-Hochberg) correction
all_data <- all_data %>%
  mutate(adj_p_value = p.adjust(values, method = "BH")) %>%
  select(names, adj_p_value, Gene)  # Select only the required columns
print("______________________________________________________")
print(all_data[1:4,])
print(dim(all_data))

# Save the full results with FDR adjustment to a new file
# full_output_file <- "/path/to//FDR_all_on_all_regions_Luekemia_RF_no_filtering_before_FDR_W2_with_gene_names.txt"
# full_output_file <- "/path/to//Lima_FDR_all_on_all_regions_Luekemia_RF_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt"
full_output_file <- "path/to/FDR_all_on_all_regions_Luekemia_CNN_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt"
# Apply the function to the 'names' column
# all_data$names <- sapply(all_data$names, transform_values)

write_tsv(all_data, full_output_file)

# here I need to pick unique regions with smallest p-value
unique_data_after_FDR <- all_data %>%
  group_by(names) %>%
  filter(adj_p_value == min(adj_p_value)) %>%
  ungroup()

# Filter the data to keep only entries with adjusted p-values < 0.05
filtered_data <- unique_data_after_FDR %>%
  filter(adj_p_value < 0.05)
print("______________________________________________________")

# Save the filtered data to a new file
# filtered_output_file <- "/path/to//sig_unique_regions_after_FDR_smaller_0.05_Luekemia_RF_W2_with_gene_names.txt"
filtered_output_file <- "path/to/sig_unique_regions_after_FDR_smaller_0.05_Luekemia_CNN_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt"
# 
# filtered_data$names <- sapply(filtered_data$names, transform_values)
write_tsv(filtered_data, filtered_output_file)


filtered_output_file2 <- "path/to/JUST_sig_unique_regions_after_FDR_smaller_0.05_Luekemia_CNN_gens_byDSEQ_COV_SEX_no_filtering_before_FDR.txt"

writeLines(filtered_data$names, filtered_output_file2)
print("done")




# regions_modified <- gsub("\\.", "-", a)
