install.packages("VennDiagram")
library(VennDiagram)

# union_CNN_real = union(CNN_genes, real_genes)
# RF_extra = setdiff(RF_genes,union_CNN_real )
# common <- Reduce(intersect, list(RF_genes,real_genes, CNN_genes))

# Load the gene lists from the text files
CNN_genes <- readLines("path/to/data/New_running_after_correctin_CNN/sig_genes_CNN.txt")
RF_genes <- readLines("path/to/data/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
real_genes <- readLines("path/to/data/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(RF = RF_genes, CNN = CNN_genes, real = real_genes),
  category.names = c("RF", "CNN", "real"),
  filename = NULL,
  col = c("#0a63a5", "#ab0303", "black"),  # Circle colors
  fill = NA,  # No fill for circles
  alpha = 0,  # Fully transparent circles
  cex = 0.8,  # Smaller numbers inside the diagram
  cat.cex = 1.2,  # Category label size
  cat.col = c("#0a63a5", "#ab0303", "black"),  # Category label colors
  main = "Venn Diagram of Significant Gene Sets"
)

# Plot the Venn diagram
grid.draw(venn.plot)




# #same as upper one, just no text lable on the plot
# library(VennDiagram)
# 
# # Load the gene lists from the text files
# CNN_genes <- readLines("path/to/data/New_running_after_correctin_CNN/sig_genes_CNN.txt")
# RF_genes <- readLines("path/to/data/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_RF.txt")
# real_genes <- readLines("path/to/data/sig_genes_deseq_step_by_step_NO_factor_size_adult_luekemia_SEX_covariate_NO3_outleires_removed_UNKNOWN_sample_real.txt")
# 
# # Create the Venn diagram
# venn.plot <- venn.diagram(
#   x = list(RF = RF_genes, CNN = CNN_genes, real = real_genes),
#   filename = NULL,
#   col = c("#0a63a5", "#ab0303", "black"),  # Circle colors
#   fill = NA,  # No fill for circles
#   alpha = 0,  # Fully transparent circles
#   cex = 0,  # Hide numbers inside the diagram
#   cat.cex = 0,  # Hide category label size
#   show.category = FALSE,  # Do not show category names
#   main = "Venn Diagram of Significant Gene Sets",
#   lwd = 1  # Line width for the circles
# )
# 
# # Plot the Venn diagram
# grid.draw(venn.plot)
