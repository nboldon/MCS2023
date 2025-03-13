
## Prepare data for Lasso analysis



#Load libraries
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(Seurat)
library(Signac)
library(BiocManager)
library(BiocGenerics)
library(tidyr)
library(cowplot)
library(ggplot2)

#Additional setup
setwd("/Volumes/DataBox/MCS2023/Stats/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", 
                             force = FALSE, showLogo = FALSE)

projMCS6
getAvailableMatrices(projMCS6)


##################################################

#peakMat <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
#proj <- filterPeaks(proj, min.cells = 10)

##################################################
##################################################
##################################################


## Converts ArchR Project into gene matrix - DID NOT RUN



geneMat <- getMatrixFromProject(projMCS6, useMatrix = "GeneScoreMatrix")

# Convert geneMat to a data frame
geneMat_df <- as.data.frame(assay(geneMat))

# Save to a CSV file
write.csv(geneMat_df, "gene_score_matrix.csv", row.names = TRUE)

# Save other metadata if necessary (e.g., row or column data)
write.csv(as.data.frame(rowData(geneMat)), "gene_score_rowData.csv", row.names = TRUE)
write.csv(as.data.frame(colData(geneMat)), "gene_score_colData.csv", row.names = TRUE)



## Load the matrix back in later


# Load the matrix
geneMat_df <- read.csv("gene_score_matrix.csv", row.names = 1)

# Convert back to a SummarizedExperiment object if needed
library(SummarizedExperiment)

# Load metadata (if needed)
rowData_df <- read.csv("gene_score_rowData.csv", row.names = 1)
colData_df <- read.csv("gene_score_colData.csv", row.names = 1)

# Recreate SummarizedExperiment object
geneMat <- SummarizedExperiment(
  assays = list(counts = as.matrix(geneMat_df)),
  rowData = rowData_df,
  colData = colData_df
)


##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################



## To create scatterplots using gene scores by cluster



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)


# Load the gene score data
gene_scores <- read.csv("pivoted_glut_byTx_zscores_2025-01-26.csv", stringsAsFactors = FALSE)

# Ensure all treatment columns (except 'Gene') are numeric
# Replace non-numeric values with NA
gene_scores[, -1] <- lapply(gene_scores[, -1], function(x) as.numeric(as.character(x)))

# Load the behavioral data
behavior <- read.csv("TaskA-45_AllStats.csv", stringsAsFactors = FALSE)

# Reshape gene score data to long format
gene_scores_long <- gene_scores %>%
  pivot_longer(cols = -Gene, names_to = "Treatment", values_to = "Score")

# Prepare the list of behavioral tests to iterate over
behavioral_tests <- c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")

# Define the PDF file to save all plots
pdf("GeneScores_BehavioralScatterplots.pdf", width = 12, height = 8)

# Loop through each behavioral test
for (test in behavioral_tests) {
  
  # Merge gene scores with behavioral data
  merged_data <- merge(gene_scores_long, behavior, by = "Treatment")
  
  # Create all pairwise treatment comparisons
  treatment_pairs <- combn(unique(merged_data$Treatment), 2, simplify = FALSE)
  
  # Loop through each treatment pair
  for (pair in treatment_pairs) {
    x_group <- pair[1]
    y_group <- pair[2]
    
    # Filter data for the two treatment groups
    data_x <- merged_data %>% filter(Treatment == x_group)
    data_y <- merged_data %>% filter(Treatment == y_group)
    
    # Ensure only common genes are included
    common_genes <- intersect(data_x$Gene, data_y$Gene)
    data_x <- data_x %>% filter(Gene %in% common_genes)
    data_y <- data_y %>% filter(Gene %in% common_genes)
    
    # Merge by Gene to prepare for scatterplot
    scatter_data <- merge(data_x, data_y, by = "Gene", suffixes = c("_x", "_y"))
    
    # Create scatterplot
    scatter_plot <- ggplot(scatter_data, aes(x = Score_x, y = Score_y)) +
      geom_point(aes(color = scatter_data[[test]]), alpha = 0.7, size = 3) +
      labs(
        title = paste("Scatterplot of", y_group, "vs", x_group, "for", test),
        x = paste("Gene Score for", x_group),
        y = paste("Gene Score for", y_group),
        color = test
      ) +
      scale_color_viridis(option = "C") +
      theme_minimal() +
      theme(legend.position = "right")
    
    # Print the scatterplot to the PDF
    print(scatter_plot)
  }
}

# Close the PDF device
dev.off()


##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################







# Load data
mat <- read.csv("normalized_matrix.csv", row.names = 1)
response <- read.csv("phenotype_data.csv")$phenotype

# Fit Lasso model
lasso_model <- cv.glmnet(as.matrix(mat), response, alpha = 1)

# Check selected features
selected_features <- coef(lasso_model, s = "lambda.min")
print(selected_features)
