
#Setup an interactive session
#salloc --account=eon -t 0-10:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Updated conda env 12-2023
#module load miniconda3/23.1.0#conda activate archr2023_12

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
library(tibble)

#Additional setup
setwd("/Volumes/DataBox/MCS2023/Tx_Comp")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

projMCS7

getAvailableMatrices(projMCS7)

############################################

mat<-read.csv("./C1_byTx_zscores_2024-03-21.csv",header=TRUE,row.names=1,
              check.names=FALSE)

head(mat)

################

# Makes a heatmap
#pdf(file="C1_byTx_pheatmap_2024-03-25.pdf", width=15, height=30)
#pheatmap(mat,scale="row",
#         color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

###############



# Read the CSV file
mat <- read.csv("./C1_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(mat)

# Sort the data based on the mean of the rows
sorted_mat <- mat[order(rowMeans(mat), decreasing = TRUE), ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Original treatment identifiers for mapping
original_treatment_ids <- c("t1", "t2", "t3", "t4")  # Original identifiers

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(treatment_labels, levels = treatment_labels))
rownames(annotation_df) <- colnames(sorted_mat)

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("navy", "white", "red"))(50)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "lightgreen", "Ts" = "orange", "Ts+" = "pink"))

# Create a PDF for the heatmap plot
pdf(file = "C1_byTx_sorted_heatmap_2024-10-16.pdf", width = 7, height = 32)

# Create the pheatmap with updated treatment labels
pheatmap_result <- pheatmap(
  sorted_mat,
  scale = "row",  # Scale by row to normalize the data
  cluster_cols = FALSE,  # Disable column clustering
  show_colnames = TRUE,  # Show column names
  main = "Cluster 1 Heatmap of Gene Accessibility",
  color = color_palette,
  fontsize = 10, cellwidth = 10, cellheight = 10,
  annotation_col = annotation_df,  # Add updated treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE  # Disable displaying numbers in the heatmap
)

# Close the PDF device to finish plotting
dev.off()



###########################
###########################
###########################

## Ordered by treatment group


# Read the CSV file
mat <- read.csv("./C1_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(mat)

# Sort the data based on the t1 values in descending order
sorted_mat <- mat[order(mat$t1, decreasing = TRUE), ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(sorted_mat)

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("navy", "white", "red"))(50)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "lightgreen", "Ts" = "orange", "Ts+" = "pink"))

# Create a PDF for the heatmap plot
pdf(file = "C1_byTx_sorted_heatmap_by_t1_no_dendrogram_2024-10-16.pdf", width = 7, height = 32)

# Create the pheatmap with updated treatment labels without clustering
pheatmap_result <- pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Cluster 1 Heatmap of Gene Accessibility (Sorted by t1)",
  color = color_palette,
  fontsize = 10, cellwidth = 10, cellheight = 10,
  annotation_col = annotation_df,  # Add treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE  # Disable displaying numbers in the heatmap
)

# Close the PDF device to finish plotting
dev.off()



########################
########################
########################


## Sorting by 2 groups (incorporating the mean to determine significance)


# Read the CSV file
mat <- read.csv("./C2_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(mat)

# Calculate the mean of t1 and t2 for sorting
mean_values <- rowMeans(mat[, c("t1", "t2")], na.rm = TRUE)

# Sort the data based on the mean of t1 and t2 values in descending order
sorted_mat <- mat[order(mean_values, decreasing = TRUE), ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(sorted_mat)

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("navy", "white", "red"))(50)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "lightgreen", "Ts" = "orange", "Ts+" = "pink"))

# Create a PDF for the heatmap plot
pdf(file = "C2_byTx_sorted_heatmap_by_t1_t2_2024-10-16.pdf", width = 7, height = 32)

# Create the pheatmap with updated treatment labels without clustering
pheatmap_result <- pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Cluster 2 Heatmap of Gene Accessibility (Sorted by t1 and t2)",
  color = color_palette,
  fontsize = 10, cellwidth = 10, cellheight = 10,
  annotation_col = annotation_df,  # Add treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE  # Disable displaying numbers in the heatmap
)

# Close the PDF device to finish plotting
dev.off()




##################################
##################################
##################################
##################################


## Unable to run
## Do not have multiple accessibility values



## To highlight significant differences between groups, ex: t-tests: 


# Read the CSV file that contains raw accessibility values
# Assuming each row represents a gene and each column represents samples within treatment groups.
# You might need to reshape your data or combine data from multiple samples if they're in separate files.
# Example: Load raw data where each treatment group has multiple observations (this is hypothetical)
raw_data <- read.csv("./C1_raw_accessibility_data.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(raw_data)

# Initialize a vector for p-values
p_values <- numeric(nrow(raw_data))

# Perform t-tests for each gene between t1 and t2
for (i in 1:nrow(raw_data)) {
  # Assuming raw_data has separate columns for each sample in t1 and t2
  # Example column names: t1_sample1, t1_sample2, ..., t2_sample1, t2_sample2, ...
  
  # Extract samples for t1 and t2
  t1_samples <- raw_data[i, grep("^t1_", colnames(raw_data))]
  t2_samples <- raw_data[i, grep("^t2_", colnames(raw_data))]
  
  # Perform t-test on samples
  p_values[i] <- t.test(t1_samples, t2_samples)$p.value
}

# Adjust p-values for multiple testing using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create a data frame to store accessibility values and significance
results <- data.frame(
  Gene = rownames(raw_data),
  t1 = rowMeans(raw_data[, grep("^t1_", colnames(raw_data))]),  # Mean across samples in t1
  t2 = rowMeans(raw_data[, grep("^t2_", colnames(raw_data))]),  # Mean across samples in t2
  p_value = p_values,
  adjusted_p_value = adjusted_p_values
)

# Identify significant genes (e.g., adjusted p-value < 0.05)
results$sig <- results$adjusted_p_value < 0.05

# Sort the data based on the t1 values in descending order
sorted_results <- results[order(results$t1, decreasing = TRUE), ]

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("navy", "white", "red"))(50)

# Create a PDF for the heatmap plot
pdf(file = "C1_byTx_sorted_heatmap_with_significance_2024-10-16.pdf", width = 7, height = 32)

# Create the pheatmap, including the gene names
pheatmap_result <- pheatmap(
  raw_data[rownames(sorted_results), ],  # Use the original matrix, sorted by t1
  scale = "none",  # No scaling to maintain original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Cluster 1 Heatmap of Gene Accessibility (Sorted by t1)",
  color = color_palette,
  fontsize = 10, cellwidth = 10, cellheight = 10,
  display_numbers = sorted_results$sig,  # Display significance markers
  number_color = ifelse(sorted_results$sig, "red", "black")  # Change number color for significance
)

# Overlay significance markers on the heatmap
text(x = seq(1, ncol(raw_data)), 
     y = rep(nrow(sorted_results) + 1, ncol(raw_data)), 
     labels = ifelse(sorted_results$sig, "*", ""),  # Place asterisk for significant genes
     cex = 1.5, col = "red")

# Close the PDF device to finish plotting
dev.off()



##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################


## Remove boxes and sort based on t3 top 50%

setwd("/Volumes/DataBox/MCS2023/Tx_Comp")


# Read the CSV file
mat <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(mat)

# Sort the data based on the t3 values in descending order
top_50_percent_genes <- mat[order(mat$t3, decreasing = TRUE), ]

# Calculate the number of genes in the top 50%
num_top_genes <- floor(0.5 * nrow(top_50_percent_genes))

# Subset the matrix to include only the top 10% genes based on the 't3' column
top_genes_mat <- top_50_percent_genes[1:num_top_genes, ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(top_genes_mat)

# Create a continuous color gradient (smooth)
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)  # Increase number of colors for smooth gradient

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Create a PDF for the heatmap plot
pdf(file = "C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-05.pdf", width = 9, height = 32)

# Create the pheatmap with updated treatment labels without clustering
pheatmap_result <- pheatmap(
  top_genes_mat,  # Use only the top 50% genes
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Cluster 3 Heatmap of Gene Accessibility (Top 50% Genes for Ts)",
  color = color_palette,  # Use the continuous color gradient
  fontsize = 10, cellwidth = 30, cellheight = 10,
  annotation_col = annotation_df,  # Add treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  border_color = NA  # Remove the borders around the cells
)

# Close the PDF device to finish plotting
dev.off()


##################################
##################################
##################################
##################################

## Hierarchical clustering with dendrogram - 1 group


# Read the CSV file
mat <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of your data
str(mat)

# Sort the data based on the t3 values in descending order
top_50_percent_genes <- mat[order(mat$t3, decreasing = TRUE), ]

# Calculate the number of genes in the top 50%
num_top_genes <- floor(0.5 * nrow(top_50_percent_genes))

# Subset the matrix to include only the top 50% genes based on the 't3' column
top_genes_mat <- top_50_percent_genes[1:num_top_genes, ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(top_genes_mat)

# Create a continuous color gradient (smooth)
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)  # Increase number of colors for smooth gradient

# Create a color mapping for treatment groups
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Create a PDF for the heatmap plot
pdf(file = "C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-06_with_dendrograms.pdf", width = 9, height = 32)

# Create the pheatmap with clustering enabled
pheatmap_result <- pheatmap(
  top_genes_mat,  # Use only the top 50% genes
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = TRUE,  # Enable column clustering for hierarchical dendrogram
  cluster_rows = TRUE,  # Enable row clustering for hierarchical dendrogram
  show_colnames = TRUE,  # Show column names
  main = "Cluster 3 Heatmap of Gene Accessibility (Top 50% Genes for Ts)",
  color = color_palette,  # Use the continuous color gradient
  fontsize = 10, cellwidth = 30, cellheight = 10,
  annotation_col = annotation_df,  # Add treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  border_color = NA  # Remove the borders around the cells
)

# Close the PDF device to finish plotting
dev.off()




##################################
##################################
##################################
##################################

## Hierarchical clustering with NO dendrogram - 1 group, top 50%

# Read the CSV file
mat <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Sort the data based on the t3 values in descending order
top_50_percent_genes <- mat[order(mat$t3, decreasing = TRUE), ]

# Calculate the number of genes in the top 50%
num_top_genes <- floor(0.5 * nrow(top_50_percent_genes))

# Subset the matrix to include only the top 50% genes based on the 't3' column
top_genes_mat <- top_50_percent_genes[1:num_top_genes, ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(top_genes_mat)

# Create a continuous color gradient (smooth)
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)  # Smooth gradient

# Create a color mapping for treatment groups
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Perform hierarchical clustering without displaying the dendrogram
row_order <- hclust(dist(top_genes_mat))$order
col_order <- hclust(dist(t(top_genes_mat)))$order

# Reorder the matrix based on the clustering results
top_genes_mat <- top_genes_mat[row_order, col_order]

# Create a PDF for the heatmap plot
pdf(file = "C3_byTx_hierarchical_heatmap_1-Comp_Top50Percent_noDendrogramGenes_2024-11-10.pdf.pdf", width = 9, height = 32)

# Create the pheatmap without dendrograms and row names
pheatmap_result <- pheatmap(
  top_genes_mat,  # Use only the top 50% genes
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering to remove dendrogram display
  cluster_rows = FALSE,  # Disable row clustering to remove dendrogram display
  show_colnames = TRUE,  # Show column names
  show_rownames = FALSE,  # Hide gene names (row names)
  main = "Cluster 3 Heatmap of Gene Accessibility (Top 50% Genes for Ts)",
  color = color_palette,  # Use the continuous color gradient
  fontsize = 10, cellwidth = 30, cellheight = 10,
  annotation_col = annotation_df,  # Add treatment group labels
  annotation_colors = annotation_colors,  # Use colors for the annotations
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  border_color = NA  # Remove the borders around the cells
)

# Close the PDF device to finish plotting
dev.off()



##################################
##################################
##################################
##################################

## Hierarchical clustering with NO dendrogram - multiple groups, top 50% of variance

# Load required libraries
library(pheatmap)

# Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Select the top 50% most variable genes based on variance for each dataset
select_top_50_percent <- function(mat) {
  gene_variances <- apply(mat, 1, var)
  top_genes <- names(sort(gene_variances, decreasing = TRUE)[1:round(0.5 * length(gene_variances))])
  mat[top_genes, ]
}

top_genes1 <- select_top_50_percent(mat1)
top_genes2 <- select_top_50_percent(mat2)

# Apply hierarchical clustering on each dataset separately
row_order1 <- hclust(dist(top_genes1))$order
row_order2 <- hclust(dist(top_genes2))$order
col_order1 <- hclust(dist(t(top_genes1)))$order
col_order2 <- hclust(dist(t(top_genes2)))$order

# Reorder each matrix based on the hierarchical clustering
top_genes1 <- top_genes1[row_order1, col_order1]
top_genes2 <- top_genes2[row_order2, col_order2]

# Combine the matrices side by side with unique column names
combined_matrix <- cbind(top_genes1, top_genes2)
colnames(combined_matrix) <- c(paste0(colnames(top_genes1), "_C3"), paste0(colnames(top_genes2), "_C6"))

# Define treatment annotations for the columns
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Adjust color palette and breaks for better contrast and full range
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)
color_breaks <- seq(min(combined_matrix, na.rm = TRUE), max(combined_matrix, na.rm = TRUE), length.out = 101)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Plot the heatmap with adjusted settings
jpeg("hierarchical_heatmap_C3-C6_top50byVariance_2024-11-10.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  cluster_rows = FALSE,     # No row clustering to remove dendrogram display
  cluster_cols = FALSE,     # No column clustering to remove dendrogram display
  color = color_palette,    # Custom color palette
  breaks = color_breaks,    # Set color breaks to cover full range of data
  show_rownames = FALSE,    # Hide gene names
  show_colnames = FALSE,    # Hide column names for a cleaner plot
  main = "",                # Remove title
  legend = FALSE,           # Remove legend for a cleaner plot
  annotation_colors = annotation_colors,
  gaps_col = ncol(top_genes1)  # Add space between the two matrices
)
dev.off()



##################################
##################################
##################################
##################################

## Hierarchical clustering with NO dendrogram - multiple groups, top 50 of T3 genes


# Load required libraries
library(pheatmap)

# Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Step 1: Select top 50 genes based on their expression in T3 (from 't3' column)
select_top_50_t3 <- function(mat, treatment_col) {
  # Sort genes by their expression values in the 't3' column and select top 50
  top_genes <- mat[order(mat[, treatment_col], decreasing = TRUE), ][1:50, ]
  return(top_genes)
}

top_genes1_t3 <- select_top_50_t3(mat1, "t3")
top_genes2_t3 <- select_top_50_t3(mat2, "t3")

# Step 2: Apply hierarchical clustering on the selected top genes based on gene expression across all treatments
row_order1 <- hclust(dist(top_genes1_t3))$order  # Clustering rows based on their expression across all treatments
row_order2 <- hclust(dist(top_genes2_t3))$order  # Same for the second dataset

# Step 3: Reorder each matrix based on hierarchical clustering of genes (rows)
top_genes1_t3 <- top_genes1_t3[row_order1, ]
top_genes2_t3 <- top_genes2_t3[row_order2, ]

# Step 4: Combine the matrices side by side with unique column names
combined_matrix <- cbind(top_genes1_t3, top_genes2_t3)
colnames(combined_matrix) <- c(paste0(colnames(top_genes1_t3), "_C3"), paste0(colnames(top_genes2_t3), "_C6"))

# Step 5: Define treatment annotations for the columns
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Step 6: Adjust color palette and breaks for better contrast and full range
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)
color_breaks <- seq(min(combined_matrix, na.rm = TRUE), max(combined_matrix, na.rm = TRUE), length.out = 101)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Step 7: Plot the heatmap with adjusted settings
jpeg("hierarchical_heatmap_C3-C6_top50T3_2024-11-10.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  cluster_rows = FALSE,     # Perform row clustering to display dendrogram (based on gene expression)
  cluster_cols = FALSE,    # No column clustering (since we are keeping the column order fixed)
  color = color_palette,   # Custom color palette
  breaks = color_breaks,   # Set color breaks to cover full range of data
  show_rownames = FALSE,    # Show gene names (for the top 50 genes)
  show_colnames = FALSE,   # Hide column names for a cleaner plot
  main = "Top 50 Genes in T3 - Hierarchical Clustering",  # Title
  legend = TRUE,   # Display the legend
  annotation_colors = annotation_colors,
  gaps_col = ncol(top_genes1_t3)  # Add space between the two datasets
)
dev.off()



##################################
##################################
##################################
##################################

## Hierarchical clustering with NO dendrogram - multiple groups, all genes

## Due to differing numbers of genes in each comparison, 
## this code limits the genes displayed to include only the genes that are common in both datasets


# Step 1: Load required libraries
library(pheatmap)

# Step 2: Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Step 3: Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Step 4: Find common genes in both datasets
common_genes <- intersect(rownames(mat1), rownames(mat2))

# Step 5: Subset both matrices to include only the common genes
mat1 <- mat1[common_genes, ]
mat2 <- mat2[common_genes, ]

# Step 6: Apply hierarchical clustering on each subsetted dataset
row_order1 <- hclust(dist(mat1))$order
row_order2 <- hclust(dist(mat2))$order
col_order1 <- hclust(dist(t(mat1)))$order
col_order2 <- hclust(dist(t(mat2)))$order

# Reorder each matrix based on the hierarchical clustering
mat1 <- mat1[row_order1, col_order1]
mat2 <- mat2[row_order2, col_order2]

# Step 7: Combine the matrices side by side with unique column names
combined_matrix <- cbind(mat1, mat2)
colnames(combined_matrix) <- c(paste0(colnames(mat1), "_C3"), paste0(colnames(mat2), "_C6"))

# Step 8: Define treatment annotations for the columns
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Step 9: Plot the heatmap with hierarchical ordering for common genes, no dendrogram, and no gene names
jpeg("combined_heatmap_common_genes_no_dendrogram.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  cluster_rows = FALSE,  # No row clustering to remove dendrogram display
  cluster_cols = FALSE,  # No column clustering to remove dendrogram display
  show_rownames = FALSE,  # Hide gene names
  show_colnames = FALSE,  # Hide column names for a cleaner plot
  main = "",  # Remove title
  legend = FALSE,  # Remove legend for a cleaner plot
  gaps_col = ncol(mat1)  # Add space between the two matrices
)
dev.off()



##################################
##################################
##################################
##################################


## To display multiple heatmaps on the same page


# Load required libraries
library(pheatmap)
library(gridExtra)  # For arranging plots side by side

# Read the first CSV file
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Read the second CSV file
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Check the structure of both datasets
str(mat1)
str(mat2)

# Sort the data based on the t3 values in descending order for both datasets
top_50_percent_genes1 <- mat1[order(mat1$t3, decreasing = TRUE), ]
top_50_percent_genes2 <- mat2[order(mat2$t3, decreasing = TRUE), ]

# Calculate the number of genes in the top 50% for both datasets
num_top_genes1 <- floor(0.5 * nrow(top_50_percent_genes1))
num_top_genes2 <- floor(0.5 * nrow(top_50_percent_genes2))

# Subset the matrices to include only the top 50% genes
top_genes_mat1 <- top_50_percent_genes1[1:num_top_genes1, ]
top_genes_mat2 <- top_50_percent_genes2[1:num_top_genes2, ]

# Define updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels for both datasets
annotation_df1 <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df1) <- colnames(top_genes_mat1)

annotation_df2 <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df2) <- colnames(top_genes_mat2)

# Create a continuous color gradient (smooth)
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)  # Smooth gradient for heatmap colors

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Create the PDF device for the heatmaps
pdf(file = "Comparison_Heatmaps_2024-11-06.pdf", width = 18, height = 32)

# Create the first heatmap (for C3 dataset)
pheatmap_result1 <- pheatmap(
  top_genes_mat1,  # Use the top 50% genes from mat1
  scale = "none",  # No scaling to maintain original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering
  show_colnames = TRUE,  # Show column names
  main = "Cluster 3 Heatmap of Gene Accessibility (Top 50% Genes for Ts)",
  color = color_palette,  # Continuous color gradient
  fontsize = 10, cellwidth = 30, cellheight = 10,
  annotation_col = annotation_df1,  # Add treatment group labels for mat1
  annotation_colors = annotation_colors,  # Treatment colors
  display_numbers = FALSE,  # Disable displaying numbers
  border_color = NA  # Remove borders
)

# Create the second heatmap (for C6 dataset)
pheatmap_result2 <- pheatmap(
  top_genes_mat2,  # Use the top 50% genes from mat2
  scale = "none",  # No scaling to maintain original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering
  show_colnames = TRUE,  # Show column names
  main = "Cluster 6 Heatmap of Gene Accessibility (Top 50% Genes for Ts)",
  color = color_palette,  # Continuous color gradient
  fontsize = 10, cellwidth = 30, cellheight = 10,
  annotation_col = annotation_df2,  # Add treatment group labels for mat2
  annotation_colors = annotation_colors,  # Treatment colors
  display_numbers = FALSE,  # Disable displaying numbers
  border_color = NA  # Remove borders
)

# Arrange the two heatmaps side by side with no space between them
grid.arrange(
  pheatmap_result1$gtable, pheatmap_result2$gtable,
  ncol = 2,  # Number of columns
  widths = c(1, 1)  # Equal width for both heatmaps
)

# Close the PDF device to finish plotting
dev.off()
  
  
##################################
##################################
##################################
##################################


## To lay the heatmaps side by side with no spaces, gene names, titles, or legends
## Using the top 50 genes from each comparison

# Step 1: Load required libraries
library(pheatmap)

# Step 2: Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Step 3: Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Step 4: Take the top 50 genes based on 't3' column for both datasets
top_50_genes1 <- mat1[order(mat1$t3, decreasing = TRUE), ][1:50, ]
top_50_genes2 <- mat2[order(mat2$t3, decreasing = TRUE), ][1:50, ]

# Step 5: Create full matrices for both datasets with unique column names
full_mat1 <- matrix(NA, nrow = 50, ncol = ncol(mat1))
rownames(full_mat1) <- rownames(top_50_genes1)
colnames(full_mat1) <- paste0(colnames(mat1), "_C3")  # Add a suffix for uniqueness
full_mat1[rownames(top_50_genes1), ] <- as.matrix(top_50_genes1)

full_mat2 <- matrix(NA, nrow = 50, ncol = ncol(mat2))
rownames(full_mat2) <- rownames(top_50_genes2)
colnames(full_mat2) <- paste0(colnames(mat2), "_C6")  # Add a different suffix for uniqueness
full_mat2[rownames(top_50_genes2), ] <- as.matrix(top_50_genes2)

# Step 6: Combine the matrices side by side
combined_matrix <- cbind(full_mat1, full_mat2)

# Step 7: Convert combined_matrix to numeric, filling NA with 0 if necessary for visualization
combined_matrix[is.na(combined_matrix)] <- 0

# Step 8: Define desired column order based on actual column names in combined_matrix
desired_order <- c("t1_C3", "t2_C3", "t3_C3", "t4_C3", "t1_C6", "t2_C6", "t3_C6", "t4_C6")

# Reorder combined_matrix according to the desired order
combined_matrix <- combined_matrix[, desired_order]

# Step 9: Create treatment annotations for the reordered combined matrix
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Step 10: Plot the heatmap without clustering, gene names, and save as JPEG
jpeg("combined_heatmap_no_genes.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  main = "Combined Heatmap of C3 and C6",
  cluster_rows = FALSE,  # Disable row clustering to remove row dendrogram
  cluster_cols = FALSE,  # Disable column clustering to remove column dendrogram
  labels_row = FALSE     # Completely remove row labels (gene names)
)
dev.off()


##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################


## Hierarchical clustering with NO dendrogram - multiple groups, top 90 of T1 genes


# Load required libraries
library(pheatmap)

# Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Step 1: Select top 90 genes based on their expression in T1 (from 't1' column)
select_top_90_t1 <- function(mat, treatment_col) {
  # Sort genes by their expression values in the 't1' column and select top 200
  top_genes <- mat[order(mat[, treatment_col], decreasing = TRUE), ][1:90, ]
  return(top_genes)
}

top_genes1_t1 <- select_top_90_t1(mat1, "t1")
top_genes2_t1 <- select_top_90_t1(mat2, "t1")

# Step 2: Apply hierarchical clustering on the selected top genes based on gene expression across all treatments
row_order1 <- hclust(dist(top_genes1_t1))$order  # Clustering rows based on their expression across all treatments
row_order2 <- hclust(dist(top_genes2_t1))$order  # Same for the second dataset

# Step 3: Reorder each matrix based on hierarchical clustering of genes (rows)
top_genes1_t1 <- top_genes1_t1[row_order1, ]
top_genes2_t1 <- top_genes2_t1[row_order2, ]

# Step 4: Combine the matrices side by side with unique column names
combined_matrix <- cbind(top_genes1_t1, top_genes2_t1)
colnames(combined_matrix) <- c(paste0(colnames(top_genes1_t1), "_C3"), paste0(colnames(top_genes2_t1), "_C6"))

# Step 5: Define treatment annotations for the columns
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Step 6: Adjust color palette and breaks for better contrast and full range
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)
color_breaks <- seq(min(combined_matrix, na.rm = TRUE), max(combined_matrix, na.rm = TRUE), length.out = 101)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Step 7: Plot the heatmap with adjusted settings
jpeg("hierarchical_heatmap_C3-C6_top200-T1_2024-11-10.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  cluster_rows = FALSE,     # Perform row clustering to display dendrogram (based on gene expression)
  cluster_cols = FALSE,    # No column clustering (since we are keeping the column order fixed)
  color = color_palette,   # Custom color palette
  breaks = color_breaks,   # Set color breaks to cover full range of data
  show_rownames = FALSE,    # Show gene names (for the top 50 genes)
  show_colnames = FALSE,   # Hide column names for a cleaner plot
  main = "Top 200 of DA Genes in T1 - Hierarchical Clustering",  # Title
  legend = TRUE,   # Display the legend
  annotation_colors = annotation_colors,
  gaps_col = ncol(top_genes1_t1)  # Add space between the two datasets
)
dev.off()



##################################
##################################
##################################


## Hierarchical clustering with NO dendrogram - multiple groups, top 90% of T1 genes
## Cannot run based on differing number of rows between clusters

# Load required libraries
library(pheatmap)

# Define treatment labels
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")

# Read the datasets
mat1 <- read.csv("./C3_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)
mat2 <- read.csv("./C6_byTx_zscores_2024-03-21.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Step 1: Select top 90% genes based on their expression in T1 (from 't1' column)
select_top_90_percent_t1 <- function(mat, treatment_col) {
  # Calculate the 10th percentile of the 't1' column (exclude bottom 10%)
  percentile_cutoff <- quantile(mat[, treatment_col], 0.10, na.rm = TRUE)
  # Select genes with expression >= 10th percentile
  top_genes <- mat[mat[, treatment_col] >= percentile_cutoff, ]
  return(top_genes)
}

top_genes1_t1 <- select_top_90_percent_t1(mat1, "t1")
top_genes2_t1 <- select_top_90_percent_t1(mat2, "t1")

# Step 2: Apply hierarchical clustering on the selected top genes based on gene expression across all treatments
row_order1 <- hclust(dist(top_genes1_t1))$order  # Clustering rows based on their expression across all treatments
row_order2 <- hclust(dist(top_genes2_t1))$order  # Same for the second dataset

# Step 3: Reorder each matrix based on hierarchical clustering of genes (rows)
top_genes1_t1 <- top_genes1_t1[row_order1, ]
top_genes2_t1 <- top_genes2_t1[row_order2, ]

# Step 4: Combine the matrices side by side with unique column names
combined_matrix <- cbind(top_genes1_t1, top_genes2_t1)
colnames(combined_matrix) <- c(paste0(colnames(top_genes1_t1), "_C3"), paste0(colnames(top_genes2_t1), "_C6"))

# Step 5: Define treatment annotations for the columns
combined_treatments <- c(rep(c("2N", "2N+", "Ts", "Ts+"), 2))  # 2 for two datasets (C3 and C6)
combined_annotation <- data.frame(Treatment = factor(combined_treatments, levels = treatment_labels))
rownames(combined_annotation) <- colnames(combined_matrix)

# Step 6: Adjust color palette and breaks for better contrast and full range
color_palette <- colorRampPalette(c("darkgreen", "white", "purple"))(100)
color_breaks <- seq(min(combined_matrix, na.rm = TRUE), max(combined_matrix, na.rm = TRUE), length.out = 101)

# Create a color mapping for treatment groups (optional, adjust colors as needed)
annotation_colors <- list(Treatment = c("2N" = "lightblue", "2N+" = "yellow", "Ts" = "lightgreen", "Ts+" = "orange"))

# Step 7: Plot the heatmap with adjusted settings
jpeg("hierarchical_heatmap_C3-C6_top90percent-T1_2024-11-10.jpg", width = 800, height = 800, quality = 100)
pheatmap(
  combined_matrix,
  annotation_col = combined_annotation,
  cluster_rows = FALSE,     # Perform row clustering to display dendrogram (based on gene expression)
  cluster_cols = FALSE,    # No column clustering (since we are keeping the column order fixed)
  color = color_palette,   # Custom color palette
  breaks = color_breaks,   # Set color breaks to cover full range of data
  show_rownames = FALSE,    # Show gene names (for the top 50 genes)
  show_colnames = FALSE,   # Hide column names for a cleaner plot
  main = "Top 90% of DA Genes in T1 - Hierarchical Clustering",  # Title
  legend = TRUE,   # Display the legend
  annotation_colors = annotation_colors,
  gaps_col = ncol(top_genes1_t1)  # Add space between the two datasets
)
dev.off()




##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
################################## 
##################################
##################################
##################################

