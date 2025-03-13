

## GO Genes of Interest


#################################################################


## Barplots for genes of interest count by cluster


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Define the path to your single comparison CSV file
file_path <- "/Volumes/DataBox/GO_Analysis/GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv"

# Get the current date for the output filename
date_today <- Sys.Date()

# Load the CSV file and replace empty strings and "NULL" values with NA
data <- read.csv(file_path, check.names = FALSE, na.strings = c("", "NULL"))

# Clean column names (only for clusters, not the gene names column)
colnames(data)[-1] <- trimws(colnames(data)[-1])  # Remove leading/trailing whitespace from cluster columns
colnames(data)[-1] <- gsub("^$", "Empty_Column", colnames(data)[-1])  # Replace empty cluster columns with a placeholder

# Create a long-format dataframe to better track gene presence across clusters
data_long <- data %>%
  pivot_longer(cols = -1, names_to = "Cluster", values_to = "Gene") %>%
  filter(!is.na(Gene))  # Remove NA values for clean comparison

# Remove rows where Cluster is NA (important for proper clustering)
data_long <- data_long %>% filter(!is.na(Cluster))

# Ensure correct sorting by extracting the numeric part of the cluster labels, and sort numerically
data_long$ClusterNum <- as.numeric(gsub("C", "", data_long$Cluster))  # Extract numeric values after "C"
valid_clusters <- data_long %>% filter(!is.na(ClusterNum))  # Keep only valid clusters with numbers

# Re-sort the clusters numerically based on the numeric part
sorted_clusters <- sort(unique(valid_clusters$ClusterNum))

# Reassign the correct order of the Cluster levels based on the sorted numeric order
data_long$Cluster <- factor(data_long$Cluster, levels = paste0("C", sorted_clusters))

# Identify duplicate genes
duplicate_genes <- data_long %>%
  group_by(Gene) %>%
  filter(n() > 1) %>%
  ungroup()

# Add a new column to mark duplicates
data_long <- data_long %>%
  mutate(DuplicateFlag = ifelse(Gene %in% duplicate_genes$Gene, "Duplicate", "Unique"))

# Count unique genes for each cluster (excluding duplicates)
gene_counts <- data_long %>%
  filter(DuplicateFlag == "Unique") %>%
  group_by(Cluster) %>%
  summarise(Gene_Count = n_distinct(Gene)) %>%
  ungroup()

# Display the final data structure before plotting
print("Data for plotting:")
print(gene_counts)

# Save the list of duplicate genes for later inspection
write.csv(duplicate_genes, file = paste0("GO-BP-GenesOfInterest_Duplicate_Genes_t1-t2_t3-t4_", date_today, ".csv"), row.names = FALSE)

# Plot the data and save as a JPG
output_file <- paste0("Cluster_Gene_Counts_t1-t2_t3-t4_", date_today, ".jpg")
plot <- ggplot(gene_counts, aes(x = Cluster, y = Gene_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "GO BP Genes of Interest Counts by Cluster - 2N/2N+ vs. Ts/Ts+",
       x = "Cluster",
       y = "Number of Unique Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 10, 10, 10))

# Save the plot as a JPG file
ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)


#####################################################
#####################################################
#####################################################


## GO Barplots for Genes of Interest


# Load necessary libraries
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  # Mouse genome annotations
library(gridExtra)  # For combining multiple plots
library(ggplot2)  # For theme function

# Load the CSV file containing genes of interest
file_path <- "/Volumes/DataBox/GO_Analysis/GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv"
gene_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Prepare a list to store the plots
all_go_plots <- list()

# Loop over each column in the gene_data (each cluster)
for (col_name in colnames(gene_data)) {
  # Extract the gene symbols for this cluster (remove any NA values)
  cluster_genes <- na.omit(gene_data[[col_name]])
  
  # Map gene symbols to Entrez gene IDs using org.Mm.eg.db
  entrez_genes <- tryCatch({
    bitr(cluster_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) { NULL })
  
  # Check if we successfully mapped any genes
  if (!is.null(entrez_genes) && nrow(entrez_genes) > 0) {
    
    # Remove any rows with NA in the Entrez ID column (if present)
    entrez_genes <- entrez_genes[!is.na(entrez_genes$ENTREZID), ]
    
    # If there are still valid genes, run GO enrichment
    if (nrow(entrez_genes) > 0) {
      go_enrichment <- enrichGO(gene = entrez_genes$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP",  # Biological Process
                                pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                                qvalueCutoff = 0.05)
      
      # Check if GO enrichment results are available
      if (!is.null(go_enrichment) && nrow(go_enrichment) > 0) {
        # Create a barplot for the GO enrichment result
        go_plot <- barplot(go_enrichment, showCategory = 10) + 
          ggtitle(paste("2N/2N+ vs Ts/Ts+ GO Genes of Interest - Cluster", col_name)) + 
          theme(
            plot.margin = margin(10, 40, 10, 10)  # Increase right margin for space
          )
        
        # Save the plot as a JPG file with adjusted size (widening the plot)
        plot_file <- paste0(col_name,"_GO_GenesOfInterest_byCluster_Barplots_t1-t2_t3-t4_2024-11-08.jpg")
        ggsave(plot_file, plot = go_plot, width = 12, height = 8, dpi = 300, device = "jpeg")  # Save as JPG
      }
    } else {
      message("No valid Entrez IDs found for cluster: ", col_name)
    }
  } else {
    message("Gene mapping failed for cluster: ", col_name)
  }
}

# t1-t2_t3-t4 gene mapping failed for C4, C6, C10, C12, C15, C16, C18, C19, C20, C21, C22, C25
# t2-t4_t1-t3 gene mapping failed for C5, C6, C10, C12, C13, C14, C16, C17, C18, C19, C20, C21, C22, C24, C25
# t3_t4 gene mapping failed for C4, C5, C6, C12, C14, C16, C17, C18, C19, C20, C21, C22, C24, C25


#####################################################
#####################################################
#####################################################


## GO Bubble Plots for Genes of Interest


# Load necessary libraries
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  # Mouse genome annotations
library(gridExtra)  # For combining multiple plots
library(ggplot2)  # For theme function

# Load the CSV file containing genes of interest
file_path <- "/Volumes/DataBox/GO_Analysis/GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv"
gene_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Prepare a list to store the plots
all_go_plots <- list()

# Loop over each column in the gene_data (each cluster)
for (col_name in colnames(gene_data)) {
  # Extract the gene symbols for this cluster (remove any NA values)
  cluster_genes <- na.omit(gene_data[[col_name]])
  
  # Map gene symbols to Entrez gene IDs using org.Mm.eg.db
  entrez_genes <- tryCatch({
    bitr(cluster_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) { NULL })
  
  # Check if we successfully mapped any genes
  if (!is.null(entrez_genes) && nrow(entrez_genes) > 0) {
    
    # Remove any rows with NA in the Entrez ID column (if present)
    entrez_genes <- entrez_genes[!is.na(entrez_genes$ENTREZID), ]
    
    # If there are still valid genes, run GO enrichment
    if (nrow(entrez_genes) > 0) {
      go_enrichment <- enrichGO(gene = entrez_genes$ENTREZID, 
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP",  # Biological Process
                                pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg
                                qvalueCutoff = 0.05)
      
      # Check if GO enrichment results are available
      if (!is.null(go_enrichment) && nrow(go_enrichment) > 0) {
        # Create a plot for the GO enrichment result
        go_plot <- dotplot(go_enrichment, showCategory = 10) + 
          ggtitle(paste("2N/2N+ vs. Ts/Ts+ GO Genes Of Interest - Cluster", col_name)) + 
          theme(
            plot.margin = margin(10, 40, 10, 10)  # Increase right margin for space
          )
        
        # Save the plot with adjusted size (widening the plot)
        plot_file <- paste0(col_name, "_GO_GenesOfInterest_Dotplot_t1-t2_t3-t4_2024-11-08.png")
        ggsave(plot_file, plot = go_plot, width = 12, height = 8, dpi = 300)  # Increase width to 12
      }
    } else {
      message("No valid Entrez IDs found for cluster: ", col_name)
    }
  } else {
    message("Gene mapping failed for cluster: ", col_name)
  }
}

# t1-t2_t3-t4 gene mapping failed for C4, C6, C10, C12, C15, C16, C18, C19, C20, C21, C22, C25
# t2-t4_t1-t3 gene mapping failed for C5, C6, C10, C12, C13, C14, C16, C17, C18, C19, C20, C21, C22, C24, C25
# t3_t4 gene mapping failed for C4, C5, C6, C12, C14, C16, C17, C18, C19, C20, C21, C22, C24, C25
