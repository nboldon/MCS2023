library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(viridis)
library(enrichplot)  

setwd("/Volumes/DataBox/GO_Analysis")


### Code adapted from 11_GO_Barplot_2.R
## Continues from code: /Volumes/DataBox/ProjMCS6/CellTypes_TvsTComp_2025-01-18.R
## also uses results from: /Volumes/DataBox/MCS2023/Tx_Comp/projMCS7_GeneMarkers_byCellType-TxGroup_2025-01-16.R

# Also see cellType level GO Analysis (instead of pairwise comparison as analyzed below):
## GO Genes of Interest generated from: /Volumes/DataBox/GO_Analysis/GO_Analysis_2025-01-18.R
##    saved as: /Volumes/DataBox/GO_Analysis/GO_Enrichment_Results_2025-01-18.csv


# GO terms
#BP = Biological Process
#MP = Molecular function
#CC = Cellular component


####################################
####################################
####################################


## Complete for each cell type and treatment comparison:


library(clusterProfiler)
library(org.Mm.eg.db)  # Required for OrgDb
library(viridis)       # For viridis color palette
library(enrichplot)    # For enhanced visualization
library(ggplot2)       # For custom plotting

# Load the data
glut_T1vsT4 <- read.csv("/Volumes/DataBox/MCS2023/Tx_Comp/glut_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv")

# Filter genes based on Log2FC
genes_to_test <- glut_T1vsT4[abs(glut_T1vsT4$Log2FC) >= 0.5, "name"]

# Perform GO enrichment analysis
GO_results <- enrichGO(
  gene = genes_to_test,
  OrgDb = "org.Mm.eg.db",
  keyType = "SYMBOL",
  ont = "BP"
)

# Save the GO results as a DataFrame
GO_results_df <- as.data.frame(GO_results)
write.csv(GO_results_df, "glut_T1vsT4_GO_results_2025-01-19.csv", row.names = FALSE)

# Extract top 20 categories for plotting
top_GO <- GO_results_df[1:20, ]  # Assuming GO_results_df is ordered by significance

# Define variables for dynamic title
cell_type <- "Glut"                # Cell type
treatment_comparison <- "T1 vs T4" # Treatment group comparison

# Create a barplot with ggplot2 using the viridis palette
ggplot(top_GO, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "D", name = "Adjusted p-value") +
  coord_flip() +
  labs(
    title = paste("Top 20 GO Terms (Biological Process) for", cell_type, "Cells in", treatment_comparison),
    x = "GO Term",
    y = "Gene Count"
  ) +
  theme_minimal()

# Save the plot as a PNG file
ggsave(
  paste0(cell_type, "_", treatment_comparison, "_top20_2025-01-19.png"),
  width = 10, height = 8, dpi = 300
)





####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################


## Loop to run through all cell types and treatment group comparisons



# Define cell types and treatment comparisons
cell_types <- c("glut", "gaba", "astro", "microglia", "oligo", "endo")  # Example cell types
treatment_comparisons <- c("T1vsT2", "T1vsT3", "T1vsT4", "T2vsT1", "T2vsT3", "T2vsT4",
                           "T3vsT1", "T3vsT2", "T3vsT4", "T4vsT1", "T4vsT2", "T4vsT3")  # Example treatment comparisons

# Loop through each combination of cell type and treatment comparison
for (cell_type in cell_types) {
  for (treatment_comparison in treatment_comparisons) {
    # Construct file path dynamically based on cell type and treatment comparison
    file_name <- paste0(
      "/Volumes/DataBox/MCS2023/Tx_Comp/",
      tolower(cell_type), "_",
      gsub(" ", "", treatment_comparison), "Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv"
    )
    
    # Check if the file exists
    if (!file.exists(file_name)) {
      message("File not found: ", file_name)
      next
    }
    
    # Load the data
    data <- read.csv(file_name)
    
    # Filter genes based on Log2FC
    genes_to_test <- data[abs(data$Log2FC) >= 0.5, "name"]
    
    # Skip if no genes meet the threshold
    if (length(genes_to_test) == 0) {
      message("No genes meet the Log2FC threshold for ", cell_type, " in ", treatment_comparison)
      next
    }
    
    # Perform GO enrichment analysis
    GO_results <- enrichGO(
      gene = genes_to_test,
      OrgDb = "org.Mm.eg.db",
      keyType = "SYMBOL",
      ont = "BP"
    )
    
    # Check if GO enrichment returned any results
    if (is.null(GO_results) || nrow(as.data.frame(GO_results)) == 0) {
      message("No GO terms found for ", cell_type, " in ", treatment_comparison)
      next
    }
    
    # Save the GO results as a DataFrame
    GO_results_df <- as.data.frame(GO_results)
    output_csv <- paste0(
      tolower(cell_type), "_", 
      gsub(" ", "", treatment_comparison), "_GO_results_2025-01-19.csv"
    )
    write.csv(GO_results_df, output_csv, row.names = FALSE)
    message("Saved GO results for ", cell_type, " in ", treatment_comparison, " to ", output_csv)
    
    # Extract top 20 categories for plotting
    top_GO <- GO_results_df[1:min(20, nrow(GO_results_df)), ]
    
    # Create a barplot with ggplot2 using the viridis palette
    plot <- ggplot(top_GO, aes(x = reorder(Description, -Count), y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis(option = "D", name = "Adjusted p-value") +
      coord_flip() +
      labs(
        title = paste("Top 20 GO Terms (Biological Process) for", cell_type, "Cells in", treatment_comparison),
        x = "GO Term",
        y = "Gene Count"
      ) +
      theme_minimal()
    
    # Save the plot as a PNG file
    output_png <- paste0(
      tolower(cell_type), "_", 
      gsub(" ", "", treatment_comparison), "_top20_viridis_2025-01-19.png"
    )
    ggsave(output_png, plot, width = 10, height = 8, dpi = 300)
    message("Saved plot for ", cell_type, " in ", treatment_comparison, " to ", output_png)
  }
}


######################################################
######################################################
######################################################


## Returned results

No genes meet the Log2FC threshold for glut in T1vsT2
No GO terms found for glut in T1vsT3
Saved GO results for glut in T1vsT4 to glut_T1vsT4_GO_results_2025-01-19.csv
Saved plot for glut in T1vsT4 to glut_T1vsT4_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for glut in T2vsT1
No GO terms found for glut in T2vsT3
No GO terms found for glut in T2vsT4
Saved GO results for glut in T3vsT1 to glut_T3vsT1_GO_results_2025-01-19.csv
Saved plot for glut in T3vsT1 to glut_T3vsT1_top20_viridis_2025-01-19.png
Saved GO results for glut in T3vsT2 to glut_T3vsT2_GO_results_2025-01-19.csv
Saved plot for glut in T3vsT2 to glut_T3vsT2_top20_viridis_2025-01-19.png
Saved GO results for glut in T3vsT4 to glut_T3vsT4_GO_results_2025-01-19.csv
Saved plot for glut in T3vsT4 to glut_T3vsT4_top20_viridis_2025-01-19.png
Saved GO results for glut in T4vsT1 to glut_T4vsT1_GO_results_2025-01-19.csv
Saved plot for glut in T4vsT1 to glut_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for glut in T4vsT2 to glut_T4vsT2_GO_results_2025-01-19.csv
Saved plot for glut in T4vsT2 to glut_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for glut in T4vsT3
No genes meet the Log2FC threshold for gaba in T1vsT2
No GO terms found for gaba in T1vsT3
Saved GO results for gaba in T1vsT4 to gaba_T1vsT4_GO_results_2025-01-19.csv
Saved plot for gaba in T1vsT4 to gaba_T1vsT4_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for gaba in T2vsT1
No GO terms found for gaba in T2vsT3
Saved GO results for gaba in T2vsT4 to gaba_T2vsT4_GO_results_2025-01-19.csv
Saved plot for gaba in T2vsT4 to gaba_T2vsT4_top20_viridis_2025-01-19.png
No GO terms found for gaba in T3vsT1
No GO terms found for gaba in T3vsT2
No genes meet the Log2FC threshold for gaba in T3vsT4
Saved GO results for gaba in T4vsT1 to gaba_T4vsT1_GO_results_2025-01-19.csv
Saved plot for gaba in T4vsT1 to gaba_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for gaba in T4vsT2 to gaba_T4vsT2_GO_results_2025-01-19.csv
Saved plot for gaba in T4vsT2 to gaba_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for gaba in T4vsT3
No genes meet the Log2FC threshold for astro in T1vsT2
Saved GO results for astro in T1vsT3 to astro_T1vsT3_GO_results_2025-01-19.csv
Saved plot for astro in T1vsT3 to astro_T1vsT3_top20_viridis_2025-01-19.png
Saved GO results for astro in T1vsT4 to astro_T1vsT4_GO_results_2025-01-19.csv
Saved plot for astro in T1vsT4 to astro_T1vsT4_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for astro in T2vsT1
Saved GO results for astro in T2vsT3 to astro_T2vsT3_GO_results_2025-01-19.csv
Saved plot for astro in T2vsT3 to astro_T2vsT3_top20_viridis_2025-01-19.png
Saved GO results for astro in T2vsT4 to astro_T2vsT4_GO_results_2025-01-19.csv
Saved plot for astro in T2vsT4 to astro_T2vsT4_top20_viridis_2025-01-19.png
Saved GO results for astro in T3vsT1 to astro_T3vsT1_GO_results_2025-01-19.csv
Saved plot for astro in T3vsT1 to astro_T3vsT1_top20_viridis_2025-01-19.png
No GO terms found for astro in T3vsT2
No genes meet the Log2FC threshold for astro in T3vsT4
Saved GO results for astro in T4vsT1 to astro_T4vsT1_GO_results_2025-01-19.csv
Saved plot for astro in T4vsT1 to astro_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for astro in T4vsT2 to astro_T4vsT2_GO_results_2025-01-19.csv
Saved plot for astro in T4vsT2 to astro_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for astro in T4vsT3
No genes meet the Log2FC threshold for microglia in T1vsT2
Saved GO results for microglia in T1vsT3 to microglia_T1vsT3_GO_results_2025-01-19.csv
Saved plot for microglia in T1vsT3 to microglia_T1vsT3_top20_viridis_2025-01-19.png
No GO terms found for microglia in T1vsT4
No genes meet the Log2FC threshold for microglia in T2vsT1
Saved GO results for microglia in T2vsT3 to microglia_T2vsT3_GO_results_2025-01-19.csv
Saved plot for microglia in T2vsT3 to microglia_T2vsT3_top20_viridis_2025-01-19.png
No GO terms found for microglia in T2vsT4
Saved GO results for microglia in T3vsT1 to microglia_T3vsT1_GO_results_2025-01-19.csv
Saved plot for microglia in T3vsT1 to microglia_T3vsT1_top20_viridis_2025-01-19.png
Saved GO results for microglia in T3vsT2 to microglia_T3vsT2_GO_results_2025-01-19.csv
Saved plot for microglia in T3vsT2 to microglia_T3vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for microglia in T3vsT4
Saved GO results for microglia in T4vsT1 to microglia_T4vsT1_GO_results_2025-01-19.csv
Saved plot for microglia in T4vsT1 to microglia_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for microglia in T4vsT2 to microglia_T4vsT2_GO_results_2025-01-19.csv
Saved plot for microglia in T4vsT2 to microglia_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for microglia in T4vsT3
Saved GO results for oligo in T1vsT2 to oligo_T1vsT2_GO_results_2025-01-19.csv
Saved plot for oligo in T1vsT2 to oligo_T1vsT2_top20_viridis_2025-01-19.png
Saved GO results for oligo in T1vsT3 to oligo_T1vsT3_GO_results_2025-01-19.csv
Saved plot for oligo in T1vsT3 to oligo_T1vsT3_top20_viridis_2025-01-19.png
No GO terms found for oligo in T1vsT4
Saved GO results for oligo in T2vsT1 to oligo_T2vsT1_GO_results_2025-01-19.csv
Saved plot for oligo in T2vsT1 to oligo_T2vsT1_top20_viridis_2025-01-19.png
No GO terms found for oligo in T2vsT3
Saved GO results for oligo in T2vsT4 to oligo_T2vsT4_GO_results_2025-01-19.csv
Saved plot for oligo in T2vsT4 to oligo_T2vsT4_top20_viridis_2025-01-19.png
Saved GO results for oligo in T3vsT1 to oligo_T3vsT1_GO_results_2025-01-19.csv
Saved plot for oligo in T3vsT1 to oligo_T3vsT1_top20_viridis_2025-01-19.png
Saved GO results for oligo in T3vsT2 to oligo_T3vsT2_GO_results_2025-01-19.csv
Saved plot for oligo in T3vsT2 to oligo_T3vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for oligo in T3vsT4
Saved GO results for oligo in T4vsT1 to oligo_T4vsT1_GO_results_2025-01-19.csv
Saved plot for oligo in T4vsT1 to oligo_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for oligo in T4vsT2 to oligo_T4vsT2_GO_results_2025-01-19.csv
Saved plot for oligo in T4vsT2 to oligo_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for oligo in T4vsT3
No genes meet the Log2FC threshold for endo in T1vsT2
Saved GO results for endo in T1vsT3 to endo_T1vsT3_GO_results_2025-01-19.csv
Saved plot for endo in T1vsT3 to endo_T1vsT3_top20_viridis_2025-01-19.png
Saved GO results for endo in T1vsT4 to endo_T1vsT4_GO_results_2025-01-19.csv
Saved plot for endo in T1vsT4 to endo_T1vsT4_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for endo in T2vsT1
Saved GO results for endo in T2vsT3 to endo_T2vsT3_GO_results_2025-01-19.csv
Saved plot for endo in T2vsT3 to endo_T2vsT3_top20_viridis_2025-01-19.png
Saved GO results for endo in T2vsT4 to endo_T2vsT4_GO_results_2025-01-19.csv
Saved plot for endo in T2vsT4 to endo_T2vsT4_top20_viridis_2025-01-19.png
Saved GO results for endo in T3vsT1 to endo_T3vsT1_GO_results_2025-01-19.csv
Saved plot for endo in T3vsT1 to endo_T3vsT1_top20_viridis_2025-01-19.png
Saved GO results for endo in T3vsT2 to endo_T3vsT2_GO_results_2025-01-19.csv
Saved plot for endo in T3vsT2 to endo_T3vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for endo in T3vsT4
Saved GO results for endo in T4vsT1 to endo_T4vsT1_GO_results_2025-01-19.csv
Saved plot for endo in T4vsT1 to endo_T4vsT1_top20_viridis_2025-01-19.png
Saved GO results for endo in T4vsT2 to endo_T4vsT2_GO_results_2025-01-19.csv
Saved plot for endo in T4vsT2 to endo_T4vsT2_top20_viridis_2025-01-19.png
No genes meet the Log2FC threshold for endo in T4vsT3




####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################


## Generate GO terms and plots for significant DA genes in all cell types and grouped treatment comparisons (not pairwise comps)


library(clusterProfiler)
library(org.Mm.eg.db)
library(viridis)
library(enrichplot)
library(ggplot2)




# Load the data
file_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/Combined_CellType_significant_t1-t2_t3-t4_comparison_2025-01-18.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE)

# Loop through each cell type column
for (cell_type in colnames(data)) {
  # Sanitize the cell type name to replace '+' with 'Plus' and '/' or spaces with '_'
  sanitized_cell_type <- gsub("[/ ]", "_", cell_type)  # Replace '/' and spaces with '_'
  sanitized_cell_type <- gsub("\\+", "Plus", sanitized_cell_type)  # Replace '+' with 'Plus'
  
  # Extract genes and their Log2FC values for the current cell type
  genes_to_test <- rownames(data)[!is.na(data[[cell_type]])]
  log2fc_values <- data[[cell_type]][!is.na(data[[cell_type]])]
  gene_log2fc <- setNames(log2fc_values, genes_to_test)
  
  # Skip if there are no genes to analyze
  if (length(genes_to_test) == 0) {
    message("No significant genes for ", cell_type)
    next
  }
  
  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    gene = names(gene_log2fc),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  # Skip if no GO terms are found
  if (nrow(as.data.frame(GO_results)) == 0) {
    message("No GO terms found for ", cell_type)
    next
  }
  
  # Prepare the GO results for plotting
  GO_results_df <- as.data.frame(GO_results)
  GO_results_df$Log2FC <- sapply(GO_results_df$geneID, function(genes) {
    # Extract and average Log2FC values for genes in each GO term
    gene_list <- unlist(strsplit(genes, "/"))
    mean(gene_log2fc[gene_list], na.rm = TRUE)
  })
  
  # Save the GO results to a CSV file
  output_csv_file <- paste0(sanitized_cell_type, "_GO_results.csv")
  write.csv(GO_results_df, output_csv_file, row.names = FALSE)
  message("Saved GO terms for ", cell_type, " as ", output_csv_file)
  
  # Create the bar plot
  p <- ggplot(GO_results_df, aes(
    x = Count, 
    y = reorder(Description, Count), 
    fill = Log2FC
  )) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "D", name = "Log2FC") +
    labs(
      title = paste("GO Enrichment for", cell_type),
      x = "Number of Genes",
      y = "GO Term"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # Save the plot
  output_file <- paste0(sanitized_cell_type, "_GO_enrichment_plot.pdf")
  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Saved plot for ", cell_type, " as ", output_file)
}


# No GO terms found for glut_DA_2N_2N+_Ts_Ts+
# No GO terms found for gaba_DA_2N_2N+_Ts_Ts+


####################################
####################################
####################################


# Load the data
file_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/Combined_CellType_significant_t2-t4_t1-t3_comparison_2025-01-18.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE)

# Loop through each cell type column
for (cell_type in colnames(data)) {
  # Sanitize the cell type name to replace '+' with 'Plus' and '/' or spaces with '_'
  sanitized_cell_type <- gsub("[/ ]", "_", cell_type)  # Replace '/' and spaces with '_'
  sanitized_cell_type <- gsub("\\+", "Plus", sanitized_cell_type)  # Replace '+' with 'Plus'
  
  # Extract genes and their Log2FC values for the current cell type
  genes_to_test <- rownames(data)[!is.na(data[[cell_type]])]
  log2fc_values <- data[[cell_type]][!is.na(data[[cell_type]])]
  gene_log2fc <- setNames(log2fc_values, genes_to_test)
  
  # Skip if there are no genes to analyze
  if (length(genes_to_test) == 0) {
    message("No significant genes for ", cell_type)
    next
  }
  
  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    gene = names(gene_log2fc),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  # Skip if no GO terms are found
  if (nrow(as.data.frame(GO_results)) == 0) {
    message("No GO terms found for ", cell_type)
    next
  }
  
  # Prepare the GO results for plotting
  GO_results_df <- as.data.frame(GO_results)
  GO_results_df$Log2FC <- sapply(GO_results_df$geneID, function(genes) {
    # Extract and average Log2FC values for genes in each GO term
    gene_list <- unlist(strsplit(genes, "/"))
    mean(gene_log2fc[gene_list], na.rm = TRUE)
  })
  
  # Save the GO results to a CSV file
  output_csv_file <- paste0(sanitized_cell_type, "_GO_results.csv")
  write.csv(GO_results_df, output_csv_file, row.names = FALSE)
  message("Saved GO terms for ", cell_type, " as ", output_csv_file)
  
  # Create the bar plot
  p <- ggplot(GO_results_df, aes(
    x = Count, 
    y = reorder(Description, Count), 
    fill = Log2FC
  )) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "D", name = "Log2FC") +
    labs(
      title = paste("GO Enrichment for", cell_type),
      x = "Number of Genes",
      y = "GO Term"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # Save the plot
  output_file <- paste0(sanitized_cell_type, "_GO_enrichment_plot.pdf")
  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Saved plot for ", cell_type, " as ", output_file)
}


# No GO terms found for glut_DA_2N+/Ts+_2N/Ts
# No GO terms found for oligo_DA_2N+/Ts+_2N/Ts
# No GO terms found for endo_DA_2N+/Ts+_2N/Ts
# No GO terms found for astro_DA_2N+/Ts+_2N/Ts


####################################
####################################
####################################


# Load the data
file_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/Combined_CellType_significant_t3_t4_comparison_2025-01-18.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE)

# Loop through each cell type column
for (cell_type in colnames(data)) {
  # Sanitize the cell type name to replace '+' with 'Plus' and '/' or spaces with '_'
  sanitized_cell_type <- gsub("[/ ]", "_", cell_type)  # Replace '/' and spaces with '_'
  sanitized_cell_type <- gsub("\\+", "Plus", sanitized_cell_type)  # Replace '+' with 'Plus'
  
  # Extract genes and their Log2FC values for the current cell type
  genes_to_test <- rownames(data)[!is.na(data[[cell_type]])]
  log2fc_values <- data[[cell_type]][!is.na(data[[cell_type]])]
  gene_log2fc <- setNames(log2fc_values, genes_to_test)
  
  # Skip if there are no genes to analyze
  if (length(genes_to_test) == 0) {
    message("No significant genes for ", cell_type)
    next
  }
  
  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    gene = names(gene_log2fc),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  # Skip if no GO terms are found
  if (nrow(as.data.frame(GO_results)) == 0) {
    message("No GO terms found for ", cell_type)
    next
  }
  
  # Prepare the GO results for plotting
  GO_results_df <- as.data.frame(GO_results)
  GO_results_df$Log2FC <- sapply(GO_results_df$geneID, function(genes) {
    # Extract and average Log2FC values for genes in each GO term
    gene_list <- unlist(strsplit(genes, "/"))
    mean(gene_log2fc[gene_list], na.rm = TRUE)
  })
  
  # Save the GO results to a CSV file
  output_csv_file <- paste0(sanitized_cell_type, "_GO_results.csv")
  write.csv(GO_results_df, output_csv_file, row.names = FALSE)
  message("Saved GO terms for ", cell_type, " as ", output_csv_file)
  
  # Create the bar plot
  p <- ggplot(GO_results_df, aes(
    x = Count, 
    y = reorder(Description, Count), 
    fill = Log2FC
  )) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "D", name = "Log2FC") +
    labs(
      title = paste("GO Enrichment for", cell_type),
      x = "Number of Genes",
      y = "GO Term"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # Save the plot
  output_file <- paste0(sanitized_cell_type, "_GO_enrichment_plot.pdf")
  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Saved plot for ", cell_type, " as ", output_file)
}


# No GO terms found for glut_DA_Ts-Ts+
# No GO terms found for astro_DA_Ts-Ts+


####################################
####################################
####################################


# Load the data
file_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/Combined_CellType_significant_t3_vs_all_2025-01-19.csv"
data <- read.csv(file_path, row.names = 1, check.names = FALSE)

# Loop through each cell type column
for (cell_type in colnames(data)) {
  # Sanitize the cell type name to replace '+' with 'Plus' and '/' or spaces with '_'
  sanitized_cell_type <- gsub("[/ ]", "_", cell_type)  # Replace '/' and spaces with '_'
  sanitized_cell_type <- gsub("\\+", "Plus", sanitized_cell_type)  # Replace '+' with 'Plus'
  
  # Extract genes and their Log2FC values for the current cell type
  genes_to_test <- rownames(data)[!is.na(data[[cell_type]])]
  log2fc_values <- data[[cell_type]][!is.na(data[[cell_type]])]
  gene_log2fc <- setNames(log2fc_values, genes_to_test)
  
  # Skip if there are no genes to analyze
  if (length(genes_to_test) == 0) {
    message("No significant genes for ", cell_type)
    next
  }
  
  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    gene = names(gene_log2fc),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  # Skip if no GO terms are found
  if (nrow(as.data.frame(GO_results)) == 0) {
    message("No GO terms found for ", cell_type)
    next
  }
  
  # Prepare the GO results for plotting
  GO_results_df <- as.data.frame(GO_results)
  GO_results_df$Log2FC <- sapply(GO_results_df$geneID, function(genes) {
    # Extract and average Log2FC values for genes in each GO term
    gene_list <- unlist(strsplit(genes, "/"))
    mean(gene_log2fc[gene_list], na.rm = TRUE)
  })
  
  # Save the GO results to a CSV file
  output_csv_file <- paste0(sanitized_cell_type, "_GO_results.csv")
  write.csv(GO_results_df, output_csv_file, row.names = FALSE)
  message("Saved GO terms for ", cell_type, " as ", output_csv_file)
  
  # Create the bar plot
  p <- ggplot(GO_results_df, aes(
    x = Count, 
    y = reorder(Description, Count), 
    fill = Log2FC
  )) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "D", name = "Log2FC") +
    labs(
      title = paste("GO Enrichment for", cell_type),
      x = "Number of Genes",
      y = "GO Term"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
  
  # Save the plot
  output_file <- paste0(sanitized_cell_type, "_GO_enrichment_plot.pdf")
  ggsave(output_file, plot = p, width = 10, height = 6)
  message("Saved plot for ", cell_type, " as ", output_file)
}



# No GO terms found for glut_DA_T3_vs_T1
# No GO terms found for endo_DA_T3_vs_T1
# No GO terms found for glut_DA_T3_vs_T2
# No GO terms found for gaba_DA_T3_vs_T2
# No GO terms found for glut_DA_T3_vs_T4
# No GO terms found for astro_DA_T3_vs_T4




##################################################################
##################################################################
##################################################################


## To combine spreadsheets from loop:


# Load necessary library
library(dplyr)

# Set the directory where the result files are stored
output_dir <- "/Volumes/DataBox/GO_Analysis/CellTypes"

# List all the .csv files in the output directory
csv_files <- list.files(output_dir, pattern = "\\_GO_results.csv$", full.names = TRUE)

# Create an empty list to store the data frames
all_data <- list()

# Loop through each .csv file and read it into a data frame
for (file in csv_files) {
  # Read the data from the CSV file
  data <- read.csv(file)
  
  # Ensure 'geneID' is always treated as a character (to handle data type inconsistencies)
  data$geneID <- as.character(data$geneID)
  
  # Assign the full file name to the new file_name column
  data$file_name <- basename(file)
  
  # Append the data frame to the list
  all_data <- append(all_data, list(data))
  
  # Debug: Check the number of rows being added
  message(paste("Rows in file", basename(file), ":", nrow(data)))
}

# Combine all the data frames into one
combined_data <- bind_rows(all_data)

# Debug: Check the number of rows after combining
message(paste("Total rows after combining:", nrow(combined_data)))

# Optionally, save the combined data to a new CSV file
output_file <- "/Volumes/DataBox/GO_Combined_Results.csv"
write.csv(combined_data, output_file, row.names = FALSE)

# Check if the file is saved successfully
if(file.exists(output_file)) {
  message("The combined file has been saved to: ", output_file)
} else {
  message("Failed to save the combined file.")
}

# Check the first few rows of the combined data
head(combined_data)







##################################################################
##################################################################
##################################################################
