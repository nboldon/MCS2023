

## clusterProfiler GO with Simplify function applied


# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(readr)         # For reading CSV files
library(dplyr)         # For data manipulation

# Step 1: Load data
file_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/Combined_Cluster_significant_t2-t4_t1-t3_comparison_2024-10-18.csv"
data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) # Preserve original column names

# Step 2: Filter genes for each column
gene_lists <- list()
for (col in colnames(data)[-1]) {  # Skip the first column (Gene names)
  filtered_genes <- data %>%
    filter(!is.na(.data[[col]])) %>%  # Select rows where the column is not NA
    pull(Gene)  # Extract the gene names
  
  cat("Column:", col, "- Number of genes:", length(filtered_genes), "\n")  # Debugging output
  gene_lists[[col]] <- filtered_genes  # Store in a list with column name as key
}

# Step 3: Run GO analysis and simplify results for each comparison
go_results <- list()
simplified_go_results <- list()  # To store simplified results
for (comparison in names(gene_lists)) {
  genes <- gene_lists[[comparison]]
  
  if (length(genes) > 0) {  # Ensure there are genes to analyze
    go <- tryCatch({
      enrichGO(
        gene = genes,
        OrgDb = org.Mm.eg.db,  # Use the correct organism database
        keyType = "SYMBOL",    # Replace with the appropriate key type if needed
        ont = "ALL",           # Can be "BP", "CC", or "MF"
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
    }, error = function(e) {
      cat("Error in enrichGO for", comparison, ":", conditionMessage(e), "\n")
      NULL
    })
    
    go_results[[comparison]] <- go
    
    # Simplify the GO results
    if (!is.null(go)) {
      simplified_go <- tryCatch({
        simplify(
          go,
          cutoff = 0.7,  # Adjust similarity cutoff as needed (default: 0.7)
          by = "p.adjust",  # Simplify based on adjusted p-value
          select_fun = min  # Select the GO term with the minimum p-value
        )
      }, error = function(e) {
        cat("Error in simplify for", comparison, ":", conditionMessage(e), "\n")
        NULL
      })
      simplified_go_results[[comparison]] <- simplified_go
    } else {
      simplified_go_results[[comparison]] <- NULL
    }
  } else {
    cat("No genes to analyze for", comparison, "\n")
    go_results[[comparison]] <- NULL  # Handle cases with no genes
    simplified_go_results[[comparison]] <- NULL
  }
}

# Step 4: Save or visualize the results
output_folder <- "GO_Dotplots"  # Define an output folder for dotplots
dir.create(output_folder, showWarnings = FALSE)  # Create the folder if it doesn't exist

# Helper function to sanitize column names for file names
sanitize_name <- function(name) {
  name <- gsub("\\/", "-", name)  # Replace forward slashes with dashes
  name <- gsub("\\+", "plus", name)  # Replace plus signs with "plus"
  name <- gsub("\\.+", "_", name)  # Replace multiple dots with underscores
  name
}

# Save results and plots
for (comparison in names(simplified_go_results)) {
  if (!is.null(simplified_go_results[[comparison]]) && nrow(as.data.frame(simplified_go_results[[comparison]])) > 0) {
    # Original column name for plot title
    plot_title <- paste("Simplified GO Enrichment for", comparison)
    
    # Sanitize column name for file names
    sanitized_name <- sanitize_name(comparison)
    
    # Save simplified GO results to CSV
    output_file <- file.path(output_folder, paste0(sanitized_name, "_Simplified_GO_results_2024-11-23.csv"))
    write.csv(as.data.frame(simplified_go_results[[comparison]]), file = output_file, row.names = FALSE)
    
    # Generate the dotplot
    dotplot_obj <- tryCatch({
      dotplot(simplified_go_results[[comparison]], showCategory = 20) +
        ggtitle(plot_title) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    }, error = function(e) {
      cat("Error in dotplot for", comparison, ":", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(dotplot_obj)) {
      # Save the plot as a PNG file
      plot_file <- file.path(output_folder, paste0(sanitized_name, "_Simplified_dotplot_2024-11-23.png"))
      ggsave(filename = plot_file, plot = dotplot_obj, width = 10, height = 7, dpi = 300)
    }
  } else {
    cat("No simplified GO results for", comparison, "\n")
  }
}
