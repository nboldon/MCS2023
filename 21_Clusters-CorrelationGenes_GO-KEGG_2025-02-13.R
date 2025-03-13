

# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_Cluster_zscore-Stats")


# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(viridis)
#library(KEGG.db)
library(KEGGREST)

# Read the correlation data
corr_data <- read.csv("C3_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv")

# Filter for significant genes
sig_genes <- corr_data %>%
  filter(Significant == TRUE) %>%
  pull(Gene)

# Convert gene symbols to ENTREZ IDs
gene_ids <- bitr(sig_genes, 
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

# Run GO enrichment analysis for all three GO categories
go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

# Run KEGG pathway analysis
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,
                          organism = 'mmu',  # 'mmu' for mouse, 'hsa' for human
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# Save results to CSV files
write.csv(as.data.frame(go_bp), "GO_biological_process_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "GO_molecular_function_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "GO_cellular_component_results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_result), "KEGG_pathway_results.csv", row.names = FALSE)

# Function to create and save GO/KEGG plots
create_enrichment_plot <- function(enrichment_data, title, filename) {
  if (nrow(as.data.frame(enrichment_data)) > 0) {
    p <- barplot(enrichment_data,
                 showCategory = 15,  # Show top 15 terms
                 x = "Count",
                 title = title) +
      scale_fill_viridis() +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
    
    ggsave(filename, p, width = 12, height = 8)
  }
}

# Create and save plots
create_enrichment_plot(go_bp, "Biological Process GO Terms", "GO_biological_process_plot.pdf")
create_enrichment_plot(go_mf, "Molecular Function GO Terms", "GO_molecular_function_plot.pdf")
create_enrichment_plot(go_cc, "Cellular Component GO Terms", "GO_cellular_component_plot.pdf")
create_enrichment_plot(kegg_result, "KEGG Pathway Enrichment", "KEGG_pathway_plot.pdf")

# Create KEGG pathway visualization
# Only create if there are significant pathways
if (nrow(as.data.frame(kegg_result)) > 0) {
  # Create a dotplot for KEGG pathways
  kegg_dot <- dotplot(kegg_result, 
                      showCategory = 15,
                      title = "KEGG Pathway Enrichment") +
    scale_color_viridis() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
  
  ggsave("KEGG_pathway_dotplot.pdf", kegg_dot, width = 12, height = 8)
}

# Print summary of results
cat("Number of enriched BP terms:", nrow(as.data.frame(go_bp)), "\n")
cat("Number of enriched MF terms:", nrow(as.data.frame(go_mf)), "\n")
cat("Number of enriched CC terms:", nrow(as.data.frame(go_cc)), "\n")
cat("Number of enriched KEGG pathways:", nrow(as.data.frame(kegg_result)), "\n")



################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################



## Loop for all files with the same file ending in the directory


# Function to run enrichment analysis with data validation
run_enrichment_analysis <- function(file_path) {
  # Extract task name from filename for output organization
  task_name <- gsub(".*/(.*?)_CorrelationStats.*", "\\1", file_path)
  cat("\n===========================================")
  cat("\nProcessing task:", task_name, "\n")
  
  # Read the correlation data
  corr_data <- read.csv(file_path)
  
  # Print initial data summary
  cat("\nInitial data summary:")
  cat("\nTotal rows in file:", nrow(corr_data))
  cat("\nColumns present:", paste(colnames(corr_data), collapse=", "))
  cat("\nNumber of TRUE in Significant column:", sum(corr_data$Significant))
  
  # Filter for significant genes
  sig_genes <- corr_data %>%
    filter(Significant == TRUE) %>%
    pull(Gene)
  
  cat("\nNumber of significant genes found:", length(sig_genes))
  if(length(sig_genes) > 0) {
    cat("\nSignificant genes:", paste(sig_genes, collapse=", "))
  }
  
  # Skip if no significant genes
  if(length(sig_genes) == 0) {
    cat("\nNo significant genes found in", file_path, "\n")
    return()
  }
  
  # Create output directory for this task
  dir_name <- paste0("enrichment_results_", task_name)
  dir.create(dir_name, showWarnings = FALSE)
  
  # Convert gene symbols to ENTREZ IDs with error checking
  gene_ids <- tryCatch({
    bitr(sig_genes, 
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)
  }, error = function(e) {
    cat("\nError in gene ID conversion:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if(is.null(gene_ids)) {
    cat("\nGene ID conversion failed\n")
    return()
  }
  
  cat("\nNumber of genes successfully mapped to ENTREZ IDs:", nrow(gene_ids))
  unmapped_genes <- setdiff(sig_genes, gene_ids$SYMBOL)
  if(length(unmapped_genes) > 0) {
    cat("\nGenes not mapped to ENTREZ IDs:", paste(unmapped_genes, collapse=", "))
  }
  
  # Run GO enrichment analyses
  go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
  
  go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
  
  go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
  
  # Test KEGG separately before full analysis
  cat("\nTesting KEGG database connection...")
  kegg_test <- tryCatch({
    enrichKEGG(gene = head(gene_ids$ENTREZID, n=10),
               organism = 'mmu',
               pvalueCutoff = 1,  # Set high for testing
               qvalueCutoff = 1)
  }, error = function(e) {
    cat("\nKEGG test error:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if(!is.null(kegg_test)) {
    cat("\nKEGG database connection successful")
  }
  
  # Run KEGG pathway analysis
  cat("\nRunning full KEGG analysis...\n")
  kegg_result <- tryCatch({
    enrichKEGG(gene = gene_ids$ENTREZID,
               organism = 'mmu',
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)
  }, error = function(e) {
    cat("KEGG analysis error:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Add explicit KEGG result checking
  cat("\nRunning full KEGG analysis...\n")
  kegg_result <- tryCatch({
    enrichKEGG(gene = gene_ids$ENTREZID,
               organism = 'mmu',
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)
  }, error = function(e) {
    cat("KEGG analysis error:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Detailed KEGG results inspection
  if(!is.null(kegg_result)) {
    cat("KEGG analysis completed\n")
    kegg_df <- as.data.frame(kegg_result)
    cat("Number of total KEGG pathways tested:", length(kegg_result@universe), "\n")
    cat("Number of genes used in KEGG analysis:", length(unique(gene_ids$ENTREZID)), "\n")
    
    if(nrow(kegg_df) > 0) {
      cat("Number of enriched KEGG pathways:", nrow(kegg_df), "\n")
      cat("Top enriched pathways:\n")
      print(head(kegg_df[, c("ID", "Description", "pvalue", "Count")]))
    } else {
      # Try running with relaxed p-value to see if there are any hits at all
      kegg_relaxed <- enrichKEGG(gene = gene_ids$ENTREZID,
                                 organism = 'mmu',
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1)
      cat("\nNo significant KEGG pathways at p<0.05\n")
      cat("Number of pathways with any hits (p<1):", nrow(as.data.frame(kegg_relaxed)), "\n")
      if(nrow(as.data.frame(kegg_relaxed)) > 0) {
        cat("Top pathways (even if not significant):\n")
        print(head(as.data.frame(kegg_relaxed)[, c("ID", "Description", "pvalue", "Count")]))
      }
    }
  } else {
    cat("KEGG analysis failed\n")
  }
  
  # Save KEGG results even if not significant
  if(!is.null(kegg_result)) {
    kegg_df <- as.data.frame(kegg_result)
    write.csv(kegg_df, 
              file.path(dir_name, "KEGG_pathway_results.csv"), 
              row.names = FALSE)
    
    # Also save relaxed results for reference
    kegg_relaxed <- enrichKEGG(gene = gene_ids$ENTREZID,
                               organism = 'mmu',
                               pvalueCutoff = 1,
                               qvalueCutoff = 1)
    write.csv(as.data.frame(kegg_relaxed), 
              file.path(dir_name, "KEGG_pathway_results_all.csv"), 
              row.names = FALSE)
  }
  
  # Save results to CSV files
  write.csv(as.data.frame(go_bp), 
            file.path(dir_name, "GO_biological_process_results.csv"), 
            row.names = FALSE)
  write.csv(as.data.frame(go_mf), 
            file.path(dir_name, "GO_molecular_function_results.csv"), 
            row.names = FALSE)
  write.csv(as.data.frame(go_cc), 
            file.path(dir_name, "GO_cellular_component_results.csv"), 
            row.names = FALSE)
  
  # Save KEGG results if they exist
  if(!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    write.csv(as.data.frame(kegg_result), 
              file.path(dir_name, "KEGG_pathway_results.csv"), 
              row.names = FALSE)
  }
  
  # Function to create and save enrichment plots
  create_enrichment_plot <- function(enrichment_data, title, filename) {
    if(nrow(as.data.frame(enrichment_data)) > 0) {
      p <- barplot(enrichment_data,
                   showCategory = 20,
                   x = "Count",
                   title = paste(task_name, "-", title)) +
        scale_fill_viridis() +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10))
      
      ggsave(file.path(dir_name, filename), p, width = 12, height = 8)
    }
  }
  
  # Create and save GO plots
  create_enrichment_plot(go_bp, "Biological Process GO Terms", "GO_biological_process_plot.pdf")
  create_enrichment_plot(go_mf, "Molecular Function GO Terms", "GO_molecular_function_plot.pdf")
  create_enrichment_plot(go_cc, "Cellular Component GO Terms", "GO_cellular_component_plot.pdf")
  
  # Create KEGG plots if results exist
  if(!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
    kegg_bar <- barplot(kegg_result,
                        showCategory = 20,
                        title = paste(task_name, "- KEGG Pathway Enrichment")) +
      scale_fill_viridis() +
      theme_minimal()
    
    kegg_dot <- dotplot(kegg_result, 
                        showCategory = 20,
                        title = paste(task_name, "- KEGG Pathway Enrichment")) +
      scale_color_viridis() +
      theme_minimal()
    
    ggsave(file.path(dir_name, "KEGG_pathway_barplot.pdf"), kegg_bar, width = 12, height = 8)
    ggsave(file.path(dir_name, "KEGG_pathway_dotplot.pdf"), kegg_dot, width = 12, height = 8)
  }
  
  # Save the gene ID mapping for reference
  write.csv(gene_ids, 
            file.path(dir_name, "gene_id_mapping.csv"), 
            row.names = FALSE)
  
  # Create summary file
  sink(file.path(dir_name, "analysis_summary.txt"))
  cat("Analysis Summary for", task_name, "\n")
  cat("Number of significant genes:", length(sig_genes), "\n")
  cat("Number of genes mapped to ENTREZ IDs:", nrow(gene_ids), "\n")
  cat("Number of enriched BP terms:", nrow(as.data.frame(go_bp)), "\n")
  cat("Number of enriched MF terms:", nrow(as.data.frame(go_mf)), "\n")
  cat("Number of enriched CC terms:", nrow(as.data.frame(go_cc)), "\n")
  if(!is.null(kegg_result)) {
    cat("Number of enriched KEGG pathways:", nrow(as.data.frame(kegg_result)), "\n")
  }
  sink()
  
  cat("\nAnalysis completed for", task_name, "\n")
}

# Get list of files and process them
files <- list.files(pattern = "*CorrelationStats.*withFDR_2025-02-12\\.csv$", full.names = TRUE)
cat("Found", length(files), "files to process\n")

# Create a directory for all results
main_results_dir <- "enrichment_results_all"
dir.create(main_results_dir, showWarnings = FALSE)

# Process all files
for(file in files) {
  tryCatch({
    run_enrichment_analysis(file)
  }, error = function(e) {
    cat("\nError processing", file, ":", conditionMessage(e), "\n")
  })
}

cat("\n===========================================")
cat("\nProcessing complete!")
cat("\n===========================================\n")
