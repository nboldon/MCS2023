# Load necessary libraries
library(ggplot2)
library(topGO)
library(goseq)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(clusterProfiler)

# Set your working directory where the CSV files are located
setwd("/Volumes/DataBox/C18_Subset/C18_Tx_Comps")  # Change this to your folder path

# List all CSV files in the directory
files <- list.files(pattern = "*.csv")

# Loop through each file and perform the analysis
for (file in files) {
  
  # Read the data
  sample_data <- read.csv(file)
  
  # Filter genes with abs(Log2FC) >= 0.5
  genes_to_test <- sample_data[abs(sample_data$Log2FC) >= 0.5, "name"]
  
  if (length(genes_to_test) == 0) {
    print(paste0("No significant genes found in ", file))
    next  # Skip this iteration if no genes pass the filter
  }
  
  ####################################
  # Gene Ontology Enrichment using clusterProfiler
  ####################################
  
  tryCatch({
    GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
    
    # Convert the results to a data frame
    GO_df <- as.data.frame(GO_results)
    
    if (nrow(GO_df) > 0) {
      # Create bar plot for top 20 GO terms
      fit <- plot(barplot(GO_results, showCategory = 20))
      
      # Save the plot as a PNG file, using the sample name from the file name
      png(paste0(sub(".csv", "", file), "_top20_GO.png"), res = 250, width = 1200, height = 2750)
      print(fit)
      dev.off()
    } else {
      print(paste0("No GO terms found for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in GO enrichment for ", file, ": ", e$message))
  })
  
  ####################################
  # TopGO Analysis
  ####################################
  
  tryCatch({
    # Prepare the gene list (a named vector of scores)
    geneList <- sample_data$Log2FC
    names(geneList) <- sample_data$name
    
    # Create topGOdata object
    GOdata <- new("topGOdata", 
                  ontology = "BP", 
                  allGenes = geneList, 
                  geneSel = function(x) abs(x) >= 0.5, 
                  nodeSize = 10, 
                  annot = annFUN.org, 
                  mapping = "org.Mm.eg.db", 
                  ID = "SYMBOL")
    
    # Run enrichment analysis
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    
    # Get top 10 GO terms
    top_results <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
    
    if (nrow(top_results) > 0) {
      # Plot the bar plot using ggplot2
      ggplot(top_results, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = "Top 10 GO Biological Process Terms", x = "GO Term", y = "-log10(Fisher p-value)") +
        theme_minimal() +
        ggsave(paste0(sub(".csv", "", file), "_TopGO_Top10_Barplot.png"))
    } else {
      print(paste0("No significant GO terms in TopGO for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in TopGO for ", file, ": ", e$message))
  })
  
  ####################################
  # GOseq Analysis
  ####################################
  
  tryCatch({
    # Create a binary vector indicating which genes are differentially expressed
    geneVector <- as.integer(abs(sample_data$Log2FC) >= 0.5)
    names(geneVector) <- sample_data$name
    
    # Perform GO analysis
    pwf <- nullp(geneVector, "mm10", "geneSymbol")  # Adjust for length bias
    GO_results_goseq <- goseq(pwf, "mm10", "geneSymbol")
    
    # Filter significant results
    significantGO <- GO_results_goseq[GO_results_goseq$over_represented_pvalue < 0.05, ]
    
    if (nrow(significantGO) > 0) {
      # Select top 10 terms
      top_results_goseq <- head(significantGO[order(significantGO$over_represented_pvalue), ], 10)
      
      # Plot bar plot of top 10 GO terms
      ggplot(top_results_goseq, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = "Top 10 GO Biological Process Terms (GOseq)", x = "GO Term", y = "-log10(p-value)") +
        theme_minimal() +
        ggsave(paste0(sub(".csv", "", file), "_GOseq_Top10_Barplot.png"))
    } else {
      print(paste0("No significant GO terms in GOseq for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in GOseq for ", file, ": ", e$message))
  })
  
  ####################################
  # BioMart GO Term Mapping
  ####################################
  
  tryCatch({
    # Connect to Ensembl
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # Retrieve GO terms for the genes
    go_data <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), 
                     filters = "external_gene_name", 
                     values = genes_to_test, 
                     mart = ensembl)
    
    if (nrow(go_data) > 0) {
      # Count the number of genes associated with each GO term
      go_counts <- as.data.frame(table(go_data$name_1006))
      
      # Plot top 10 GO terms by gene count
      top_go <- head(go_counts[order(go_counts$Freq, decreasing = TRUE), ], 10)
      
      ggplot(top_go, aes(x = reorder(Var1, Freq), y = Freq)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = "Top 10 GO Terms by Gene Count", x = "GO Term", y = "Gene Count") +
        theme_minimal() +
        ggsave(paste0(sub(".csv", "", file), "_BioMart_Top10_GO_Count.png"))
    } else {
      print(paste0("No GO terms found in BioMart for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in BioMart analysis for ", file, ": ", e$message))
  })
  
  # Print message to indicate completion for this sample
  print(paste0("Finished processing ", file))
}
