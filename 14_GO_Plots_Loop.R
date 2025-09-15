# Load necessary libraries

library(ggplot2)
library(topGO)
library(goseq)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(clusterProfiler)



## Similar code to below, but includes saving a .csv file of Cluster Profiler results and today's date (1/23/25)



# Set your working directory where the CSV files are located
setwd("/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps") 

# List all CSV files in the directory
files <- list.files(pattern = "*.csv")

# Get today's date
today <- Sys.Date()

# Loop through each file and perform the analysis
for (file in files) {
  
  # Read the data
  sample_data <- read.csv(file)
  
  # Extract the treatment comparison from the file name (modify this if your filename structure is different)
  treatment_comparison <- sub(".csv", "", file)
  
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
      fit <- plot(barplot(GO_results, showCategory = 20, title = paste("Top 20 GO Terms -", treatment_comparison)))
      
      # Save the plot as a PNG file, using the sample name from the file name
      png(paste0(treatment_comparison, "_top20_GO.png"), res = 250, width = 1600, height = 3000)
      print(fit)
      dev.off()
      
      # Save the GO results to a CSV file with today's date
      write.csv(GO_df, paste0(treatment_comparison, "_GO_results_", today, ".csv"), row.names = FALSE)
      
    } else {
      print(paste0("No GO terms found for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in GO enrichment for ", file, ": ", e$message))
  })
}



####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################

## Results summary:

[1] "No significant genes found in C1_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C1_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C1_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C1_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C1_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C10_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C10_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"

[1] "No GO terms found for C10_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C10_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C11_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C11_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C11_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No GO terms found for C11_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C11_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C11_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C11_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C12_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C12_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C12_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C12_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C12_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C12_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C12_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C12_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C12_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C13_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C13_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C13_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C13_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C13_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C13_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C13_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C14_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C14_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
No gene sets have size between 10 and 500 ...
--> return NULL...
[1] "No GO terms found for C14_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C14_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C14_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C14_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C15_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C15_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C15_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C15_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C15_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C16_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C16_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C16_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C16_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C16_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C16_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C16_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No GO terms found for C16_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C17_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C17_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C17_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C17_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
No gene sets have size between 10 and 500 ...
--> return NULL...
[1] "No GO terms found for C17_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C17_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv"
[1] "No GO terms found for C18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No GO terms found for C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C18_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C19_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C19_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C2_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C2_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C2_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C2_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C2_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C2_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C2_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C2_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
No gene sets have size between 10 and 500 ...
--> return NULL...
[1] "No GO terms found for C20_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C20_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C20_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C20_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C20_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C20_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C21_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C21_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C21_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C21_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
No gene sets have size between 10 and 500 ...
--> return NULL...
[1] "No GO terms found for C21_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C22_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C22_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C22_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C22_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C22_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C22_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C22_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C23_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C23_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C23_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C23_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C23_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C24_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C24_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C24_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C24_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C24_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C24_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C24_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C25_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C25_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C3_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C3_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C3_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No GO terms found for C3_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C3_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C4_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C4_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
No gene sets have size between 10 and 500 ...
--> return NULL...
[1] "No GO terms found for C4_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C4_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C4_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C4_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C5_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C5_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C5_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C5_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C5_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C5_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C5_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C6_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C6_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No GO terms found for C6_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No GO terms found for C6_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C6_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C6_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C6_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C7_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C7_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C7_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C7_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C7_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C7_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C7_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C8_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C8_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C8_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No GO terms found for C8_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C8_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C8_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"
[1] "No significant genes found in C9_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"
[1] "No significant genes found in C9_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv"
[1] "No significant genes found in C9_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C9_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C9_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"
[1] "No significant genes found in C9_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
[1] "No significant genes found in C9_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"


##################################################################
##################################################################
##################################################################


## To combine spreadsheets from loop:


# Load necessary library
library(dplyr)

# Set the directory where the result files are stored
output_dir <- "/Volumes/DataBox/GO_Analysis/Clusters"

# List all the .csv files in the output directory, excluding temporary files
csv_files <- list.files(output_dir, pattern = "GO_results_2025-01-23.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("^~\\$", basename(csv_files))]  # Remove temporary files starting with ~$

# Create an empty list to store the data frames
all_data <- list()

# Loop through each .csv file and read it into a data frame
for (file in csv_files) {
  # Try to read the file with error handling
  tryCatch({
    # Read the data from the CSV file
    data <- read.csv(file)
    
    # Check if data was actually read
    if(nrow(data) > 0) {
      # Ensure 'geneID' is always treated as a character (if this column exists)
      if("geneID" %in% colnames(data)) {
        data$geneID <- as.character(data$geneID)
      }
      
      # Assign the full file name to the new file_name column
      data$file_name <- basename(file)
      
      # Append the data frame to the list
      all_data <- append(all_data, list(data))
      
      # Debug: Check the number of rows being added
      message(paste("Rows in file", basename(file), ":", nrow(data)))
    } else {
      message(paste("Warning: No data in file", basename(file)))
    }
  }, error = function(e) {
    message(paste("Error processing file", basename(file), ":", e$message))
  })
}

# Combine all the data frames into one, but first check if any data was loaded
if(length(all_data) > 0) {
  library(dplyr)
  combined_data <- bind_rows(all_data)
  
  # Debug: Check the number of rows after combining
  message(paste("Total rows after combining:", nrow(combined_data)))
  
  # Optionally, save the combined data to a new CSV file
  output_file <- "/Volumes/DataBox/GO_Analysis/Clusters/GO_Combined_Cluster_Results.csv"
  write.csv(combined_data, output_file, row.names = FALSE)
  
  # Check if the file is saved successfully
  if(file.exists(output_file)) {
    message("The combined file has been saved to: ", output_file)
  } else {
    message("Failed to save the combined file.")
  }
  
  # Check the first few rows of the combined data
  print(head(combined_data))
} else {
  message("No data was successfully loaded from any of the files.")
}






##################################################################
##################################################################
##################################################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################


1/23/2025

sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Denver
tzcode source: internal

attached base packages:
  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] clusterProfiler_4.12.6                    biomaRt_2.60.1                           
[3] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0 GenomicFeatures_1.56.0                   
[5] GenomicRanges_1.56.1                      GenomeInfoDb_1.40.1                      
[7] goseq_1.56.0                              geneLenDataBase_1.40.1                   
[9] BiasedUrn_2.0.12                          topGO_2.56.0                             
[11] SparseM_1.84-2                            GO.db_3.19.1                             
[13] AnnotationDbi_1.66.0                      IRanges_2.38.1                           
[15] S4Vectors_0.42.1                          Biobase_2.64.0                           
[17] graph_1.82.0                              BiocGenerics_0.50.0                      
[19] ggplot2_3.5.1                            

loaded via a namespace (and not attached):
  [1] splines_4.4.0               BiocIO_1.14.0               bitops_1.0-8               
[4] ggplotify_0.1.2             filelock_1.0.3              tibble_3.2.1               
[7] R.oo_1.26.0                 polyclip_1.10-7             XML_3.99-0.17              
[10] lifecycle_1.0.4             httr2_1.0.3                 doParallel_1.0.17          
[13] lattice_0.22-6              MASS_7.3-61                 magrittr_2.0.3             
[16] rmarkdown_2.28              yaml_2.3.10                 cowplot_1.1.3              
[19] DBI_1.2.3                   RColorBrewer_1.1-3          abind_1.4-5                
[22] zlibbioc_1.50.0             purrr_1.0.2                 R.utils_2.12.3             
[25] ggraph_2.2.1                RCurl_1.98-1.16             yulab.utils_0.1.7          
[28] tweenr_2.0.3                rappdirs_0.3.3              circlize_0.4.16            
[31] GenomeInfoDbData_1.2.12     enrichplot_1.24.4           ggrepel_0.9.6              
[34] tidytree_0.4.6              pheatmap_1.0.12             codetools_0.2-20           
[37] DelayedArray_0.30.1         DOSE_3.30.5                 xml2_1.3.6                 
[40] ggforce_0.4.2               tidyselect_1.2.1            shape_1.4.6.1              
[43] aplot_0.2.3                 UCSC.utils_1.0.0            farver_2.1.2               
[46] viridis_0.6.5               matrixStats_1.4.1           BiocFileCache_2.12.0       
[49] GenomicAlignments_1.40.0    jsonlite_1.8.8              GetoptLong_1.0.5           
[52] tidygraph_1.3.1             iterators_1.0.14            foreach_1.5.2              
[55] tools_4.4.0                 progress_1.2.3              treeio_1.28.0              
[58] Rcpp_1.0.13                 glue_1.7.0                  gridExtra_2.3              
[61] SparseArray_1.4.8           xfun_0.47                   mgcv_1.9-1                 
[64] qvalue_2.36.0               MatrixGenerics_1.16.0       dplyr_1.1.4                
[67] withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                
[70] digest_0.6.37               R6_2.5.1                    gridGraphics_0.5-1         
[73] colorspace_2.1-1            RSQLite_2.3.7               R.methodsS3_1.8.2          
[76] utf8_1.2.4                  tidyr_1.3.1                 generics_0.1.3             
[79] data.table_1.16.0           rtracklayer_1.64.0          prettyunits_1.2.0          
[82] graphlayouts_1.1.1          httr_1.4.7                  S4Arrays_1.4.1             
[85] scatterpie_0.2.4            org.Mm.eg.db_3.19.1         pkgconfig_2.0.3            
[88] gtable_0.3.5                blob_1.2.4                  ComplexHeatmap_2.20.0      
[91] XVector_0.44.0              shadowtext_0.1.4            htmltools_0.5.8.1          
[94] fgsea_1.30.0                clue_0.3-65                 scales_1.3.0               
[97] png_0.1-8                   ggfun_0.1.6                 knitr_1.48                 
[100] rstudioapi_0.16.0           reshape2_1.4.4              rjson_0.2.22               
[103] nlme_3.1-166                curl_5.2.2                  cachem_1.1.0               
[106] GlobalOptions_0.1.2         stringr_1.5.1               parallel_4.4.0             
[109] restfulr_0.0.15             pillar_1.9.0                grid_4.4.0                 
[112] vctrs_0.6.5                 dbplyr_2.5.0                cluster_2.1.6              
[115] evaluate_0.24.0             cli_3.6.3                   compiler_4.4.0             
[118] Rsamtools_2.20.0            rlang_1.1.4                 crayon_1.5.3               
[121] labeling_0.4.3              plyr_1.8.9                  fs_1.6.4                   
[124] stringi_1.8.4               viridisLite_0.4.2           BiocParallel_1.38.0        
[127] txdbmaker_1.0.1             munsell_0.5.1               Biostrings_2.72.1          
[130] lazyeval_0.2.2              pacman_0.5.1                GOSemSim_2.30.2            
[133] Matrix_1.7-0                hms_1.1.3                   patchwork_1.2.0            
[136] bit64_4.0.5                 KEGGREST_1.44.1             SummarizedExperiment_1.34.0
[139] igraph_2.0.3                memoise_2.0.1               ggtree_3.12.0              
[142] fastmatch_1.1-4             bit_4.0.5                   ape_5.8                    
[145] gson_0.1.0          




####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################


## Loop to run multiple different GO functions:
# List all CSV files in the directory
files <- list.files(pattern = "*.csv")

# Loop through each file and perform the analysis
for (file in files) {
  
  # Read the data
  sample_data <- read.csv(file)
  
  # Extract the treatment comparison from the file name (modify this if your filename structure is different)
  treatment_comparison <- sub(".csv", "", file)
  
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
      fit <- plot(barplot(GO_results, showCategory = 20, title = paste("Top 20 GO Terms -", treatment_comparison)))
      
      # Save the plot as a PNG file, using the sample name from the file name
      png(paste0(treatment_comparison, "_top20_GO.png"), res = 250, width = 1600, height = 3000)
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
                  nodeSize = 20, 
                  annot = annFUN.org, 
                  mapping = "org.Mm.eg.db", 
                  ID = "SYMBOL")
    
    # Run enrichment analysis
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    
    # Get top 20 GO terms
    top_results <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
    
    if (nrow(top_results) > 0) {
      # Plot the bar plot using ggplot2
      ggplot(top_results, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("Top 20 GO Biological Process Terms -", treatment_comparison), x = "GO Term", y = "-log10(Fisher p-value)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
        ggsave(paste0(treatment_comparison, "_TopGO_Top20_Barplot.png"), width = 10, height = 7)
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
    # Create a binary vector indicating which genes are differentially accessible
    geneVector <- as.integer(abs(sample_data$Log2FC) >= 0.5)
    names(geneVector) <- sample_data$name
    
    # Perform GO analysis
    pwf <- nullp(geneVector, "mm10", "geneSymbol")  # Adjust for length bias
    GO_results_goseq <- goseq(pwf, "mm10", "geneSymbol")
    
    # Filter significant results
    significantGO <- GO_results_goseq[GO_results_goseq$over_represented_pvalue < 0.05, ]
    
    if (nrow(significantGO) > 0) {
      # Select top 20 terms
      top_results_goseq <- head(significantGO[order(significantGO$over_represented_pvalue), ], 20)
      
      # Plot bar plot of top 20 GO terms
      ggplot(top_results_goseq, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("Top 20 GO Biological Process Terms (GOseq) -", treatment_comparison), x = "GO Term", y = "-log10(p-value)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
        ggsave(paste0(treatment_comparison, "_GOseq_Top20_Barplot.png"), width = 10, height = 7)
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
      
      # Plot top 20 GO terms by gene count
      top_go <- head(go_counts[order(go_counts$Freq, decreasing = TRUE), ], 20)
      
      ggplot(top_go, aes(x = reorder(Var1, Freq), y = Freq)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste("Top 20 GO Terms by Gene Count -", treatment_comparison), x = "GO Term", y = "Gene Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
        ggsave(paste0(treatment_comparison, "_BioMart_Top20_GO_Count.png"), width = 10, height = 7)
    } else {
      print(paste0("No GO terms found in BioMart for ", file))
    }
  }, error = function(e) {
    print(paste0("Error in BioMart analysis for ", file, ": ", e$message))
  })
  
  # Print message to indicate completion for this sample
  print(paste0("Finished processing ", file))
}
