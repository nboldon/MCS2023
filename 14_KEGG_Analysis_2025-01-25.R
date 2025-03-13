# Load necessary libraries
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)


setwd("/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps") 



##################################################################



## Run for 1 file only:



# Load the data
file <- "C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"

# Check if the file exists and is not empty
if (!file.exists(file)) {
  stop(paste("The file", file, "does not exist. Please check the path."))
}

data <- read.csv(file)
if (nrow(data) == 0) {
  stop(paste("The file", file, "is empty. Please check the input data."))
}

# Extract upregulated and downregulated genes based on Log2FC
data_up <- dplyr::filter(data, Log2FC > 0)   # Upregulated genes
data_down <- dplyr::filter(data, Log2FC < 0) # Downregulated genes

# Extract gene names and clean up missing or invalid entries
genes_up <- data_up$name[!is.na(data_up$name) & data_up$name != ""]
genes_down <- data_down$name[!is.na(data_down$name) & data_down$name != ""]

# Check if there are any genes to process
if (length(genes_up) == 0 && length(genes_down) == 0) {
  stop("No upregulated or downregulated genes found in the dataset.")
}

# Convert gene names to Entrez IDs
entrez_up <- if (length(genes_up) > 0) {
  na.omit(mapIds(org.Mm.eg.db, genes_up, "ENTREZID", "SYMBOL"))
} else {
  NULL
}

entrez_down <- if (length(genes_down) > 0) {
  na.omit(mapIds(org.Mm.eg.db, genes_down, "ENTREZID", "SYMBOL"))
} else {
  NULL
}

# Perform KEGG pathway analysis for upregulated and downregulated genes
kegg_up <- if (!is.null(entrez_up) && length(entrez_up) > 0) {
  enrichKEGG(
    gene = entrez_up,
    organism = "mmu",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else {
  warning("No valid upregulated genes for KEGG analysis.")
  NULL
}

kegg_down <- if (!is.null(entrez_down) && length(entrez_down) > 0) {
  enrichKEGG(
    gene = entrez_down,
    organism = "mmu",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else {
  warning("No valid downregulated genes for KEGG analysis.")
  NULL
}

# Save KEGG results to CSV files
if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
  write.csv(as.data.frame(kegg_up), file = "C3_T1vsT4Comp_KEGG_Upregulated.csv", row.names = FALSE)
  message("KEGG results for upregulated genes saved to 'C3_T1vsT4Comp_KEGG_Upregulated.csv'")
}

if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
  write.csv(as.data.frame(kegg_down), file = "C3_T1vsT4Comp_KEGG_Downregulated.csv", row.names = FALSE)
  message("KEGG results for downregulated genes saved to 'C3_T1vsT4Comp_KEGG_Downregulated.csv'")
}

# Generate KEGG Cytoscape-style network plots and save to a PDF
pdf_file <- paste0("C3_T1vsT4Comp_Cytoscape_Network_", Sys.Date(), ".pdf")
pdf_opened <- FALSE

if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
  pdf(pdf_file)
  pdf_opened <- TRUE
  cnetplot(kegg_up, showCategory = 5, title = "Cytoscape Upregulated")
} else {
  message("No KEGG results for upregulated genes.")
}

if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
  if (!pdf_opened) {
    pdf(pdf_file)
    pdf_opened <- TRUE
  }
  cnetplot(kegg_down, showCategory = 5, title = "Cytoscape Downregulated")
} else {
  message("No KEGG results for downregulated genes.")
}

# Close the PDF device if it was opened
if (pdf_opened) {
  dev.off()
  message(paste("PDF saved as", pdf_file))
} else {
  message("No plots were generated, so no PDF file was created.")
}

# Summary messages
if (!is.null(entrez_up)) {
  message(paste("Number of upregulated genes mapped to Entrez IDs:", length(entrez_up)))
}
if (!is.null(entrez_down)) {
  message(paste("Number of downregulated genes mapped to Entrez IDs:", length(entrez_down)))
}


print(head(genes_up))
print(entrez_up)
print(head(genes_down))
print(entrez_down)





##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################


## To loop through all files in directory:


# Load necessary libraries
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(RCy3)

# Set up the Ensembl mart for mouse (mm10)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Set the directory containing the .csv files
input_dir <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps"
output_dir <- "/Volumes/DataBox/KEGG_Cytoscape"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# List all .csv files in the input directory
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Function to convert gene symbols to Entrez IDs using biomaRt
convert_gene_symbols_to_entrez <- function(gene_symbols) {
  conversion <- getBM(attributes = c("external_gene_name", "entrezgene_id"), 
                      filters = "external_gene_name", 
                      values = gene_symbols, 
                      mart = ensembl)
  
  # Return the Entrez IDs as a named vector
  return(setNames(conversion$entrezgene_id, conversion$external_gene_name))
}

# Loop through each .csv file
for (file in csv_files) {
  # Extract the base name of the file (without extension) for naming outputs
  file_base <- tools::file_path_sans_ext(basename(file))
  
  # Load the data and handle potential errors
  tryCatch({
    message(paste("Processing file:", file))
    data <- read.csv(file)
    
    if (nrow(data) == 0) {
      warning(paste("The file", file, "is empty. Skipping."))
      next
    }
    
    # Extract upregulated and downregulated genes based on Log2FC
    data_up <- dplyr::filter(data, Log2FC > 0)   # Upregulated genes
    data_down <- dplyr::filter(data, Log2FC < 0) # Downregulated genes
    
    # Extract gene symbols and clean up missing or invalid entries
    genes_up <- data_up$name[!is.na(data_up$name) & data_up$name != ""]
    genes_down <- data_down$name[!is.na(data_down$name) & data_down$name != ""]
    
    if (length(genes_up) == 0 && length(genes_down) == 0) {
      warning(paste("No upregulated or downregulated genes found in", file, ". Skipping."))
      next
    }
    
    # Convert gene symbols to Entrez IDs using biomaRt
    entrez_up <- if (length(genes_up) > 0) {
      convert_gene_symbols_to_entrez(genes_up)
    } else {
      NULL
    }
    
    entrez_down <- if (length(genes_down) > 0) {
      convert_gene_symbols_to_entrez(genes_down)
    } else {
      NULL
    }
    
    # Perform KEGG pathway analysis
    kegg_up <- if (!is.null(entrez_up) && length(entrez_up) > 0) {
      enrichKEGG(
        gene = entrez_up,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    } else {
      warning(paste("No valid upregulated genes for KEGG analysis in", file))
      NULL
    }
    
    kegg_down <- if (!is.null(entrez_down) && length(entrez_down) > 0) {
      enrichKEGG(
        gene = entrez_down,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    } else {
      warning(paste("No valid downregulated genes for KEGG analysis in", file))
      NULL
    }
    
    # Save KEGG results to CSV files
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      write.csv(as.data.frame(kegg_up), file = file.path(output_dir, paste0(file_base, "_KEGG_Upregulated.csv")), row.names = FALSE)
      message(paste("KEGG results for upregulated genes saved for", file))
    }
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      write.csv(as.data.frame(kegg_down), file = file.path(output_dir, paste0(file_base, "_KEGG_Downregulated.csv")), row.names = FALSE)
      message(paste("KEGG results for downregulated genes saved for", file))
    }
    
    # Generate Cytoscape network and export using Rcy3
    if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
      # Convert KEGG results to a data frame for network creation
      up_network <- data.frame(
        from = kegg_up$ID,  # KEGG pathway IDs
        to = rep("Upregulated", length(kegg_up$ID))  # Assigning upregulated gene category
      )
      createNetworkFromDataFrames(up_network, title = paste("Upregulated KEGG Network -", file_base), collection = "KEGG Networks")
    }
    
    if (!is.null(kegg_down) && nrow(as.data.frame(kegg_down)) > 0) {
      # Convert KEGG results to a data frame for network creation
      down_network <- data.frame(
        from = kegg_down$ID,  # KEGG pathway IDs
        to = rep("Downregulated", length(kegg_down$ID))  # Assigning downregulated gene category
      )
      createNetworkFromDataFrames(down_network, title = paste("Downregulated KEGG Network -", file_base), collection = "KEGG Networks")
    }
    
    # Export the network to Cytoscape
    exportNetworkToCytoscape(file.path(output_dir, paste0(file_base, "_Cytoscape_Network_", Sys.Date(), ".xgmml")))
    message(paste("Network for", file, "exported to Cytoscape"))
    
  }, error = function(e) {
    # Catch any errors during the processing of each file
    message(paste("Error processing file", file, ":", e$message))
  })
}



##################################################################
##################################################################
##################################################################


## To combine spreadsheets from loop:


# Load necessary library
library(dplyr)

# Set the directory where the result files are stored
output_dir <- "/Volumes/DataBox/KEGG_Cytoscape/KEGG_Output"

# List all the .csv files in the output directory
csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)

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
output_file <- "/Volumes/DataBox/KEGG_Cytoscape/KEGG_Combined_Results.csv"
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
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################



## Loop return:


Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C1_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Reading KEGG annotation online: "https://rest.kegg.jp/link/mmu/pathway"...
Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/mmu"...
--> No gene can be mapped....
--> Expected input gene ID: 11529,394430,14377,79459,226413,94215
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 18642,171282,66945,72094,14979,20916
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: 353204,11997,30963,56421,56348,230163
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
--> No gene can be mapped....
--> Expected input gene ID: 78070,269951,171282,270198,331026,18534
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 11370,18293,319625,67680,18534,243085
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
--> No gene can be mapped....
--> Expected input gene ID: 14751,72141,11409,71147,74559,11532
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C10_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
--> No gene can be mapped....
--> Expected input gene ID: 14380,110460,230163,231396,71832,100559
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
--> No gene can be mapped....
--> Expected input gene ID: 14121,72094,14194,11671,26922,18642
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 15926,11428,171281,110695,20917,72094
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C11_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C12_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C13_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 80911,29858,109729085,17330,53418,19895
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 319625,100040843,75778,225913,66646,26876
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
--> No gene can be mapped....
--> Expected input gene ID: 110695,30963,71519,69123,20322,11529
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C14_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 110460,72141,19895,18655,12183,16833
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 243996,74559,18293,270076,18597,14194
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C15_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 11997,18598,15926,74551,52538,17449
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C16_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 26897,18648,110119,218138,13117,74637
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C17_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 171282,69080,18563,433256,14187,14377
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 226413,394430,104112,75540,112417,54397
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
--> No gene can be mapped....
--> Expected input gene ID: 103988,21881,110446,171281,54397,394430
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-30.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C18_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 66198,16770,79459,328099,394434,29858
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C19_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C2_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 21991,14187,66775,14447,15929,230163
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 26922,115487111,70757,14194,56348,319801
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 67689,216739,235582,11997,11529,115487111
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C20_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 232714,99035,14433,54325,13382,319801
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 216019,12895,11364,277753,56421,74551
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 16832,52538,170718,66198,14377,100040843
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C21_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 75578,16828,235674,12686,13177,74637
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 67880,232714,218138,435802,72157,56727
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 11669,20917,80911,74419,11605,14194
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C22_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 54325,76051,71336,68801,353204,319625
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 16591,319625,277753,11669,56012,83553
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 66775,435802,74147,11429,54128,13807
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C23_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 16832,270076,56012,14595,69983,23986
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 11997,110460,18746,231396,110208,12895
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
--> No gene can be mapped....
--> Expected input gene ID: 30963,12686,216019,106529,269951,231396
--> return NULL...
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C24_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C25_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 75778,268756,18639,23986,328845,230163
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT4Comp_KEGG_Upregulated.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T1vsT4Comp_KEGG_Upregulated.csv : In argument: `Log2FC > 0`.
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
--> No gene can be mapped....
--> Expected input gene ID: 106529,353204,230163,14194,103988,328845
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 74147,11522,223722,18655,71773,234309
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 50790,268756,54325,71832,56012,218138
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 11430,230639,18746,69080,234309,68263
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C5_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
--> No gene can be mapped....
--> Expected input gene ID: 14378,19139,72535,14433,14120,231396
--> return NULL...
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 69080,14380,110208,66945,94284,16591
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
--> No gene can be mapped....
--> Expected input gene ID: 13382,99035,60525,50790,74559,26897
--> return NULL...
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C6_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C7_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for upregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
KEGG results for downregulated genes saved for /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : undefined columns selected
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Error processing file /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv : could not find function "exportNetworkToCytoscape"
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C8_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
Processing file: /Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C9_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
There were 50 or more warnings (use warnings() to see the first 50)
> 
  > sessionInfo()
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
  [1] RCy3_2.24.0            DOSE_3.30.5            clusterProfiler_4.12.6 org.Mm.eg.db_3.19.1    AnnotationDbi_1.66.0  
[6] IRanges_2.38.1         S4Vectors_0.42.1       Biobase_2.64.0         BiocGenerics_0.50.0    biomaRt_2.60.1        
[11] BiocManager_1.30.25   

loaded via a namespace (and not attached):
  [1] splines_4.4.0           pbdZMQ_0.3-13           bitops_1.0-9            ggplotify_0.1.2        
[5] filelock_1.0.3          tibble_3.2.1            R.oo_1.27.0             polyclip_1.10-7        
[9] graph_1.82.0            XML_3.99-0.18           lifecycle_1.0.4         httr2_1.0.7            
[13] doParallel_1.0.17       globals_0.16.3          lattice_0.22-6          MASS_7.3-64            
[17] backports_1.5.0         magrittr_2.0.3          rmarkdown_2.29          yaml_2.3.10            
[21] cowplot_1.1.3           DBI_1.2.3               RColorBrewer_1.1-3      zlibbioc_1.50.0        
[25] purrr_1.0.2             R.utils_2.12.3          RCurl_1.98-1.16         ggraph_2.2.1           
[29] yulab.utils_0.1.9       tweenr_2.0.3            rappdirs_0.3.3          circlize_0.4.16        
[33] GenomeInfoDbData_1.2.12 enrichplot_1.24.4       ggrepel_0.9.6           listenv_0.9.1          
[37] tidytree_0.4.6          pheatmap_1.0.12         parallelly_1.41.0       codetools_0.2-20       
[41] xml2_1.3.6              ggforce_0.4.2           tidyselect_1.2.1        shape_1.4.6.1          
[45] aplot_0.2.4             UCSC.utils_1.0.0        farver_2.1.2            viridis_0.6.5          
[49] base64enc_0.1-3         matrixStats_1.5.0       BiocFileCache_2.12.0    jsonlite_1.8.9         
[53] GetoptLong_1.0.5        tidygraph_1.3.1         iterators_1.0.14        foreach_1.5.2          
[57] tools_4.4.0             progress_1.2.3          treeio_1.28.0           Rcpp_1.0.14            
[61] glue_1.8.0              gridExtra_2.3           xfun_0.50               qvalue_2.36.0          
[65] IRdisplay_1.1           GenomeInfoDb_1.40.1     dplyr_1.1.4             withr_3.0.2            
[69] fastmap_1.2.0           fansi_1.0.6             caTools_1.18.3          digest_0.6.37          
[73] R6_2.5.1                gridGraphics_0.5-1      colorspace_2.1-1        GO.db_3.19.1           
[77] gtools_3.9.5            RSQLite_2.3.9           R.methodsS3_1.8.2       utf8_1.2.4             
[81] tidyr_1.3.1             generics_0.1.3          data.table_1.16.4       prettyunits_1.2.0      
[85] graphlayouts_1.2.1      httr_1.4.7              scatterpie_0.2.4        RJSONIO_1.3-1.9        
[89] pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4              ComplexHeatmap_2.20.0  
[93] XVector_0.44.0          shadowtext_0.1.4        htmltools_0.5.8.1       base64url_1.4          
[97] fgsea_1.30.0            clue_0.3-66             scales_1.3.0            png_0.1-8              
[101] ggfun_0.1.8             knitr_1.49              rstudioapi_0.17.1       uuid_1.2-1             
[105] reshape2_1.4.4          rjson_0.2.23            nlme_3.1-166            curl_6.1.0             
[109] repr_1.1.7              cachem_1.1.0            GlobalOptions_0.1.2     stringr_1.5.1          
[113] KernSmooth_2.23-26      parallel_4.4.0          pillar_1.10.1           grid_4.4.0             
[117] vctrs_0.6.5             gplots_3.2.0            dbplyr_2.5.0            cluster_2.1.8          
[121] evaluate_1.0.3          cli_3.6.3               compiler_4.4.0          rlang_1.1.4            
[125] crayon_1.5.3            plyr_1.8.9              fs_1.6.5                stringi_1.8.4          
[129] viridisLite_0.4.2       BiocParallel_1.38.0     munsell_0.5.1           Biostrings_2.72.1      
[133] lazyeval_0.2.2          pacman_0.5.1            GOSemSim_2.30.2         Matrix_1.7-1           
[137] IRkernel_1.3.2          hms_1.1.3               patchwork_1.3.0         bit64_4.6.0-1          
[141] future_1.34.0           ggplot2_3.5.1           KEGGREST_1.44.1         igraph_2.1.2           
[145] memoise_2.0.1           ggtree_3.12.0           fastmatch_1.1-6         bit_4.5.0.1            
[149] gson_0.1.0              ape_5.8-1              
>



