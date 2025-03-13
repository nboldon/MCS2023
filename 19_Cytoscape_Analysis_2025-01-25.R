


## DOES NOT RUN



# Load necessary libraries
library(RCy3)
library(biomaRt)
library(dplyr)

# Set up the Ensembl mart for mouse (mm10)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

# Set the file paths
input_file <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"
output_dir <- "/Volumes/DataBox/KEGG_Cytoscape"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to convert gene symbols to Entrez IDs using biomaRt
convert_gene_symbols_to_entrez <- function(gene_symbols) {
  # Clean the gene symbols by removing non-alphanumeric characters
  gene_symbols <- gsub("[^[:alnum:]_]", "", gene_symbols)  # Remove any non-alphanumeric characters
  
  # Get the conversion from gene symbol to Entrez ID using biomaRt
  conversion <- tryCatch({
    getBM(attributes = c("external_gene_name", "entrezgene_id"),
          filters = "external_gene_name", 
          values = gene_symbols, 
          mart = ensembl)
  }, error = function(e) {
    message(paste("Error in biomaRt query:", e$message))
    return(NULL)
  })
  
  # Filter out NAs from the conversion
  if (is.null(conversion) || nrow(conversion) == 0) {
    return(NULL)  # Return NULL if no valid Entrez IDs found
  }
  
  # Remove rows with missing Entrez IDs
  conversion <- conversion[!is.na(conversion$entrezgene_id), ]
  
  # Return Entrez IDs mapped to gene symbols
  return(setNames(conversion$entrezgene_id, conversion$external_gene_name))
}

# Load the data for the specific file
tryCatch({
  message(paste("Processing file:", input_file))
  data <- read.csv(input_file)
  
  if (nrow(data) == 0) {
    warning(paste("The file", input_file, "is empty. Skipping."))
  } else {
    # Ensure column names match expected values
    if (!("name" %in% colnames(data))) {
      stop("The 'name' column is missing in the input file.")
    }
    if (!("Log2FC" %in% colnames(data))) {
      stop("The 'Log2FC' column is missing in the input file.")
    }
    
    # Ensure Log2FC is numeric
    data$Log2FC <- as.numeric(data$Log2FC)
    
    # Filter out rows with NA in Log2FC or name
    data <- data %>% filter(!is.na(Log2FC) & !is.na(name) & name != "")
    
    # Extract upregulated and downregulated genes based on Log2FC
    data_up <- filter(data, Log2FC > 0)  # Upregulated genes
    data_down <- filter(data, Log2FC < 0) # Downregulated genes
    
    # Extract gene symbols for upregulated and downregulated genes
    genes_up <- data_up$name
    genes_down <- data_down$name
    
    if (length(genes_up) == 0 && length(genes_down) == 0) {
      warning(paste("No upregulated or downregulated genes found in", input_file, ". Skipping."))
    } else {
      # Convert gene symbols to Entrez IDs using biomaRt
      entrez_up <- if (length(genes_up) > 0) {
        convert_gene_symbols_to_entrez(genes_up)
      } else { NULL }
      
      entrez_down <- if (length(genes_down) > 0) {
        convert_gene_symbols_to_entrez(genes_down)
      } else { NULL }
      
      # Check results before creating networks
      if (!is.null(entrez_up) && length(entrez_up) > 0) {
        print("Upregulated Entrez IDs:")
        print(entrez_up)
      } else {
        message("No valid Entrez IDs for upregulated genes")
      }
      
      if (!is.null(entrez_down) && length(entrez_down) > 0) {
        print("Downregulated Entrez IDs:")
        print(entrez_down)
      } else {
        message("No valid Entrez IDs for downregulated genes")
      }
      
      # Check if Cytoscape is running
      if (!cytoscapePing()) {
        stop("Cytoscape is not running. Please start Cytoscape.")
      }
      
      # Generate the Cytoscape network if there is data
      # Generate the Cytoscape network if there is data
      if (!is.null(entrez_up) && length(entrez_up) > 0) {
        up_network <- data.frame(
          from = as.character(entrez_up),  # Ensure Entrez IDs are characters
          to = rep("Upregulated", length(entrez_up))  # Ensure category is character
        )
        createNetworkFromDataFrames(up_network, title = "Upregulated Genes Network", collection = "Gene Networks")
        exportNetworkAsImage(file.path(output_dir, paste("upregulated_network.png")))
      }
      
      if (!is.null(entrez_down) && length(entrez_down) > 0) {
        down_network <- data.frame(
          from = as.character(entrez_down),  # Ensure Entrez IDs are characters
          to = rep("Downregulated", length(entrez_down))  # Ensure category is character
        )
        createNetworkFromDataFrames(down_network, title = "Downregulated Genes Network", collection = "Gene Networks")
        exportNetworkAsImage(file.path(output_dir, paste("downregulated_network.png")))
      }
    }
  }
}, error = function(e) {
  # Catch any errors during the processing of the file
  message(paste("Error processing file", input_file, ":", e$message))
})





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




## LOOP DOES NOT RUN



## Open Cytoscape desk app


cytoscapePing()
cytoscapeDisconnect()
cytoscapeConnect()



# Load necessary libraries
library(RCy3)
library(biomaRt)

# Set up the Ensembl mart for mouse (mm10)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

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
    
    # Check if Cytoscape is running
    if (!cytoscapePing()) {
      stop("Cytoscape is not running. Please start Cytoscape.")
    }
    
    # Generate the Cytoscape network if there is data
    if (!is.null(entrez_up) && length(entrez_up) > 0) {
      up_network <- data.frame(
        from = entrez_up,  # Entrez IDs of upregulated genes
        to = rep("Upregulated", length(entrez_up))  # Assigning upregulated gene category
      )
      createNetworkFromDataFrames(up_network, title = paste("Upregulated Genes Network -", file_base), collection = "Gene Networks")
      exportNetworkAsImage(file.path(output_dir, paste(file_base, "_upregulated_network.png", sep = "")))
    }
    
    if (!is.null(entrez_down) && length(entrez_down) > 0) {
      down_network <- data.frame(
        from = entrez_down,  # Entrez IDs of downregulated genes
        to = rep("Downregulated", length(entrez_down))  # Assigning downregulated gene category
      )
      createNetworkFromDataFrames(down_network, title = paste("Downregulated Genes Network -", file_base), collection = "Gene Networks")
      exportNetworkAsImage(file.path(output_dir, paste(file_base, "_downregulated_network.png", sep = "")))
    }
    
  }, error = function(e) {
    # Catch any errors during the processing of each file
    message(paste("Error processing file", file, ":", e$message))
  })
}
