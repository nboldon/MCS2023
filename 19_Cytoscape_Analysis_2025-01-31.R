

# First, make sure Cytoscape is running
# Go to View -> Show Command Panel to enable the Command Panel
# Then open R and run the script.



######################################


## Test code to ensure connection to Cytoscape


# Clear workspace and reload RCy3
rm(list = ls())
if ("RCy3" %in% (.packages())) {
  detach("package:RCy3", unload = TRUE)
}
library(RCy3)

# Basic connection test function
test_cytoscape_connection <- function() {
  message("Step 1: Testing basic ping...")
  if (cytoscapePing()) {
    message("✓ Basic ping successful")
    
    message("\nStep 2: Testing version info...")
    tryCatch({
      version_info <- cytoscapeVersionInfo()
      message("✓ Version info retrieved successfully")
      message(paste("Cytoscape version:", version_info$cytoscapeVersion))
      
      message("\nStep 3: Testing simple network creation...")
      createNetworkFromDataFrames(
        nodes = data.frame(id = c("A", "B"), 
                           name = c("Node A", "Node B")),
        edges = data.frame(source = "A", 
                           target = "B", 
                           interaction = "test")
      )
      message("✓ Test network created successfully")
      
      return(TRUE)
    }, error = function(e) {
      message(paste("Error:", e$message))
      return(FALSE)
    })
  } else {
    message("× Failed to connect to Cytoscape")
    return(FALSE)
  }
}

# Run the test
message("Starting Cytoscape connection test...\n")
message("Please ensure:")
message("1. Cytoscape is running")
message("2. Command Panel is visible (View -> Show Command Panel)")
message("3. No other processes are using port 1234\n")

test_result <- test_cytoscape_connection()

if (test_result) {
  message("\nAll tests passed! You can now run your main analysis script.")
} else {
  message("\nConnection test failed. Please try these troubleshooting steps:")
  message("1. Restart Cytoscape")
  message("2. Wait 30 seconds after Cytoscape fully loads")
  message("3. Make sure Command Panel is enabled")
  message("4. Restart R Studio")
  message("5. Run this test script again")
}


###################################################
###################################################
###################################################


## Process files for input into Cytoscape


process_gene_file <- function(input_file, output_dir) {
  # First check Cytoscape connection
  message("Checking Cytoscape connection...")
  check_cytoscape()
  
  # Check packages and setup
  message("Setting up Ensembl connection...")
  ensembl <- setup_mart()
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read and process data
  message("Reading input file...")
  data <- read.csv(input_file)
  
  # Ensure necessary columns exist
  if (!all(c("name", "Log2FC", "FDR") %in% colnames(data))) {
    stop("Required columns 'name', 'Log2FC', and 'FDR' not found in input file.")
  }
  
  # Apply filtering criteria
  filtered_data <- data %>%
    filter((Log2FC >= 0.5 | Log2FC <= -0.5) & FDR < 0.05)
  
  message(paste("Filtered dataset contains", nrow(filtered_data), "genes"))
  
  # Split into up/down regulated groups
  data_up <- filter(filtered_data, Log2FC >= 0.5)
  data_down <- filter(filtered_data, Log2FC <= -0.5)
  
  message(paste("Found", nrow(data_up), "upregulated and", 
                nrow(data_down), "downregulated genes"))
  
  # Process each group
  for (type in c("Upregulated", "Downregulated")) {
    message(paste("\nProcessing", type, "genes..."))
    gene_set <- if(type == "Upregulated") data_up$name else data_down$name
    if (length(gene_set) > 0) {
      gene_data <- convert_genes(gene_set, ensembl)
      if (nrow(gene_data) > 0) {
        create_enhanced_network(gene_data, type, output_dir)
      } else {
        message(paste("No matching genes found for", type, "group"))
      }
    }
  }
}


##################################



setwd("/Volumes/DataBox/MCS2023/Tx_Comp/Cytoscape")

# Load necessary libraries
library(RCy3)
library(biomaRt)
library(dplyr)
library(igraph)

# Function to check Cytoscape connection
check_cytoscape <- function() {
  tryCatch({
    if (!cytoscapePing()) {
      stop("Cytoscape is not running. Please start Cytoscape and try again.")
    }
    
    # Test the connection more thoroughly
    cytoscapeVersionInfo()
    message("Successfully connected to Cytoscape")
    return(TRUE)
  }, error = function(e) {
    stop(paste("Cytoscape connection error:", e$message, 
               "\nPlease ensure Cytoscape is running and try these steps:",
               "\n1. Start Cytoscape",
               "\n2. Wait for it to fully load",
               "\n3. Ensure the Command Panel is enabled (View -> Show Command Panel)",
               "\n4. Try running the script again"))
  })
}

# Modified create_enhanced_network function with better error handling
create_enhanced_network <- function(gene_data, network_type, output_dir) {
  tryCatch({
    # First create a basic network to ensure connectivity
    message("Creating network...")
    
    # Convert data frame to ensure proper format
    nodes_df <- gene_data %>%
      select(entrezgene_id, external_gene_name, description) %>%
      rename(id = entrezgene_id, name = external_gene_name) %>%
      as.data.frame()
    
    # Create edges dataframe (even if empty)
    edges_df <- data.frame(source = character(), target = character())
    
    # Create network with explicit error checking
    network_name <- paste(network_type, "Gene Network")
    createNetworkFromDataFrames(
      nodes = nodes_df,
      edges = edges_df,
      title = network_name,
      collection = "Gene Networks"
    )
    
    # Verify network creation
    current_network <- getNetworkList()
    if (length(current_network) == 0) {
      stop("Network creation failed")
    }
    
    message("Applying visual style...")
    # Apply visual style with error checking
    style_name <- paste(network_type, "Style")
    createVisualStyle(style_name, defaults = list(
      NODE_SHAPE = "ellipse",
      NODE_SIZE = 60,
      NODE_FILL_COLOR = if(network_type == "Upregulated") "purple" else "darkgreen",
      EDGE_TRANSPARENCY = 120
    ))
    
    setVisualStyle(style_name)
    
    message("Exporting network...")
    # Export with full path specification
    export_path <- file.path(output_dir, 
                             paste0(tolower(network_type), "_network.png"))
    exportImage(export_path, type = "PNG")
    
    message(paste("Network exported to:", export_path))
    
  }, error = function(e) {
    stop(paste("Error in network creation:", e$message))
  })
}

# Modified main function
process_gene_file <- function(input_file, output_dir) {
  # First check Cytoscape connection
  message("Checking Cytoscape connection...")
  check_cytoscape()
  
  # Check packages and setup
  message("Setting up Ensembl connection...")
  ensembl <- setup_mart()
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read and process data
  message("Reading input file...")
  data <- read.csv(input_file)
  if (!"name" %in% colnames(data) || !"Log2FC" %in% colnames(data)) {
    stop("Required columns 'name' and 'Log2FC' not found in input file.")
  }
  
  # Split into up/down regulated
  data_up <- filter(data, Log2FC > 0)
  data_down <- filter(data, Log2FC < 0)
  
  message(paste("Found", nrow(data_up), "upregulated and", 
                nrow(data_down), "downregulated genes"))
  
  # Process each group
  for (type in c("Upregulated", "Downregulated")) {
    message(paste("\nProcessing", type, "genes..."))
    gene_set <- if(type == "Upregulated") data_up$name else data_down$name
    if (length(gene_set) > 0) {
      gene_data <- convert_genes(gene_set, ensembl)
      if (nrow(gene_data) > 0) {
        create_enhanced_network(gene_data, type, output_dir)
      } else {
        message(paste("No matching genes found for", type, "group"))
      }
    }
  }
}

# Usage
message("Starting analysis...")
input_file <- "/Volumes/DataBox/MCS2023/Tx_Comp/CellType_PairwiseComps/gaba_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv"
output_dir <- "/Volumes/DataBox/MCS2023/Tx_Comp/Cytoscape/"
process_gene_file(input_file, output_dir)




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
