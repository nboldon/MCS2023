
# Set working directory (change this to your directory if not already set)
setwd("/Volumes/DataBox/Heatmap_Comps")

# Print the working directory to confirm it's set correctly
cat("Current working directory:", getwd(), "\n")

# Load necessary libraries
library(dplyr)
library(pheatmap)


############################################

mat<-read.csv("/Volumes/DataBox/MCS2023/Tx_Comp/C3_Combined_GeneMarkers_Comps_2024-11-07.csv",header=TRUE,row.names=1,
              check.names=FALSE)

head(mat)

############################################
############################################
############################################

## Makes a heatmap
pdf(file="C3_Combined_GeneMarkers_Comps_Heatmap_2024-12-19.pdf", width=15, height=80)
pheatmap(mat,scale="row",
         color=colorRampPalette(c("#440154FF", "white", "#1F9E89FF"))(50))
dev.off()

###############

## Loop to make heatmaps for all clusters

# Replace all NA values with zero


# Loop through clusters 1 to 25
for (cluster in 1:25) {
  
  # Define the input CSV file for the current cluster within the working directory
  input_file <- sprintf("/Volumes/DataBox/MCS2023/Tx_Comp/C%d_Combined_GeneMarkers_Comps_2024-11-07.csv", cluster)
  
  # Check if the file exists before proceeding
  if (file.exists(input_file)) {
    
    # Read the data
    data <- read.csv(input_file, row.names = "name")
    
    # Select only the Log2FC columns for the heatmap
    mat <- data %>% select(contains("Log2FC")) %>% as.matrix()
    
    # Replace NA values with 0 (or another preferred value)
    mat[is.na(mat)] <- 0
    
    # Check if any rows are entirely zero, and filter them out
    mat <- mat[rowSums(mat != 0) > 0, ]
    
    # Skip if matrix is empty after filtering
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      cat("Skipping cluster", cluster, "- no data to plot.\n")
      next
    }
    
    # Define the output PDF file name for the current cluster within the working directory
    pdf_file <- sprintf("C%d_Combined_GeneMarkers_Comps_Heatmap_2024-12-19.pdf", cluster)
    
    # Create and save the heatmap
    pdf(file = pdf_file, width = 15, height = 180)
    pheatmap(mat, scale = "row",
             color = colorRampPalette(c("#440154FF", "white", "#1F9E89FF"))(50))
    dev.off()
    
    cat("Processed and saved heatmap for cluster:", cluster, "\n")
  } else {
    cat("Skipping cluster", cluster, "- file does not exist in:", getwd(), "\n")
  }
}

## No data to plot for clusters 5, 7, 9

###############

## Loop to make heatmaps for all clusters

# Remove rows with NA values

# Loop through clusters 1 to 25
for (cluster in 1:25) {
  
  # Define the input CSV file for the current cluster
  input_file <- sprintf("/Volumes/DataBox/MCS2023/Tx_Comp/C%d_Combined_GeneMarkers_Comps_2024-11-07.csv", cluster)
  
  # Check if the file exists before proceeding
  if (file.exists(input_file)) {
    
    # Read the data
    data <- read.csv(input_file, row.names = "name")
    
    # Select only the Log2FC columns for the heatmap (assuming these are for the treatment groups)
    mat <- data %>% select(contains("Log2FC")) %>% as.matrix()
    
    # Filter genes (rows) that have no missing values in all treatment groups (columns)
    mat <- mat[rowSums(is.na(mat)) == 0, ]  # Only keep rows with no NA values
    
    # Check if matrix is empty after filtering
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      cat("Skipping cluster", cluster, "- no genes with data across all treatment groups.\n")
      next
    }
    
    # Define the output PDF file name for the current cluster
    pdf_file <- sprintf("C%d_Combined_GeneMarkers_Comps_Heatmap_NA-removed_2024-12-19.pdf", cluster)
    
    # Create and save the heatmap
    pdf(file = pdf_file, width = 15, height = 80)
    pheatmap(mat, scale = "row",
             color = colorRampPalette(c("#440154FF", "white", "#1F9E89FF"))(50))
    dev.off()
    
    cat("Processed and saved heatmap for cluster:", cluster, "\n")
  } else {
    cat("Skipping cluster", cluster, "- file does not exist.\n")
  }
}

# No genes with data across all treatment groups for clusters 2, 4, 5, 7, 9



############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################


## Heatmaps for cell type groupings by treatment group

# Removes rows with NA values


# Define the cell types and their corresponding files
cell_types <- c(
  "Microglia",
  "Astrocyte",
  "AstrocytePrecursor",
  "Oligodendrocyte",
  "OligoPrecursor",
  "GABAergic",
  "Glutaminergic",
  "GlutPrecursor",
  "GlutOligo",
  "GlutAstro",
  "EndoVasc"
)

# Loop through each cell type
for (cell_type in cell_types) {
  
  # Define the input CSV file for the current cell type
  input_file <- sprintf("/Volumes/DataBox/MCS2023/Tx_Comp/%s_Combined_GeneMarkers_Comps_2024-12-02.csv", cell_type)
  
  # Check if the file exists before proceeding
  if (file.exists(input_file)) {
    
    # Read the data
    data <- read.csv(input_file, row.names = "name")
    
    # Select only the Log2FC columns for the heatmap (assuming these are for the treatment groups)
    mat <- data %>% select(contains("Log2FC")) %>% as.matrix()
    
    # Filter genes (rows) that have no missing values in all treatment groups (columns)
    mat <- mat[rowSums(is.na(mat)) == 0, ]  # Only keep rows with no NA values
    
    # Check if matrix is empty after filtering
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      cat("Skipping cell type", cell_type, "- no genes with data across all treatment groups.\n")
      next
    }
    
    # Define the output PDF file name for the current cell type
    pdf_file <- sprintf("%s_Combined_GeneMarkers_Comps_Heatmap_NA-removed_2024-12-19.pdf", cell_type)
    
    # Create and save the heatmap
    pdf(file = pdf_file, width = 15, height = 80)
    pheatmap(mat, scale = "row",
             color = colorRampPalette(c("#440154FF", "white", "#1F9E89FF"))(50))
    dev.off()
    
    cat("Processed and saved heatmap for cell type:", cell_type, "\n")
  } else {
    cat("Skipping cell type", cell_type, "- file does not exist.\n")
  }
}


# Skipping cell type GlutOligo - no genes with data across all treatment groups.


############################################

## Sort by treatment group of interest and print zscores for the heatmap












############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################

## Heatmaps sorted by treatment group
# Adapted from /Volumes/DataBox/Scripts/projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R



known_markers_C3 <-c("Olig1", "Mbp", "Opalin", "Spock3", "S100b", "Mag", "Mog", "Cldn11", "Ugt8a")

#Subset the markerGS object to just these known markers
se_idx <- which(rowData(markerGS)$name %in% known_markers_All127_C18)  

#Get an index of these from the summarized experiment
subset_markerGS <- markerGS[se_idx,]  

##Plot it out:

pdf(file="C18_All127_TxComp_Heatmap_2024-03-13.pdf", width=15, height=30)
plotMarkerHeatmap(   ### Copied this over from the full function reference - leaving anything not commented at defaults
  seMarker = subset_markerGS,   ## Set this to the correct object containing the markers
  cutOff = "FDR <= 1 & Log2FC >=0.01",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = known_markers_All127_C18,  #Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)
dev.off()

######################################

# To print the Z-scores for the above heatmap

z_scores <- plotMarkerHeatmap(
  seMarker = subset_markerGS,
  cutOff = "FDR <= 1 & Log2FC >=.01",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = known_markers_All127_C18,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  # Change this to TRUE
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

# Print z-scores from plotMarkerHeatmap
write.csv(z_scores, "C18_All127_byTx_zscores_2024-03-13.csv", row.names = TRUE)
