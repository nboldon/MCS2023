# Subsets by cluster, group by treatment; loop for all clusters


######################################################
######################################################


Setup an interactive session
salloc --account=eon -t 0-08:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Updated conda env 12-2023
module load miniconda3/23.1.0
conda activate archr2023_12

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(Seurat)
library(Signac)
library(BiocManager)
library(BiocGenerics)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtree)
library(aplot)
library(patchwork)
library(ggplot2)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ClusterID")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

############################################
############################################

# Define clusters vector
clusters <- paste0("C", 1:25)

# Loop through each cluster
for (cluster in clusters) {
  
  # Create subset for the current cluster
  subsetName <- paste0(cluster, "ArchRSubset")
  clusterData <- projMCS5[projMCS5$Clusters %in% cluster, ]
  
  # Specify treatment group for each sample
  treatment <- clusterData$Sample
  treatment <- gsub("C302_|C306_|C309_|C318_|C323_|C328_|C332_|C337_|C339_|C346_|C351_|C353_|C360_", "t1", treatment)
  treatment <- gsub("C304_|C308_|C312_|C349_|C315_|C321_|C324_|C355_|C327_|C330_|C333_|C358_|C336_|C342_|C348_|C362_", "t2", treatment)
  treatment <- gsub("C305_|C307_|C313_|C350_|C316_|C320_|C322_|C352_|C325_|C334_|C359_|C340_|C341_|C345_|C364_", "t3", treatment)
  treatment <- gsub("C301_|C303_|C310_|C314_|C319_|C335_|C338_|C344_|C354_|C356_|C361_|C363_", "t4", treatment)
  
  # Assign treatment to the clusterData object
  clusterData$treatment <- treatment
  
  # Get marker features for each cluster
  markerGS <- getMarkerFeatures(
    ArchRProj = clusterData,
    useMatrix = "GeneScoreMatrix",
    groupBy = "treatment",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    maxCells = 45000
  )
  
  # Get list of marker genes
  markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")
  
  # Save markerList to a CSV file
  outputFile <- paste0(cluster, "_markerList_", Sys.Date(), ".csv")
  write.csv(markerList, file = outputFile, row.names = FALSE)
  
  # Store the updated cluster data into a new object
  assign(subsetName, clusterData)
  
  # Print the current cluster to track progress
  print(paste0("Processed ", cluster))
}

