
## NOTE: code for 04_projMCS2_MarkerGenes.R was run using Log2FC>=1.25 (not abs)
## This code was run later using projMCS5 and includes abs(Log2FC)
# Initial marker gene summaries were created using FDR values (not Log2FC)


## 04_projMCS5_MarkerGenes_2025-03-06.R 

- Group by: Cluster (no subset), FDR <= 0.01, abs(Log2FC) >= 1.25
- Files saved for each cluster; ex: C1_MarkerGenes_2025-03-06.csv





#Load libraries
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

#Additional setup
setwd("/Volumes/DataBox/MCS2023/FDR_0.01-Log2FC_1.25/2024-01_Cluster_Analysis/")
addArchRGenome("mm10")
addArchRThreads (threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

######################################
######################################

##Gene scores and marker genes

markerGS <- getMarkerFeatures(
  ArchRProj = projMCS5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 4500
)

markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

#Marker list by cluster
markerList$C6

for(i in names(markerList)) {
  write.csv(markerList[[i]], file = paste(i, "_MarkerGenes_2025-03-06.csv", sep = ""))
}

#########################################
#########################################

