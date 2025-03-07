

## NOTE: The 07_projMCS5.R code did not use abs(Log2FC)
## This code accommodates for that and reruns getMarkerFeatures before applying cutoffs

#Load libraries
R
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

#Verify matrices in project
getAvailableMatrices(projMCS5)

#Identify cell count  by cluster
table(projMCS5$Clusters)

# Save the table as a .csv file locally
#write.csv(projMCS5$Clusters, file = "/project/eon/nboldon/MCS2023/fragAnalysis_TSS7/Cluster_CellCounts.csv", row.names = FALSE)
# Does not print the table correctly

############################################
############################################

# Identify marker peaks and account for differences in data quality amongst cell groups
markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS5,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# To retrieve particular slices of the SummarizedExperiment; returns a list of DataFrame objects, one for each cell group
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")
markerList

# Save markerList as a .csv file locally
write.csv(markerList, file = "/Volumes/DataBox/MCS2023/FDR_0.01-Log2FC_1.25/2024-01_Cluster_Analysis/07_Peak_markerList_2025-03-06.csv", row.names = FALSE)


###########################################
###########################################


## Plotting marker peaks

# Marker peak heatmap
heatmapPeaks <-plotMarkerHeatmap(
	seMarker = markersPeaks,
	cutOff = "FDR<=0.01 & abs(Log2FC)>=1.25",
	transpose = TRUE
)

# To plot the heatmap
draw(heatmapPeaks, 
	heatmap_legend_side="bot",
	annotation_legend_side="bot"
)

# To save an editable vectorized version of the plot
plotPDF(heatmapPeaks, 
	name = "Peak-Marker-Heatmap_2025-03-06",
	width = 8, height = 6, 
	ArchRProj = projMCS5,
	addDOC = FALSE
)

