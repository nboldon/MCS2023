## 07_projMCS5.R
- Adds peak matrix
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- Plots marker peaks
  - cutOff = "FDR<=0.01 & Log2FC>=1.25"
- Adds motif annotations
  - motifSet = "cisbp"
NOTE: The 07_projMCS5.R code did not use abs(Log2FC)
The 07_projMCS5-w-ABS-Log2FC.R code accommodates for that and reruns getMarkerFeatures before applying cutoffs





#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

############################################
############################################

# To prepare for downstream analyses, we can create a new ArchRProject called projTinCann5 and add a new matrix to it containing insertion counts within our new merged peak set
projMCS5 <- addPeakMatrix(projMCS4)

# Identify marker peaks and account for differences in data quality amongst cell groups
markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS5,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# To retrieve particular slices of the SummarizedExperiment; returns a list of DataFrame objects, one for each cell group
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList

# Save markerList as a .csv file locally
write.csv(markerList, file = "/project/eon/nboldon/MCS2023/fragAnalysis_TSS7/Peak_markerList.csv", row.names = FALSE)

## Not all of the clusters have marker peaks identified, may want to change to a different FDR or Log2FC to get more peaks
##     although I'm not convinced there's a need to really identify peaks for each cluster

saveArchRProject(ArchRProj = projMCS5, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", load = FALSE)

###########################################
###########################################

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

#Verify matrices in project
getAvailableMatrices(projMCS5)

#Identify cell count  by cluster
table(projMCS5$Clusters)

# Save the table as a .csv file locally
#write.csv(projMCS5$Clusters, file = "/project/eon/nboldon/MCS2023/fragAnalysis_TSS7/Cluster_CellCounts.csv", row.names = FALSE)
# Does not print the table correctly

###########################################
###########################################

## Plotting marker peaks

# Marker peak heatmap
heatmapPeaks <-plotMarkerHeatmap(
	seMarker = markersPeaks,
	cutOff = "FDR<=0.01 & Log2FC>=1.25",
	transpose = TRUE
)

# To plot the heatmap
draw(heatmapPeaks, 
	heatmap_legend_side="bot",
	annotation_legend_side="bot"
)

# To save an editable vectorized version of the plot
plotPDF(heatmapPeaks, 
	name = "Peak-Marker-Heatmap",
	width = 8, height = 6, 
	ArchRProj = projMCS5,
	addDOC = FALSE
)

###########################################
###########################################

# Marker peak MA & Volcano plots can be run once cell groups are established. 

###########################################
###########################################

## Motif & Feature Enrichment

# Motif enrichment in differential peaks

# Add motif annotations to the ArchRProject to look for motifs that are enriched in peaks that are up or down in various cell types.
# This creates a binary matrix where the presence of a motif in each peak is indicated numerically. 
projMCS5 <-addMotifAnnotations(
	ArchRProj = projMCS5,
	motifSet = "cisbp",
	name = "Motif",
)

###########################################
###########################################

saveArchRProject(ArchRProj = projMCS5, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", load = FALSE)
