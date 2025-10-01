## 04_projMCS2_MarkerGenes.R

- Group by: Cluster (no subset), FDR <= 0.01, abs(Log2FC) >= 1.25
- getMarkerFeatures
- Marker list by cluster (FDR & Log2FC)
- Heatmaps for marker features; cowplots for all genes
    *NOTE: projMCS2_MarkerGenes.R did not calculate abs(Log2FC)
## 04_Cluster_GeneMarkers_2025-03-06.R
- Uses 04_projMCS2_MarkerGenes.R code, but calculates abs(Log2FC)



#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=64G --nodes=1 --ntasks-per-node=16

#Load required dependencies
module load miniconda3/4.12.0
conda activate /pfs/tc1/project/eon/archr_env

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/Analysis")
addArchRGenome("mm10")
addArchRThreads (threads = 16)

#Load project
projMCS2 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS2", force = FALSE, showLogo = FALSE)

######################################
######################################

##Gene scores and marker genes

markerGS <- get MarkerFeatures(
	ArchRProj = projMCS2,
	useMatrix = "GeneScoreMatrix",
	groupBy = "Clusters",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

#Marker list by cluster
markerList$C6

Cluster.markers <- getMarkerFeatures(
	ArchRProj = projMCS2,
	useMatrix = "GeneScoreMatrix",
	groupBy = "Clusters",
	testMethod = "wilcoxon",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	maxCells = 4000
)

markerList.Cluster.markers <- data.frame(getMarkers(Cluster.markers, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.2"))

cluster.markers.test <- markerList.Cluster.markers

#To get a list of dataframe objects, one for each cluster, containing the relevant marker features
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

#Marker list by cluster
markerList$C6

for(i in names(markerList)) {
	write.csv(markerList[[i]], file = paste(i, ".csv", sep = ""))
}

#########################################
#########################################

##Heatmaps for marker features

#To visualize all of the marker features simultaneously
markerGenes <- c(
	"Olig2", "Neurod6", "Gad2", "Dio2", "C1qa", "Cd31"
)

heatmapGS <- plotMarkerHeatmap(
	seMarker = markersGS
	cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
	labelMarkers = markerGenes,
	transpose = TRUE
)

#To plot the heatmap
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#To save the heatmap
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projMCS2, addDOC = FALSE)

#############################################
#############################################

##Visualizing marker genes on an embedding
markerGenes <- c(
	"Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "GeneScoreMatrix",
	name = markerGenes,
	embedding = "UMAP",
	quantCut = c(0.01, 0.95),
	imputeWeights = NULL
)

#To plot a specific gene
p$Olig2

#To plot all genes, use cowplot to arrange the various marker genes into a single plot
p2 <- lapply(p, function(x){
	x + guides(color = FALSE, fill = FALSE) +
	theme_ArchR(baseSize = 6.5) +
	theme(plot.margin = unit(c(0,0,0,0), "cm")) +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank()
	)
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#To save an editable vectorized version of the plot
plotPDF(plotList = p, 
	name = "UMAP-MarkerGenes-WO-Imputation.pdf",
	ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

#######################################################
#######################################################

##Use MAGIC to impute gene scores by smoothing signal across nearby cells

#Impute weights to the ArchRProject
projMCS2 <- addImputeWeights(projMCS2)

#The impute weights can then be passed to plotEmbedding() when plotting gene scores overlayed on the UMAP embedding
markerGenes <- c(
	"Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "GeneScoreMatrix",
	name = markerGenes,
	embedding = "UMAP",
	imputeWeights = getImputeWeights(projMCS2)
)

#To subset the plot list to select a specific gene
p$Olig2

#To plot all the marker genes at once using cowplot
#Rearrange for grid plotting
p2<- lapply(p, function(x){
	x + guides(color = FALSE, fill = FALSE) +
	theme_ArchR(baseSize = 6.5) +
	theme(plot.margin = unit(c(0,0,0,0), "cm")) +
	theme(
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank()
	)
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#To save an editable vectorized version of the plot
plotPDF(plotList = p,
	name = "UMAP-MarkerGenes-W-Imputation.pdf",
	ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)
