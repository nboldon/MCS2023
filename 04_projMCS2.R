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
addArchRThreads(threads = 16)

#Save Project
saveArchRProject(
	ArchRProj = projMCS2,
	outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS2", load = FALSE)

#########################################
#########################################

#Create Project 2
#Filter doublets
projMCS2 <- filterDoublets(projMCS1)
projMCS2

#To adjust TSS enrichment thresholds
# idxPass <- which(projMCS2$TSSEnrichment >=10)
# cellsPass <- projMCS2$cellNames[idxPass]
# projMCS2[cellsPass, ]

#If you want to filter more cells from the ArchRProject
# projMCSTmp <- filterDoublets(projMCS1, filterRatio = 1.5)
#To remove this from the R session
# rm(projMCSTmp)

#####################################
#####################################

##Dimensionality reduction with ArchR

#Create a reducedDims object called "IterativeLSI"
# projMCS2 <- addIterativeLSI(
#	ArchRProj = projMCS2,
#	useMatrix = "TileMatrix",
#	name = "IterativeLSI",
#	iterations = 2,
#	clusterParams = list( #See Seurat::FindClusters
#		resolution = c(0.2),
#		sampleCells = 10000,
#		n.start = 10
#	),
#	varFeatures = 25000,
#	dimsToUse = 1:30,
#	seed = 123
# )

#If you see downstream that you have subtle batch effects, another option is to add more LSI iterations and to start from a lower initial clustering resolution as shown below. Additionally, the number of variable features can be lowered to increase focus on the more variable features
projMCS2 <- addIterativeLSI(
	ArchRProj = projMCS2,
	useMatrix = "TileMatrix",
	name = "IterativeLSI2",
	iterations = 4,
	clusterParams = list ( 
		resolution = c(0.1, 0.2, 0.4),
		sampleCells = 10000,
		n.start = 10
	),
	varFeatures = 15000,
	dimsToUse = 1:30,
	seed = 123
)

#When the iterative LSI approach isn't enough of a correction for strong batch effect differences, the batch effect correction tool called Harmony can be used. This creates a new reducedDims object called "Harmony" in the MCS2 object
projMCS2 <- addHarmony(
	ArchRProj = projMCS2,
	reducedDims = "IterativeLSI2",
	name = "Harmony",
	groupBy = "Sample"
)

########################################
########################################

#Clustering using Seurat's FindClusters() function

projMCS2 <- addClusters(
	input = projMCS2,
	reducedDims = "IterativeLSI2",
	method = "Seurat",
	name = "Clusters",
	resolution = 0.8,
	seed = 123
)

#To access these clusters, use the $ accessor which shows the cluster ID for each single cell
head(projMCS2$Clusters)

#To tabulate the number of cells present in each cluster
table(projMCS2$Clusters)

#To better understand which samples reside in which clusters, create a cluster confusion matrix across each sample
cM <- confusionMatrix(paste0(projMCS2$Clusters), paste0(projMCS2$Sample))
cM

#To save the matrix as a .csv file
write.csv(cM, "SeuratClusters.csv")

#To plot the confusion matrix as a heatmap
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)
p

#Clustering using Scran
projMCS2 <- addClusters(
	input = projMCS2,
	reducedDims = "IterativeLSI2",
	method = "scran",
	name = "ScranClusters",
	k = 15,
	seed = 123
)

######################################
######################################

##Single Cell Embeddings - see 04_projMCS2_scEmbeddings.R
