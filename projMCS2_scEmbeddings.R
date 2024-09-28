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

##Single Cell Embeddings - UMAP

#To run a Uniform Manifold Approximation and Projection (UMAP)
projMCS2 <- addUMAP(
	ArchRProj = projMCS2,
	reducedDims = "IterativeLSI2",
	name = "UMAP",
	nNeighbors = 30,
	minDist = 0.5,
	metric = "cosine",
	seed = 123
)

#To plot the UMAP results by sample
p1 <- plotEmbedding(ArchRProj = projMCS2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#To plot the UMAP results by cluster
p2 <- plotEmbedding(ArchRProj = projMCS2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#To visualize these 2 plots side by side
ggAlignPlots(p1, p2, type = "h")

#To save an editable vectorized version of the plots:
plotPDF(p1,p2, name = "UMAP-Sample-Clusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

#plotEmbedding() can also be used to visualize the results from clustering using scran

p1 <- plotEmbedding(ArchRProj = projMCS2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = projMCS2, colorBy = "cellColData", name = "ScranClusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

#To save an editable vectorized version of this plot
plotPDF(p1,p2, name = "UMAP-Sample-ScranClusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)


#################################################
#################################################

##Single cell embeddings - tSNE

#To run a t-Stocastic Neighbor Embedding (tSNE)
projMCS2 <- addTSNE(
	ArchRProj = projMCS2,
	reducedDims = "IterativeLSI2",
	name = "TSNE",
	perplexity = 30,
	seed = 123
)

#To plot the tSNE (the same parameters apply to colorBy and name regardless of the type of embedding used)
p1 <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "cellColData",
	name = "Sample",
	embedding = "TSNE"
)

p2 <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "cellColData",
	name = "Clusters",
	embedding = "TSNE"
)

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "TSNE-Sample-Clusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

#######################################
#######################################

##Comparing clustering results from Seurat and scran
#####Should this be Seurat???
#####Did not run, as scran clustering failed previously
p1 <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "cellColData",
	name = "Sample",
	embedding = "TSNE"
)

p2 <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "cellColData",
	name = "ScranClusters",
	embedding = "TSNE"
)

ggAlignPlots(p1,p2, type = "h")

plotPDF(p1,p2, name = "tSNE-Sample-ScranClusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

###########################################
###########################################

##Dimensionality reduction after Harmony for UMAP

#Used to assess the effects of Harmony by visualizing the embedding using UMAP or tSNE and comparing this to the embeddings visualized previously for iterative LSI.

#Repeat the UMAP embedding with the same parameters but for the "Harmony" reduced Dims object
projMCS2 <- addUMAP(
	ArchRProj = projMCS2,
	reducedDims = "Harmony",
	name = "UMAPHarmony",
	nNeighbors = 30,
	minDist = 0.5, 
	metric = "cosine",
	seed = 123
)

p3 <- plotEmbedding( 
	ArchRProj = projMCS2,
	colorBy = "cellColData", 
	name = "Sample",
	embedding = "UMAPHarmony"
)

p4 <- plotEmbedding(
	ArchRProj = projMCS2,
	colorBy = "cellColData",
	name = "Clusters",
	embedding = "UMAPHarmony"
)

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

##Dimensionality reduction after Harmony for tSNE

#Follow similar steps to those used for UMAP

projMCS2 <- addTSNE(
        ArchRProj = projMCS2,
        reducedDims = "Harmony",
        name = "TSNEHarmony",
        perplexity = 30,
        seed = 123
)

p3 <- plotEmbedding(
        ArchRProj = projMCS2,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "TSNEHarmony"
)

p4 <- plotEmbedding(
        ArchRProj = projMCS2,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "TSNEHarmony"
)

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "TSNE2Harmony-Sample-Clusters.pdf", ArchRProj = projMCS2, addDOC = FALSE, width = 5, height = 5)

