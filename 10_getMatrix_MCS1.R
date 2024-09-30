#Setup an interactive session
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
library(tidyr)
library(cowplot)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ProjMCS7/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS1 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS1", force = FALSE, showLogo = FALSE)

projMCS7
getAvailableMatrices(projMCS1)

########################################
########################################

#tileMatrix <- getMatrixFromProject(
#  ArchRProj = projMCS1,
#  useMatrix = "TileMatrix",
# verbose = FALSE,
#  binarize = FALSE,
#  threads = getArchRThreads(),
#  logFile = createLogFile("getMatrixFromProject")
#)
# "Error in .getMatFromArrow(ArrowFile = ArrowFile, featureDF = featureDF,  : \n  Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!\n"

tileMatrix <- getMatrixFromProject(
  ArchRProj = projMCS1,
  useMatrix = "TileMatrix",
  verbose = FALSE,
  binarize = TRUE,
  threads = getArchRThreads(),
  binarize = TRUE,
  logFile = createLogFile("getMatrixFromProject")
)
# Error in getMatrixFromProject(ArchRProj = projMCS1, useMatrix = "TileMatrix",  : 
#formal argument "binarize" matched by multiple actual arguments

geneMatrix <- getMatrixFromProject(
  ArchRProj = projMCS1,
  useMatrix = "GeneScoreMatrix",
  verbose = FALSE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

tileMatrix
geneMatrix

saveRDS(tileMatrix, file = "tileMatrix_projMCS1.rds")
saveRDS(geneMatrix, file = "geneMatrix_projMCS1.rds")

########################################
########################################

setwd("/Volumes/DataBox/ProjMCS7/")

tileMatrix <- readRDS("tileMatrix_projMCS1.rds")

########################################
########################################

tileMatrix@colData$nFrags

# Extract the data from colData
nFrags <- tileMatrix@colData@listData[["nFrags"]]
clusters <- tileMatrix@colData@listData[["Clusters"]]
sample <- tileMatrix@colData@listData[["Sample"]]
peakReads <- tileMatrix@colData@listData[["ReadsInPeaks"]]
promotorReads <- tileMatrix@colData@listData[["ReadsInPromotor"]]
tssReads <- tileMatrix@colData@listData[["ReadsInTSS"]]

geneMatrix@colData$nFrags

# Extract the data from colData
nFrags_gene <- geneMatrix@colData@listData[["nFrags"]]
clusters_gene <- geneMatrix@colData@listData[["Clusters"]]
sample_gene <- geneMatrix@colData@listData[["Sample"]]
peakReads_gene <- geneMatrix@colData@listData[["ReadsInPeaks"]]
promotorReads_gene <- geneMatrix@colData@listData[["ReadsInPromotor"]]
tssReads_gene <- geneMatrix@colData@listData[["ReadsInTSS"]]

##########################################
##########################################

# Combine the data into a data frame
combined_tile_df <- data.frame(
  nFrags = nFrags,
  Clusters = clusters,
  Sample = sample,
  #ReadsInPeaks = peakReads,
  ReadsInPromotor = promotorReads,
  ReadsInTSS = tssReads
)

# Write the data frame to a CSV file
write.csv(combined_tile_df, file = "combined_TileData_projMCS1.csv", row.names = FALSE)

# Combine the data into a data frame
combined_gene_df <- data.frame(
  nFrags = nFrags_gene,
  #Clusters = clusters_gene,
  Sample = sample_gene,
  #ReadsInPeaks = peakReads_gene,
  #ReadsInPromotor = promotorReads_gene,
  ReadsInTSS = tssReads_gene
)

# When "Clusters = clusters_gene" is included, Error in data.frame(nFrags = nFrags_gene, Clusters = clusters_gene, Sample = sample_gene,  : 
# arguments imply differing number of rows: 106456, 0

# Write the data frame to a CSV file
write.csv(combined_gene_df, file = "combined_GeneData_projMCS1.csv", row.names = FALSE)

