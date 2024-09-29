Setup an interactive session
salloc --account=eon -t 0-05:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS6)
table(projMCS6$Clusters)

############

projMCS6 <- addHarmony(
        ArchRProj = projMCS6,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
        force = TRUE
)

############################################
############################################

# Subset by Cluster
C1ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C1"]
C2ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C2"]
C3ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C3"]
C4ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C4"]
C5ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C5"]
C6ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C6"]
C7ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C7"]
C8ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C8"]
C9ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C9"]
C10ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C10"]
C11ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C11"]
C12ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C12"]
C13ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C13"]
C14ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C14"]
C15ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C15"]
C16ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C16"]
C17ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C17"]
C18ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C18"]
C19ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C19"]
C20ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C20"]
C21ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C21"]
C22ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C22"]
C23ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C23"]
C24ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C24"]
C25ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C25"]

############################################

# Specify which genotype group each sample is in:

# g1 = 2N
# g1 = 2N+
# g2 = Ts
# g2 = Ts+

genotype <- C2ArchRSubset$Sample

genotype <- gsub("C302_", "g1", genotype)
genotype <- gsub("C306_", "g1", genotype)
genotype <- gsub("C309_", "g1", genotype)
genotype <- gsub("C318_", "g1", genotype)
genotype <- gsub("C323_", "g1", genotype)
genotype <- gsub("C328_", "g1", genotype)
genotype <- gsub("C332_", "g1", genotype)
genotype <- gsub("C337_", "g1", genotype)
genotype <- gsub("C339_", "g1", genotype)
genotype <- gsub("C346_", "g1", genotype)
genotype <- gsub("C351_", "g1", genotype)
genotype <- gsub("C353_", "g1", genotype)
genotype <- gsub("C360_", "g1", genotype)
genotype <- gsub("C304_", "g1", genotype)
genotype <- gsub("C308_", "g1", genotype)
genotype <- gsub("C312_", "g1", genotype)
genotype <- gsub("C349_", "g1", genotype)
genotype <- gsub("C315_", "g1", genotype)
genotype <- gsub("C321_", "g1", genotype)
genotype <- gsub("C324_", "g1", genotype)
genotype <- gsub("C355_", "g1", genotype)
genotype <- gsub("C327_", "g1", genotype)
genotype <- gsub("C330_", "g1", genotype)
genotype <- gsub("C333_", "g1", genotype)
genotype <- gsub("C358_", "g1", genotype)
genotype <- gsub("C336_", "g1", genotype)
genotype <- gsub("C342_", "g1", genotype)
genotype <- gsub("C348_", "g1", genotype)
genotype <- gsub("C362_", "g1", genotype)
genotype <- gsub("C305_", "g2", genotype)
genotype <- gsub("C307_", "g2", genotype)
genotype <- gsub("C313_", "g2", genotype)
genotype <- gsub("C350_", "g2", genotype)
genotype <- gsub("C316_", "g2", genotype)
genotype <- gsub("C320_", "g2", genotype)
genotype <- gsub("C322_", "g2", genotype)
genotype <- gsub("C352_", "g2", genotype)
genotype <- gsub("C325_", "g2", genotype)
genotype <- gsub("C334_", "g2", genotype)
genotype <- gsub("C359_", "g2", genotype)
genotype <- gsub("C340_", "g2", genotype)
genotype <- gsub("C341_", "g2", genotype)
genotype <- gsub("C345_", "g2", genotype)
genotype <- gsub("C364_", "g2", genotype)
genotype <- gsub("C301_", "g2", genotype)
genotype <- gsub("C303_", "g2", genotype)
genotype <- gsub("C310_", "g2", genotype)
genotype <- gsub("C314_", "g2", genotype)
genotype <- gsub("C319_", "g2", genotype)
genotype <- gsub("C335_", "g2", genotype)
genotype <- gsub("C338_", "g2", genotype)
genotype <- gsub("C344_", "g2", genotype)
genotype <- gsub("C354_", "g2", genotype)
genotype <- gsub("C356_", "g2", genotype)
genotype <- gsub("C361_", "g2", genotype)
genotype <- gsub("C363_", "g2", genotype)

# Check to make sure it worked
unique(genotype)

# Assign the treatment to the actual project
C2ArchRSubset$genotype <- genotype

# Check that this worked - if not, make sure the previous line was run successfully
head(C2ArchRSubset$genotype)

############

############

C2ArchRSubset <- addHarmony(
        ArchRProj = C2ArchRSubset,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "genotype",
        force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C2ArchRSubset <- addImputeWeights(C2ArchRSubset)
getImputeWeights(C2ArchRSubset)

############################################

## Get marker features
geneMarkers <- getMarkerFeatures(
  ArchRProj = C2ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "genotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  #maxCells = 45000
)

#Get the list of markers
geneMarkerList <- getMarkers(geneMarkers, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

#############################

# Write markerList to a CSV file: 
write.csv(geneMarkerList, file = "C2_GenotypeComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv", row.names = FALSE)
