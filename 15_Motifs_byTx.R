## Compares motifs after subsetting by treatment group

#Load libraries
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
library(RColorBrewer)

#Additional setup
setwd("/Volumes/DataBox/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)


# Load the ArchR project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

# Check the project details and available matrices
projMCS7
getAvailableMatrices(projMCS7)


########################################
########################################
########################################


# Set up treatment groups by renaming sample IDs
treatment <- projMCS7$Sample

# Rename treatment groups based on sample IDs
treatment <- gsub("C302_", "t1", treatment)
treatment <- gsub("C306_", "t1", treatment)
treatment <- gsub("C309_", "t1", treatment)
treatment <- gsub("C318_", "t1", treatment)
treatment <- gsub("C323_", "t1", treatment)
treatment <- gsub("C328_", "t1", treatment)
treatment <- gsub("C332_", "t1", treatment)
treatment <- gsub("C337_", "t1", treatment)
treatment <- gsub("C339_", "t1", treatment)
treatment <- gsub("C346_", "t1", treatment)
treatment <- gsub("C351_", "t1", treatment)
treatment <- gsub("C353_", "t1", treatment)
treatment <- gsub("C360_", "t1", treatment)
treatment <- gsub("C304_", "t2", treatment)
treatment <- gsub("C308_", "t2", treatment)
treatment <- gsub("C312_", "t2", treatment)
treatment <- gsub("C349_", "t2", treatment)
treatment <- gsub("C315_", "t2", treatment)
treatment <- gsub("C321_", "t2", treatment)
treatment <- gsub("C324_", "t2", treatment)
treatment <- gsub("C355_", "t2", treatment)
treatment <- gsub("C327_", "t2", treatment)
treatment <- gsub("C330_", "t2", treatment)
treatment <- gsub("C333_", "t2", treatment)
treatment <- gsub("C358_", "t2", treatment)
treatment <- gsub("C336_", "t2", treatment)
treatment <- gsub("C342_", "t2", treatment)
treatment <- gsub("C348_", "t2", treatment)
treatment <- gsub("C362_", "t2", treatment)
treatment <- gsub("C305_", "t3", treatment)
treatment <- gsub("C307_", "t3", treatment)
treatment <- gsub("C313_", "t3", treatment)
treatment <- gsub("C350_", "t3", treatment)
treatment <- gsub("C316_", "t3", treatment)
treatment <- gsub("C320_", "t3", treatment)
treatment <- gsub("C322_", "t3", treatment)
treatment <- gsub("C352_", "t3", treatment)
treatment <- gsub("C325_", "t3", treatment)
treatment <- gsub("C334_", "t3", treatment)
treatment <- gsub("C359_", "t3", treatment)
treatment <- gsub("C340_", "t3", treatment)
treatment <- gsub("C341_", "t3", treatment)
treatment <- gsub("C345_", "t3", treatment)
treatment <- gsub("C364_", "t3", treatment)
treatment <- gsub("C301_", "t4", treatment)
treatment <- gsub("C303_", "t4", treatment)
treatment <- gsub("C310_", "t4", treatment)
treatment <- gsub("C314_", "t4", treatment)
treatment <- gsub("C319_", "t4", treatment)
treatment <- gsub("C335_", "t4", treatment)
treatment <- gsub("C338_", "t4", treatment)
treatment <- gsub("C344_", "t4", treatment)
treatment <- gsub("C354_", "t4", treatment)
treatment <- gsub("C356_", "t4", treatment)
treatment <- gsub("C361_", "t4", treatment)
treatment <- gsub("C363_", "t4", treatment)

# Assign the modified treatment names back to the project
projMCS7$treatment <- treatment


########################################
########################################

## Subset treatment groups of interest

t1Subset <- projMCS7[
  projMCS7$treatment %in% c("t1"),]
t1Subset

# numberOfCells: 24969
# medianTSS: 14.968
# medianFrags: 7696

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t1Subset),
  outputDirectory = "t1Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t1Subset)


############################################

t2Subset <- projMCS7[
  projMCS7$treatment %in% c("t2"),]
t2Subset

# numberOfCells: 30149
# medianTSS: 14.749
# medianFrags: 6860

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t2Subset),
  outputDirectory = "t2Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t2Subset)

############################################

t3Subset <- projMCS7[
  projMCS7$treatment %in% c("t3"),]
t3Subset

# numberOfCells: 27879
# medianTSS: 14.332
# medianFrags: 6533

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t3Subset),
  outputDirectory = "t3Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t3Subset)

############################################

t4Subset <- projMCS7[
  projMCS7$treatment %in% c("t4"),]
t4Subset

# numberOfCells: 21458
# medianTSS: 13.96
# medianFrags: 5522

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t4Subset),
  outputDirectory = "t4Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t4Subset)

############################################
############################################
############################################
############################################


## Repeat the following for each treatment group


t4MarkersPeaks <- getMarkerFeatures(
  ArchRProj = t4Subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

t4MarkersPeaks

# Motif enrichment in marker peaks

t4EnrichMotifs <- peakAnnoEnrichment(
  seMarker = t4MarkersPeaks,
  ArchRProj = t4Subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
)

t4EnrichMotifs

t4HeatmapEM <- plotEnrichHeatmap(t4EnrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(t4HeatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(t4HeatmapEM, name = "t4-Motifs-Enriched-Marker-Heatmap_2025-01-11", width = 8, height = 6, ArchRProj = t4Subset, addDOC = FALSE)

#############
#############
#############

## Unable to save subset projects

# Error in saveArchRProject(ArchRProj = t1Subset, outputDirectory = "/Volumes/DataBox/Save-t1Subset",  : 
# all(file.exists(zfiles)) is not TRUE

saveArchRProject(ArchRProj = t1Subset, outputDirectory = "/Volumes/DataBox/Save-t1Subset", load = FALSE)

t1Subset <- loadArchRProject(path = "/Volumes/DataBox/Save-t1Subset", force = FALSE, showLogo = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################


