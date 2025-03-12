#Setup an interactive session
salloc --account=eon -t 0-06:00 --mem=64G --nodes=1 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

############################################
############################################

projMCS6 <- addHarmony(
  ArchRProj = projMCS6,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
projMCS6 <- addImputeWeights(projMCS6)
getImputeWeights(projMCS6)

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
############################################
############################################

## Subset by cluster

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C8ArchRSubset$Sample

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

# Check to make sure it worked
unique(treatment)

# Assign the treatment to the actual project
C8ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C8ArchRSubset$treatment)

############################################
############################################
############################################

markerGenesC8 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC8

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC8, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C8-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC8 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC8

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC8, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C8-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################
############################################
############################################

markerGenesC21 <- getMarkerFeatures(
  ArchRProj = C21ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC21, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C21-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC21 <- getMarkerFeatures(
  ArchRProj = C21ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC21, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C21-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################
############################################
############################################

markerGenesC22 <- getMarkerFeatures(
  ArchRProj = C22ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC22, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C22-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC22 <- getMarkerFeatures(
  ArchRProj = C22ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC22

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC22, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C22-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)
