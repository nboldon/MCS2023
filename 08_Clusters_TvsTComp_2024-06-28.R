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

############################################
############################################

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
############################################
############################################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C1ArchRSubset$Sample

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
C1ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C1ArchRSubset$treatment)

############

C1ArchRSubset <- addHarmony(
  ArchRProj = C1ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C1ArchRSubset <- addImputeWeights(C1ArchRSubset)
getImputeWeights(C1ArchRSubset)

############################################

## Get marker features
c1MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c1markerListT3T1 <- getMarkers(c1MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T1, file = "C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c1markerListT3T2 <- getMarkers(c1MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T2, file = "C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c1markerListT3T4 <- getMarkers(cMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T4, file = "C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 2 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C2ArchRSubset$Sample

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
C2ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C2ArchRSubset$treatment)

############

C2ArchRSubset <- addHarmony(
  ArchRProj = C2ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C2ArchRSubset <- addImputeWeights(C2ArchRSubset)
getImputeWeights(C2ArchRSubset)

############################################

############################################

## Get marker features
c2MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C2ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c2markerListT3T1 <- getMarkers(c2MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c2markerListT3T1, file = "C2_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c2MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C2ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c2markerListT3T2 <- getMarkers(c2MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c2markerListT3T2, file = "C2_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c2MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C2ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c2markerListT3T4 <- getMarkers(c2MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c2markerListT3T4, file = "C2_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 3 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C3ArchRSubset$Sample

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
C3ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C3ArchRSubset$treatment)

############

C3ArchRSubset <- addHarmony(
  ArchRProj = C3ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C3ArchRSubset <- addImputeWeights(C3ArchRSubset)
getImputeWeights(C3ArchRSubset)

############################################

############################################

## Get marker features
c3MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C3ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c3markerListT3T1 <- getMarkers(c3MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c3markerListT3T1, file = "C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c3MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C3ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c3markerListT3T2 <- getMarkers(c3MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c3markerListT3T2, file = "C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c3MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C3ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c3markerListT3T4 <- getMarkers(c3MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c3markerListT3T4, file = "C3_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Clusters 4, 5, 6, 7

############################################
############################################
############################################

## Cluster 8 

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

############

C8ArchRSubset <- addHarmony(
  ArchRProj = C8ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C8ArchRSubset <- addImputeWeights(C8ArchRSubset)
getImputeWeights(C8ArchRSubset)

############################################

############################################

## Get marker features
c8MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c8markerListT3T1 <- getMarkers(c8MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c8markerListT3T1, file = "C8_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c8MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c8markerListT3T2 <- getMarkers(c8MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c8markerListT3T2, file = "C8_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c8MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c8markerListT3T4 <- getMarkers(c8MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c8markerListT3T4, file = "C8_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Clusters 9, 10

############################################
############################################
############################################

## Cluster 11 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C11ArchRSubset$Sample

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
C11ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C11ArchRSubset$treatment)

############

C11ArchRSubset <- addHarmony(
  ArchRProj = C11ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C11ArchRSubset <- addImputeWeights(C11ArchRSubset)
getImputeWeights(C11ArchRSubset)

############################################

############################################

## Get marker features
c11MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C11ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c11markerListT3T1 <- getMarkers(c11MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c11markerListT3T1, file = "C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c11MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C11ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c11markerListT3T2 <- getMarkers(c11MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c11markerListT3T2, file = "C11_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c11MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C11ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c11markerListT3T4 <- getMarkers(c11MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c11markerListT3T4, file = "C11_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 12, 13

############################################
############################################
############################################

## Cluster 14 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C14ArchRSubset$Sample

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
C14ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C14ArchRSubset$treatment)

############

C14ArchRSubset <- addHarmony(
  ArchRProj = C14ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C14ArchRSubset <- addImputeWeights(C14ArchRSubset)
getImputeWeights(C14ArchRSubset)

############################################

############################################

## Get marker features
c14MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C14ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c14markerListT3T1 <- getMarkers(c14MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c14markerListT3T1, file = "C14_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c14MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C14ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c14markerListT3T2 <- getMarkers(c14MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c14markerListT3T2, file = "C14_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c14MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C14ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c14markerListT3T4 <- getMarkers(c14MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c14markerListT3T4, file = "C14_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Clusters 15, 16, 17

############################################
############################################
############################################

## Cluster 18 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C18ArchRSubset$Sample

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
C18ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C18ArchRSubset$treatment)

############

C18ArchRSubset <- addHarmony(
  ArchRProj = C18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C18ArchRSubset <- addImputeWeights(C18ArchRSubset)
getImputeWeights(C18ArchRSubset)

#Error: Error in H5Lexists(h5loc, name) : HDF5. Links. Can't get value.
#Error in .safelapply(seq_len(nRep), function(y) { : 
#    Error Found Iteration 1 : 
#    [1] "Error in !exists : invalid argument type\n"
#  <simpleError in !exists: invalid argument type>
#    In addition: Warning message:
#    In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
#    1 function calls resulted in an error
  
#Error: Error in .safelapply(seq_len(nRep), function(y) { : 
#  Error Found Iteration 1 : 
#    [1] "Error in H5Dopen(h5loc, name) : HDF5. Dataset. Can't open object.\n"
#  <simpleError in H5Dopen(h5loc, name): HDF5. Dataset. Can't open object.>
#Error Found Iteration 2 : 
#	[1] "Error in H5Dopen(h5loc, name) : HDF5. Dataset. Can't open object.\n"
#	<simpleError in H5Dopen(h5loc, name): HDF5. Dataset. Can't open object.>
#In addition: Warning message:
#In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
#  2 function calls resulted in an error

############################################

## Get marker features
c18MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c18markerListT3T1 <- getMarkers(c18MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18markerListT3T1, file = "C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c18MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c18markerListT3T2 <- getMarkers(c18MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18markerListT3T2, file = "C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c18MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c18markerListT3T4 <- getMarkers(c18MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18markerListT3T4, file = "C18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 19 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C1ArchRSubset$Sample

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
C1ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C1ArchRSubset$treatment)

############

C1ArchRSubset <- addHarmony(
  ArchRProj = C1ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C1ArchRSubset <- addImputeWeights(C1ArchRSubset)
getImputeWeights(C1ArchRSubset)

############################################

## Get marker features
c1MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c1markerListT3T1 <- getMarkers(c1MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T1, file = "C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c1markerListT3T2 <- getMarkers(c1MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T2, file = "C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c1markerListT3T4 <- getMarkers(cMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T4, file = "C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C1ArchRSubset$Sample

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
C1ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C1ArchRSubset$treatment)

############

C1ArchRSubset <- addHarmony(
  ArchRProj = C1ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C1ArchRSubset <- addImputeWeights(C1ArchRSubset)
getImputeWeights(C1ArchRSubset)

############################################

## Get marker features
c1MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c1markerListT3T1 <- getMarkers(c1MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T1, file = "C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c1markerListT3T2 <- getMarkers(c1MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T2, file = "C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c1markerListT3T4 <- getMarkers(cMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T4, file = "C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C1ArchRSubset$Sample

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
C1ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C1ArchRSubset$treatment)

############

C1ArchRSubset <- addHarmony(
  ArchRProj = C1ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C1ArchRSubset <- addImputeWeights(C1ArchRSubset)
getImputeWeights(C1ArchRSubset)

############################################

## Get marker features
c1MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c1markerListT3T1 <- getMarkers(c1MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T1, file = "C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c1markerListT3T2 <- getMarkers(c1MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T2, file = "C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c1markerListT3T4 <- getMarkers(cMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T4, file = "C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C1ArchRSubset$Sample

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
C1ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C1ArchRSubset$treatment)

############

C1ArchRSubset <- addHarmony(
  ArchRProj = C1ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C1ArchRSubset <- addImputeWeights(C1ArchRSubset)
getImputeWeights(C1ArchRSubset)

############################################

## Get marker features
c1MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c1markerListT3T1 <- getMarkers(c1MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T1, file = "C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c1markerListT3T2 <- getMarkers(c1MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T2, file = "C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################

## Get marker features
c1MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c1markerListT3T4 <- getMarkers(cMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c1markerListT3T4, file = "C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)

############################################
############################################
############################################

