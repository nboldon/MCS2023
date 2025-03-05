#Setup an interactive session
salloc --account=eon -t 0-02:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#conda activate seurat_2024-07

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

setwd("/project/eon/nboldon/MCS2023/Subset/C18_Subset")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

##############

#Load project
C18ArchRSubset <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Subset/Save-C18ArchRSubset",
                                   force = FALSE, showLogo = FALSE)

getAvailableMatrices(C18ArchRSubset)

table(C18ArchRSubset$Clusters)


############################################
############################################
############################################
############################################
############################################
############################################

## Research Question 3:

C18ArchRSubset <- addHarmony(
  ArchRProj = C18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Subset C18_Subset by Cluster
C18.1ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C1"]
C18.2ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C2"]
C18.3ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C3"]
C18.4ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C4"]
C18.5ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C5"]
C18.6ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C6"]
C18.7ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C7"]
C18.8ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C8"]
C18.9ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C9"]
C18.10ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C10"]
C18.11ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C11"]
C18.12ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C12"]
C18.13ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C13"]
C18.14ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C14"]
C18.15ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C15"]
C18.16ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C16"]
C18.17ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C17"]
C18.18ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C18"]

#####################

# Specify which treatment group each sample is in:
# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C18.18ArchRSubset$Sample

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
C18.18ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C18.18ArchRSubset$treatment)

############################

C18.18ArchRSubset <- addHarmony(
  ArchRProj = C18.18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C18.18ArchRSubset <- addImputeWeights(C18.18ArchRSubset)
getImputeWeights(C18.18ArchRSubset)

#############################

## Get marker features
c18.18Markers <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

#Get the list of markers
c18.18markerList <- getMarkers(c18.18Markers, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file:
write.csv(c18.18markerList, file = "C18.18_TxComp_FDR-0-1_Log2FC-0-5_2024-07-25.csv", row.names = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################

## Research Question 5:

C18ArchRSubset <- addHarmony(
  ArchRProj = C18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

###############

# Subset by Cluster
C18.1ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C1"]
C18.2ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C2"]
C18.3ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C3"]
C18.4ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C4"]
C18.5ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C5"]
C18.6ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C6"]
C18.7ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C7"]
C18.8ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C8"]
C18.9ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C9"]
C18.10ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C10"]
C18.11ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C11"]
C18.12ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C12"]
C18.13ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C13"]
C18.14ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C14"]
C18.15ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C15"]
C18.16ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C16"]
C18.17ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C17"]
C18.18ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C18"]

##################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C18.15ArchRSubset$Sample

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
C18.18ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C18.18ArchRSubset$treatment)

############

C18.18ArchRSubset <- addHarmony(
  ArchRProj = C18.18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C18.18ArchRSubset <- addImputeWeights(C18.18ArchRSubset)
getImputeWeights(C18.18ArchRSubset)

############################################

## Get marker features
c18.18MarkersT1T2 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
c18.18markerListT1T2 <- getMarkers(c18.18MarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT1T2, file = "C18.18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT1T3 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
c18.18markerListT1T3 <- getMarkers(c18.18MarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT1T3, file = "C18.18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT1T4 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
c18.18markerListT1T4 <- getMarkers(c18.18MarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT1T4, file = "C18.18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c18.18MarkersT2T1 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
c18.18markerListT2T1 <- getMarkers(c18.18MarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT2T1, file = "C18.18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT2T3 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
c18.18markerListT2T3 <- getMarkers(c18.18MarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT2T3, file = "C18.18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT2T4 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
c18.18markerListT2T4 <- getMarkers(c18.18MarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT2T4, file = "C18.18_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c18.18MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c18.18markerListT3T1 <- getMarkers(c18.18MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT3T1, file = "C18.18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c18.18markerListT3T2 <- getMarkers(c18.18MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT3T2, file = "C18.18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c18.18markerListT3T4 <- getMarkers(c18.18MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT3T4, file = "C18.18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c18.18MarkersT4T1 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
c18.18markerListT4T1 <- getMarkers(c18.18MarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT4T1, file = "C18.18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT4T3 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.18markerListT4T3 <- getMarkers(c18.18MarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT4T3, file = "C18.18_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkersT4T2 <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
c18.18markerListT4T2 <- getMarkers(c18.18MarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18markerListT4T2, file = "C18.18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv", row.names = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################

# Set working directory
setwd("/Volumes/DataBox/C18_Subset/")


## Cluster 18.1

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C18.18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df2 = subset(read.csv("C18.18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df3 = subset(read.csv("C18.18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df4 = subset(read.csv("C18.18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df5 = subset(read.csv("C18.18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df6 = subset(read.csv("C18.18_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df7 = subset(read.csv("C18.18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df8 = subset(read.csv("C18.18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df9 = subset(read.csv("C18.18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df10 = subset(read.csv("C18.18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df11 = subset(read.csv("C18.18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df12 = subset(read.csv("C18.18_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv"), select = -c(group, group_name, idx, MeanDiff))
df13 = subset(read.csv("C18.18_TxComp_FDR-0-1_Log2FC-0-5_2024-07-25.csv"), select = -c(group, idx, MeanDiff))

names(df1) = c("seqnames", "start", "end", "strand", "name", "T1vsT2.Log2FC", "T1vsT2.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "T1vsT3.Log2FC", "T1vsT3.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "T1vsT4.Log2FC", "T1vsT4.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "T2vsT1.Log2FC", "T2vsT1.FDR")
names(df5) = c("seqnames", "start", "end", "strand", "name", "T2vsT3.Log2FC", "T2vsT3.FDR")
names(df6) = c("seqnames", "start", "end", "strand", "name", "T2vsT4.Log2FC", "T2vsT4.FDR")
names(df7) = c("seqnames", "start", "end", "strand", "name", "T3vsT1.Log2FC", "T3vsT1.FDR")
names(df8) = c("seqnames", "start", "end", "strand", "name", "T3vsT2.Log2FC", "T3vsT2.FDR")
names(df9) = c("seqnames", "start", "end", "strand", "name", "T3vsT4.Log2FC", "T3vsT4.FDR")
names(df10) = c("seqnames", "start", "end", "strand", "name", "T4vsT1.Log2FC", "T4vsT1.FDR")
names(df11) = c("seqnames", "start", "end", "strand", "name", "T4vsT2.Log2FC", "T4vsT2.FDR")
names(df12) = c("seqnames", "start", "end", "strand", "name", "T4vsT3.Log2FC", "T4vsT3.FDR")
names(df13) = c("seqnames", "start", "end", "strand", "name", "group_name", "C18.18.Log2FC", "C18.18.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4, df5, df6, df7, df8, df9, 
                                                                  df10, df11, df12, df13))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C18.18-TxComp_geneMarkersCombined_2024-07-25.csv", row.names = FALSE)


############################################
############################################
############################################
############################################
############################################
############################################

7/25/2024 - ArchR_2023_12

7/24/2024 - Seurat_2024-07

> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.9 (Ootpa)

Matrix products: default
BLAS/LAPACK: /pfs/tc1/home/nboldon/miniconda3/envs/seurat_2024-07/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0

Random number generation:
  RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Denver
tzcode source: system (glibc)

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] presto_1.0.0                nabor_0.5.0                
 [3] harmony_1.2.0               pheatmap_1.0.12            
 [5] enrichplot_1.20.3           clusterProfiler_4.8.3      
 [7] org.Mm.eg.db_3.17.0         AnnotationDbi_1.64.1       
 [9] rhdf5_2.44.0                SummarizedExperiment_1.32.0
[11] Biobase_2.62.0              MatrixGenerics_1.14.0      
[13] Rcpp_1.0.13                 Matrix_1.6-5               
[15] GenomicRanges_1.54.1        GenomeInfoDb_1.38.5        
[17] IRanges_2.36.0              S4Vectors_0.40.2           
[19] BiocGenerics_0.48.1         matrixStats_1.3.0          
[21] data.table_1.15.4           stringr_1.5.1              
[23] plyr_1.8.9                  magrittr_2.0.3             
[25] ggplot2_3.5.1               gtable_0.3.5               
[27] gtools_3.9.5                gridExtra_2.3              
[29] ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3                 jsonlite_1.8.8                    
  [3] farver_2.1.2                       fs_1.6.4                          
  [5] BiocIO_1.12.0                      zlibbioc_1.48.0                   
  [7] vctrs_0.6.5                        Rsamtools_2.18.0                  
  [9] memoise_2.0.1                      RCurl_1.98-1.16                   
 [11] ggtree_3.8.2                       S4Arrays_1.2.0                    
 [13] Rhdf5lib_1.22.1                    SparseArray_1.2.3                 
 [15] gridGraphics_0.5-1                 cachem_1.1.0                      
 [17] GenomicAlignments_1.38.2           igraph_2.0.3                      
 [19] lifecycle_1.0.4                    pkgconfig_2.0.3                   
 [21] R6_2.5.1                           fastmap_1.2.0                     
 [23] gson_0.1.0                         GenomeInfoDbData_1.2.11           
 [25] digest_0.6.36                      aplot_0.2.3                       
 [27] colorspace_2.1-0                   patchwork_1.2.0                   
 [29] RSQLite_2.3.7                      fansi_1.0.6                       
 [31] httr_1.4.7                         polyclip_1.10-6                   
 [33] abind_1.4-5                        compiler_4.3.3                    
 [35] bit64_4.0.5                        withr_3.0.0                       
 [37] downloader_0.4                     BiocParallel_1.36.0               
 [39] viridis_0.6.5                      DBI_1.2.3                         
 [41] ggforce_0.4.2                      MASS_7.3-60                       
 [43] DelayedArray_0.28.0                rjson_0.2.21                      
 [45] HDO.db_0.99.1                      tools_4.3.3                       
 [47] ape_5.8                            scatterpie_0.2.3                  
 [49] glue_1.7.0                         restfulr_0.0.15                   
 [51] nlme_3.1-165                       GOSemSim_2.26.1                   
 [53] rhdf5filters_1.12.1                shadowtext_0.1.4                  
 [55] reshape2_1.4.4                     fgsea_1.26.0                      
 [57] generics_0.1.3                     BSgenome_1.70.1                   
 [59] tidyr_1.3.1                        tidygraph_1.3.1                   
 [61] utf8_1.2.4                         XVector_0.42.0                    
 [63] ggrepel_0.9.5                      pillar_1.9.0                      
 [65] yulab.utils_0.1.4                  splines_4.3.3                     
 [67] dplyr_1.1.4                        tweenr_2.0.3                      
 [69] treeio_1.24.3                      lattice_0.22-6                    
 [71] rtracklayer_1.62.0                 bit_4.0.5                         
 [73] tidyselect_1.2.1                   GO.db_3.18.0                      
 [75] Biostrings_2.70.1                  RhpcBLASctl_0.23-42               
 [77] graphlayouts_1.1.1                 stringi_1.8.4                     
 [79] yaml_2.3.9                         lazyeval_0.2.2                    
 [81] ggfun_0.1.5                        codetools_0.2-20                  
 [83] BSgenome.Mmusculus.UCSC.mm10_1.4.3 ggraph_2.2.1                      
 [85] tibble_3.2.1                       qvalue_2.32.0                     
 [87] ggplotify_0.1.2                    cli_3.6.3                         
 [89] munsell_0.5.1                      png_0.1-8                         
 [91] XML_3.99-0.17                      blob_1.2.4                        
 [93] DOSE_3.26.2                        bitops_1.0-7                      
 [95] viridisLite_0.4.2                  tidytree_0.4.6                    
 [97] scales_1.3.0                       purrr_1.0.2                       
 [99] crayon_1.5.3                       rlang_1.1.4                       
[101] cowplot_1.1.3                      fastmatch_1.1-4                   
[103] KEGGREST_1.42.0       


