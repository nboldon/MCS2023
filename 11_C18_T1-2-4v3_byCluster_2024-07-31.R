#Setup an interactive session
salloc --account=eon -t 0-02:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

## Pairwise comparisons using T1, T2, T4 vs T3

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

treatment <- gsub("C302_", "t1-2-4", treatment)
treatment <- gsub("C306_", "t1-2-4", treatment)
treatment <- gsub("C309_", "t1-2-4", treatment)
treatment <- gsub("C318_", "t1-2-4", treatment)
treatment <- gsub("C323_", "t1-2-4", treatment)
treatment <- gsub("C328_", "t1-2-4", treatment)
treatment <- gsub("C332_", "t1-2-4", treatment)
treatment <- gsub("C337_", "t1-2-4", treatment)
treatment <- gsub("C339_", "t1-2-4", treatment)
treatment <- gsub("C346_", "t1-2-4", treatment)
treatment <- gsub("C351_", "t1-2-4", treatment)
treatment <- gsub("C353_", "t1-2-4", treatment)
treatment <- gsub("C360_", "t1-2-4", treatment)
treatment <- gsub("C304_", "t1-2-4", treatment)
treatment <- gsub("C308_", "t1-2-4", treatment)
treatment <- gsub("C312_", "t1-2-4", treatment)
treatment <- gsub("C349_", "t1-2-4", treatment)
treatment <- gsub("C315_", "t1-2-4", treatment)
treatment <- gsub("C321_", "t1-2-4", treatment)
treatment <- gsub("C324_", "t1-2-4", treatment)
treatment <- gsub("C355_", "t1-2-4", treatment)
treatment <- gsub("C327_", "t1-2-4", treatment)
treatment <- gsub("C330_", "t1-2-4", treatment)
treatment <- gsub("C333_", "t1-2-4", treatment)
treatment <- gsub("C358_", "t1-2-4", treatment)
treatment <- gsub("C336_", "t1-2-4", treatment)
treatment <- gsub("C342_", "t1-2-4", treatment)
treatment <- gsub("C348_", "t1-2-4", treatment)
treatment <- gsub("C362_", "t1-2-4", treatment)
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
treatment <- gsub("C301_", "t1-2-4", treatment)
treatment <- gsub("C303_", "t1-2-4", treatment)
treatment <- gsub("C310_", "t1-2-4", treatment)
treatment <- gsub("C314_", "t1-2-4", treatment)
treatment <- gsub("C319_", "t1-2-4", treatment)
treatment <- gsub("C335_", "t1-2-4", treatment)
treatment <- gsub("C338_", "t1-2-4", treatment)
treatment <- gsub("C344_", "t1-2-4", treatment)
treatment <- gsub("C354_", "t1-2-4", treatment)
treatment <- gsub("C356_", "t1-2-4", treatment)
treatment <- gsub("C361_", "t1-2-4", treatment)
treatment <- gsub("C363_", "t1-2-4", treatment)

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

############################################

## Get marker features
c18.1MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.1ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.1MarkerComps <- getMarkers(c18.1MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.1MarkerComps, file = "C18.1_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.2MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.2ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.2MarkerComps <- getMarkers(c18.2MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.2MarkerComps, file = "C18.2_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.3MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.3ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.3MarkerComps <- getMarkers(c18.3MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.3MarkerComps, file = "C18.3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.4MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.4ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.4MarkerComps <- getMarkers(c18.4MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.4MarkerComps, file = "C18.4_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.5MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.5ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.5MarkerComps <- getMarkers(c18.5MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.5MarkerComps, file = "C18.5_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.6MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.6ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.6MarkerComps <- getMarkers(c18.6MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.6MarkerComps, file = "C18.6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.7MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.7ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.7MarkerComps <- getMarkers(c18.7MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.7MarkerComps, file = "C18.7_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.8MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.8MarkerComps <- getMarkers(c18.8MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.8MarkerComps, file = "C18.8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.9MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.9ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.9MarkerComps <- getMarkers(c18.9MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.9MarkerComps, file = "C18.9_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.10MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.10ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.10MarkerComps <- getMarkers(c18.10MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.10MarkerComps, file = "C18.10_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.11MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.11ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.11MarkerComps <- getMarkers(c18.11MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.11MarkerComps, file = "C18.11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.12MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.12ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.12MarkerComps <- getMarkers(c18.12MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.12MarkerComps, file = "C18.12_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.13MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.13ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.13MarkerComps <- getMarkers(c18.13MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.13MarkerComps, file = "C18.13_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.14MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.14ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.14MarkerComps <- getMarkers(c18.14MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.14MarkerComps, file = "C18.14_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.15MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.15ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.15MarkerComps <- getMarkers(c18.15MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.15MarkerComps, file = "C18.15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.16MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.16ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.16MarkerComps <- getMarkers(c18.16MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.16MarkerComps, file = "C18.16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.17MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.17ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.17MarkerComps <- getMarkers(c18.17MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.17MarkerComps, file = "C18.17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################

## Get marker features
c18.18MarkerComps <- getMarkerFeatures(
  ArchRProj = C18.18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1-2-4",
  bgdGroups = "t3"
)

#Get the list of markers
c18.18MarkerComps <- getMarkers(c18.18MarkerComps, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c18.18MarkerComps, file = "C18.18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv", row.names = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################
############################################

## Combine Spreadsheets

# Set working directory
setwd("/Volumes/DataBox/C18_Subset/")


# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C18.1_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df2 = subset(read.csv("C18.2_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df3 = subset(read.csv("C18.3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df4 = subset(read.csv("C18.4_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df5 = subset(read.csv("C18.5_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df6 = subset(read.csv("C18.6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df7 = subset(read.csv("C18.7_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df8 = subset(read.csv("C18.8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df9 = subset(read.csv("C18.9_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df10 = subset(read.csv("C18.10_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df11 = subset(read.csv("C18.11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df12 = subset(read.csv("C18.12_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df13 = subset(read.csv("C18.13_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df14 = subset(read.csv("C18.14_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df15 = subset(read.csv("C18.15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df16 = subset(read.csv("C18.16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df17 = subset(read.csv("C18.17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))
df18 = subset(read.csv("C18.18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv"), select = -c(group, group_name, idx, MeanDiff))

names(df1) = c("seqnames", "start", "end", "strand", "name", "C18.1.Log2FC", "C18.1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "C18.2.Log2FC", "C18.2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "C18.3.Log2FC", "C18.3.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "C18.4.Log2FC", "C18.4.FDR")
names(df5) = c("seqnames", "start", "end", "strand", "name", "C18.5.Log2FC", "C18.5.FDR")
names(df6) = c("seqnames", "start", "end", "strand", "name", "C18.6.Log2FC", "C18.6.FDR")
names(df7) = c("seqnames", "start", "end", "strand", "name", "C18.7.Log2FC", "C18.7.FDR")
names(df8) = c("seqnames", "start", "end", "strand", "name", "C18.8.Log2FC", "C18.8.FDR")
names(df9) = c("seqnames", "start", "end", "strand", "name", "C18.9.Log2FC", "C18.9.FDR")
names(df10) = c("seqnames", "start", "end", "strand", "name", "C18.10.Log2FC", "C18.10.FDR")
names(df11) = c("seqnames", "start", "end", "strand", "name", "C18.11.Log2FC", "C18.11.FDR")
names(df12) = c("seqnames", "start", "end", "strand", "name", "C18.12.Log2FC", "C18.12.FDR")
names(df13) = c("seqnames", "start", "end", "strand", "name", "C18.13.Log2FC", "C18.13.FDR")
names(df14) = c("seqnames", "start", "end", "strand", "name", "C18.14.Log2FC", "C18.14.FDR")
names(df15) = c("seqnames", "start", "end", "strand", "name", "C18.15.Log2FC", "C18.15.FDR")
names(df16) = c("seqnames", "start", "end", "strand", "name", "C18.16.Log2FC", "C18.16.FDR")
names(df17) = c("seqnames", "start", "end", "strand", "name", "C18.17.Log2FC", "C18.17.FDR")
names(df18) = c("seqnames", "start", "end", "strand", "name", "C18.18.Log2FC", "C18.18.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), 
                     list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "T1-2-4vT3_CompMarkers_Combined_2024-07-31.csv", row.names = FALSE)


############################################
############################################
############################################

sessionInfo()

R version 4.3.2 (2023-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.9 (Ootpa)

Matrix products: default
BLAS/LAPACK: /pfs/tc1/home/nboldon/.conda/envs/archr2023_12/lib/libopenblasp-r0.3.25.so;  LAPACK version 3.11.0

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
[13] Rcpp_1.0.12                 Matrix_1.6-4               
[15] GenomicRanges_1.54.1        GenomeInfoDb_1.38.5        
[17] IRanges_2.36.0              S4Vectors_0.40.2           
[19] BiocGenerics_0.48.1         matrixStats_1.2.0          
[21] data.table_1.14.10          stringr_1.5.1              
[23] plyr_1.8.9                  magrittr_2.0.3             
[25] ggplot2_3.4.4               gtable_0.3.4               
[27] gtools_3.9.5                gridExtra_2.3              
[29] ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3                 jsonlite_1.8.8                    
  [3] farver_2.1.1                       fs_1.6.3                          
  [5] BiocIO_1.12.0                      zlibbioc_1.48.0                   
  [7] vctrs_0.6.5                        Rsamtools_2.18.0                  
  [9] memoise_2.0.1                      Cairo_1.6-2                       
 [11] RCurl_1.98-1.14                    ggtree_3.8.2                      
 [13] S4Arrays_1.2.0                     Rhdf5lib_1.22.1                   
 [15] SparseArray_1.2.3                  gridGraphics_0.5-1                
 [17] cachem_1.0.8                       GenomicAlignments_1.38.2          
 [19] igraph_1.5.1                       lifecycle_1.0.4                   
 [21] pkgconfig_2.0.3                    R6_2.5.1                          
 [23] fastmap_1.1.1                      gson_0.1.0                        
 [25] GenomeInfoDbData_1.2.11            digest_0.6.33                     
 [27] aplot_0.2.2                        colorspace_2.1-0                  
 [29] patchwork_1.1.3                    RSQLite_2.3.5                     
 [31] fansi_1.0.6                        httr_1.4.7                        
 [33] polyclip_1.10-6                    abind_1.4-5                       
 [35] compiler_4.3.2                     bit64_4.0.5                       
 [37] withr_3.0.0                        downloader_0.4                    
 [39] BiocParallel_1.36.0                viridis_0.6.4                     
 [41] DBI_1.2.1                          ggforce_0.4.1                     
 [43] MASS_7.3-60                        DelayedArray_0.28.0               
 [45] rjson_0.2.21                       HDO.db_0.99.1                     
 [47] tools_4.3.2                        ape_5.7-1                         
 [49] scatterpie_0.2.1                   glue_1.7.0                        
 [51] restfulr_0.0.15                    nlme_3.1-164                      
 [53] GOSemSim_2.26.1                    rhdf5filters_1.12.1               
 [55] shadowtext_0.1.2                   reshape2_1.4.4                    
 [57] fgsea_1.26.0                       generics_0.1.3                    
 [59] BSgenome_1.70.1                    tidyr_1.3.1                       
 [61] tidygraph_1.2.3                    utf8_1.2.4                        
 [63] XVector_0.42.0                     ggrepel_0.9.4                     
 [65] pillar_1.9.0                       yulab.utils_0.1.0                 
 [67] splines_4.3.2                      dplyr_1.1.4                       
 [69] tweenr_2.0.2                       treeio_1.24.3                     
 [71] lattice_0.22-5                     rtracklayer_1.62.0                
 [73] bit_4.0.5                          tidyselect_1.2.0                  
 [75] GO.db_3.18.0                       Biostrings_2.70.1                 
 [77] RhpcBLASctl_0.23-42                graphlayouts_1.0.2                
 [79] stringi_1.8.3                      yaml_2.3.8                        
 [81] lazyeval_0.2.2                     ggfun_0.1.3                       
 [83] codetools_0.2-19                   BSgenome.Mmusculus.UCSC.mm10_1.4.3
 [85] ggraph_2.1.0                       tibble_3.2.1                      
 [87] qvalue_2.32.0                      ggplotify_0.1.2                   
 [89] cli_3.6.2                          munsell_0.5.0                     
 [91] png_0.1-8                          XML_3.99-0.16.1                   
 [93] blob_1.2.4                         DOSE_3.26.2                       
 [95] bitops_1.0-7                       viridisLite_0.4.2                 
 [97] tidytree_0.4.5                     scales_1.3.0                      
 [99] purrr_1.0.2                        crayon_1.5.2                      
[101] rlang_1.1.3                        cowplot_1.1.2                     
[103] fastmatch_1.1-4                    KEGGREST_1.42.0     
