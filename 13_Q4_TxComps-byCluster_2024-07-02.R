Setup an interactive session
salloc --account=eon -t 0-04:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

treatment <- C25ArchRSubset$Sample

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
C25ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C25ArchRSubset$treatment)

############

C25ArchRSubset <- addHarmony(
  ArchRProj = C25ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "treatment",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C25ArchRSubset <- addImputeWeights(C25ArchRSubset)
getImputeWeights(C25ArchRSubset)

############################################
############################################
############################################

## Get marker features
c25MarkersT1T2 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
c25markerListT1T2 <- getMarkers(c25MarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT1T2, file = "C25_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT1T3 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
c25markerListT1T3 <- getMarkers(c25MarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT1T3, file = "C25_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT1T4 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
c25markerListT1T4 <- getMarkers(c25MarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT1T4, file = "C25_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c25MarkersT2T1 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
c25markerListT2T1 <- getMarkers(c25MarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT2T1, file = "C25_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT2T3 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
c25markerListT2T3 <- getMarkers(c25MarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT2T3, file = "C25_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT2T4 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
c25markerListT2T4 <- getMarkers(c25MarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT2T4, file = "C25_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c25MarkersT3T1 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
c25markerListT3T1 <- getMarkers(c25MarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT3T1, file = "C25_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT3T2 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
c25markerListT3T2 <- getMarkers(c25MarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT3T2, file = "C25_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT3T4 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
c25markerListT3T4 <- getMarkers(c25MarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT3T4, file = "C25_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################
############################################

## Get marker features
c25MarkersT4T1 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
c25markerListT4T1 <- getMarkers(c25MarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT4T1, file = "C25_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT4T3 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
c25markerListT4T3 <- getMarkers(c25MarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT4T3, file = "C25_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################

## Get marker features
c25MarkersT4T2 <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
c25markerListT4T2 <- getMarkers(c25MarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(c25markerListT4T2, file = "C25_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv", row.names = FALSE)

############################################
############################################
############################################

## Combine tx group cluster-specific comparisons

setwd("/Volumes/DataBox/MCS2023/ProjMCS6/ResearchQuestions/")

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df2 = subset(read.csv("C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df3 = subset(read.csv("C18_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df4 = subset(read.csv("T3_C18_GeneMarkers_2024-04-10.csv"), select = -c(X, idx, MeanDiff))
names(df1) = c("seqnames", "start", "end", "strand", "name", "T3vsT1.Log2FC", "T3vsT1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "T3vsT2.Log2FC", "T3vsT2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "T3vsT4.Log2FC", "T3vsT4.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "T3.Log2FC", "T3.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C18_All-T3vT_TxMarkers_Combined_2024-07-06.csv", row.names = FALSE)

############################################
############################################

## Cluster 1

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C1_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df2 = subset(read.csv("C1_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df3 = subset(read.csv("C1_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df4 = subset(read.csv("C1_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"), select = -c(group, group_name, idx, MeanDiff))
df5 = subset(read.csv("C1_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df6 = subset(read.csv("C1_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df7 = subset(read.csv("C1_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df8 = subset(read.csv("C1_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df9 = subset(read.csv("C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df10 = subset(read.csv("C1_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df11 = subset(read.csv("C1_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df12 = subset(read.csv("C1_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df13 = subset(read.csv("C1_projMCS6_2024-06-20.csv"), select = -c(X, idx, MeanDiff))

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
names(df13) = c("seqnames", "start", "end", "strand", "name", "C1_projMCS6.Log2FC", "C1_projMCS6.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4, df5, df6, df7, df8, df9, 
                                                                  df10, df11, df12, df13))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C1_MCS6-TxComp_geneMarkersCombined_2024-07-06.csv", row.names = FALSE)

############################################
############################################

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("T1_C25_GeneMarkers_2024-04-10.csv"), select = -c(X, idx, MeanDiff))
df2 = subset(read.csv("T2_C25_GeneMarkers_2024-04-10.csv"), select = -c(X, idx, MeanDiff))
df3 = subset(read.csv("T3_C25_GeneMarkers_2024-04-10.csv"), select = -c(X, idx, MeanDiff))
df4 = subset(read.csv("T4_C25_GeneMarkers_2024-04-10.csv"), select = -c(X, idx, MeanDiff))
df5 = subset(read.csv("C1_projMCS6_2024-06-20.csv"), select = -c(X, idx, MeanDiff))

names(df1) = c("seqnames", "start", "end", "strand", "name", "T1.Log2FC", "T1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "T2.Log2FC", "T2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "T3.Log2FC", "T3.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "T4.Log2FC", "T4.FDR")
names(df5) = c("seqnames", "start", "end", "strand", "name", "MCS6.Log2FC", "MCS6.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4, df5))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C1_MCS6-TxMarkers_Combined_2024-07-05.csv", row.names = FALSE)

############################################
############################################

## To sort rows in results file:

setwd("/Volumes/DataBox/MCS2023/Tx_Comp/")

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df2 = subset(read.csv("C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df3 = subset(read.csv("C18_projMCS6_2024-06-20.csv"), select = -c(X, idx, MeanDiff))

names(df1) = c("seqnames", "start", "end", "strand", "name", "T3vsT1.Log2FC", "T3vsT1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "T3vsT2.Log2FC", "T3vsT2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "C18_projMCS6.Log2FC", "C18_projMCS6.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3))

# Export combined data frame to a new spreadsheet
write.csv(combined_df[order(-rowSums(is.na(combined_df))), ], "C18-ordered_PairwiseCombined_2024-07-08.csv", row.names = FALSE)

############################################
############################################

## Q6.2 Combining all files for Cluster 3 pairwise comparisons:
# Spreadsheets are in descending order from the least number of "NA" column for DA pairwise comparisons in each row of gene markers

setwd("/Volumes/DataBox/MCS2023/ProjMCS6/ResearchQuestions/")

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C3_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, idx, MeanDiff))
df2 = subset(read.csv("C3_GenotypeComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv"), select = -c(group, idx, MeanDiff))
df3 = subset(read.csv("C3_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df4 = subset(read.csv("C3_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df5 = subset(read.csv("C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df6 = subset(read.csv("C3_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv"), select = -c(group, group_name, idx, MeanDiff))
df7 = subset(read.csv("C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df8 = subset(read.csv("C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df9 = subset(read.csv("C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df10 = subset(read.csv("C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df11 = subset(read.csv("C3_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv"), select = -c(group, group_name, idx, MeanDiff))
df12 = subset(read.csv("C3_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df13 = subset(read.csv("C3_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df14 = subset(read.csv("C3_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv"), select = -c(group, group_name, idx, MeanDiff))
df15 = subset(read.csv("C3_projMCS6_2024-06-20.csv"), select = -c(X, idx, MeanDiff))

names(df1) = c("Q3.group_name", "seqnames", "start", "end", "strand", "name", "Q3.Log2FC", "Q3.FDR")
names(df2) = c("Q4b.group_name", "seqnames", "start", "end", "strand", "name", "Q4b.Log2FC", "Q4b.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "T1vsT2.Log2FC", "T1vsT2.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "T1vsT3.Log2FC", "T1vsT3.FDR")
names(df5) = c("seqnames", "start", "end", "strand", "name", "T1vsT4.Log2FC", "T1vsT4.FDR")
names(df6) = c("seqnames", "start", "end", "strand", "name", "T2vsT1.Log2FC", "T2vsT1.FDR")
names(df7) = c("seqnames", "start", "end", "strand", "name", "T2vsT3.Log2FC", "T2vsT3.FDR")
names(df8) = c("seqnames", "start", "end", "strand", "name", "T2vsT4.Log2FC", "T2vsT4.FDR")
names(df9) = c("seqnames", "start", "end", "strand", "name", "T3vsT1.Log2FC", "T3vsT1.FDR")
names(df10) = c("seqnames", "start", "end", "strand", "name", "T3vsT2.Log2FC", "T3vsT2.FDR")
names(df11) = c("seqnames", "start", "end", "strand", "name", "T3vsT4.Log2FC", "T3vsT4.FDR")
names(df12) = c("seqnames", "start", "end", "strand", "name", "T4vsT1.Log2FC", "T4vsT1.FDR")
names(df13) = c("seqnames", "start", "end", "strand", "name", "T4vsT2.Log2FC", "T4vsT2.FDR")
names(df14) = c("seqnames", "start", "end", "strand", "name", "T4vsT3.Log2FC", "T4vsT3.FDR")
names(df15) = c("seqnames", "start", "end", "strand", "name", "C3_projMCS6.Log2FC", "C3_projMCS6.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4, df5, df6, df7,
              df8, df9, df10, df11, df12, df13, df14, df15))

# Export combined data frame to a new spreadsheet
write.csv(combined_df[order(-rowSums(is.na(combined_df))), ], "C3_Ordered_PairwiseCombined_2024-07-08.csv", row.names = FALSE)



