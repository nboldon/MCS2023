Setup an interactive session
salloc --account=eon -t 0-010:00:00 --mem=256G --nodes=4 --ntasks-per-node=16

# Updated conda env 12-2023
module load miniconda3/23.1.0
conda activate archr2023_12


#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
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
setwd("/project/eon/nboldon/MCS2023/ProjMCS7/peak_txComps/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS7)
table(projMCS7$Clusters)

############################################
############################################

## Loop for all clusters

# Define the list of clusters
clusters <- paste0("C", 1:25)

# Define the treatment mapping
treatment_mapping <- list(
  "t1" = c("C302_", "C306_", "C309_", "C318_", "C323_", "C328_", "C332_", "C337_", "C339_", "C346_", "C351_", "C353_", "C360_"),
  "t2" = c("C304_", "C308_", "C312_", "C349_", "C315_", "C321_", "C324_", "C355_", "C327_", "C330_", "C333_", "C358_", "C336_", "C342_", "C348_", "C362_"),
  "t3" = c("C305_", "C307_", "C313_", "C350_", "C316_", "C320_", "C322_", "C352_", "C325_", "C334_", "C359_", "C340_", "C341_", "C345_", "C364_"),
  "t4" = c("C301_", "C303_", "C310_", "C314_", "C319_", "C335_", "C338_", "C344_", "C354_", "C356_", "C361_", "C363_")
)

# Iterate through all clusters
for (cluster in clusters) {
  # Subset by Cluster
  archRSubset <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Map treatment
  treatment <- archRSubset$Sample
  for (treatment_key in names(treatment_mapping)) {
    for (pattern in treatment_mapping[[treatment_key]]) {
      treatment <- gsub(pattern, treatment_key, treatment)
    }
  }
  
  # Check to make sure it worked
  print(unique(treatment))
  
  # Assign the treatment to the actual project
  archRSubset$treatment <- treatment
  
  # Check that this worked
  print(head(archRSubset$treatment))
  
  # Use MAGIC to impute gene scores by smoothing signal across nearby cells
  archRSubset <- addImputeWeights(archRSubset, reducedDims = "IterativeLSI2")
  getImputeWeights(archRSubset)
  
  # Define treatment pairs for comparison
  treatment_pairs <- list(
    c("t1", "t2"),
    c("t1", "t3"),
    c("t1", "t4"),
    c("t2", "t1"),
    c("t2", "t3"),
    c("t2", "t4"),
    c("t3", "t1"),
    c("t3", "t2"),
    c("t3", "t4"),
    c("t4", "t1"),
    c("t4", "t2"),
    c("t4", "t3")
  )
  
  # Iterate through treatment pairs
  for (pair in treatment_pairs) {
    useGroup <- pair[1]
    bgdGroup <- pair[2]
    
    # Get marker features
    markerFeatures <- getMarkerFeatures(
      ArchRProj = archRSubset,
      useMatrix = "PeakMatrix",
      groupBy = "treatment",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = useGroup,
      bgdGroups = bgdGroup
    )
    
    # Get the list of markers
    markerList <- getMarkers(markerFeatures, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")
    
    # Define file name
    file_name <- paste0(cluster, "_", useGroup, "vs", bgdGroup, "_Peak_FDR-0-1_Log2FC-0-5_2024-09-09.csv")
    
    # Write markerList to a CSV file
    write.csv(markerList, file = file_name, row.names = FALSE)
  }
}


############################################
############################################
############################################
############################################
############################################
############################################
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



