

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
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtree)
library(aplot)
library(patchwork)
library(ggplot2)

#Additional setup
setwd("/Volumes/DataBox/MCS2023/Tx_Comp/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS7)
table(projMCS7$Clusters)


############################################
############################################


## Subset projMCS7 by cell type

glut_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C18", "C19", "C21"), ]
#glutPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C15", "C16", "C17", "C20", "C25"), ]
gaba_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C22", "C23"), ]
microglia_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C10", "C11"), ]
oligo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C2", "C3"), ]
#oligoPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C5", "C6"), ]
#glutOligo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C4"]
endo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C12", "C13", "C14"), ]
astro_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C8"]
#astroPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C1", "C7", "C9"), ]
#glutAstro_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C24"]


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


## Astro pairwise comparisons


# Specify which treatment group each sample is in:
# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+


treatment <- astro_ArchRSubset$Sample

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
astro_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(astro_ArchRSubset$treatment)


############################################
############################################
############################################


## Get marker features
astroMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
astroMarkerListT3T1 <- getMarkers(astroMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT3T1, file = "astro_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
astroMarkerListT3T2 <- getMarkers(astroMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT3T2, file = "astro_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
astroMarkerListT3T4 <- getMarkers(astroMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT3T4, file = "astro_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
astroMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
astroMarkerListT1T3 <- getMarkers(astroMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT1T3, file = "astro_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
astroMarkerListT1T2 <- getMarkers(astroMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT1T2, file = "astro_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
astroMarkerListT1T4 <- getMarkers(astroMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT1T4, file = "astro_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
astroMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
astroMarkerListT2T1 <- getMarkers(astroMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT2T1, file = "astro_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
astroMarkerListT2T3 <- getMarkers(astroMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT2T3, file = "astro_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
astroMarkerListT2T4 <- getMarkers(astroMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT2T4, file = "astro_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
astroMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
astroMarkerListT4T1 <- getMarkers(astroMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT4T1, file = "astro_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
astroMarkerListT4T2 <- getMarkers(astroMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT4T2, file = "astro_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
astroMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
astroMarkerListT4T3 <- getMarkers(astroMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(astroMarkerListT4T3, file = "astro_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################


## Glut pairwise comparisons

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- glut_ArchRSubset$Sample

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
glut_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(glut_ArchRSubset$treatment)


############################################
############################################
############################################
############################################
############################################


## Get marker features
glutMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
glutMarkerListT3T1 <- getMarkers(glutMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT3T1, file = "glut_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
glutMarkerListT3T2 <- getMarkers(glutMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT3T2, file = "glut_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
glutMarkerListT3T4 <- getMarkers(glutMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT3T4, file = "glut_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
glutMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
glutMarkerListT1T3 <- getMarkers(glutMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT1T3, file = "glut_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
glutMarkerListT1T2 <- getMarkers(glutMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT1T2, file = "glut_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
glutMarkerListT1T4 <- getMarkers(glutMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT1T4, file = "glut_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################
############################################
############################################


## Get marker features
glutMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
glutMarkerListT2T1 <- getMarkers(glutMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT2T1, file = "glut_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
glutMarkerListT2T3 <- getMarkers(glutMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT2T3, file = "glut_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
glutMarkerListT2T4 <- getMarkers(glutMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT2T4, file = "glut_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################
############################################
############################################


## Get marker features
glutMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
glutMarkerListT4T1 <- getMarkers(glutMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT4T1, file = "glut_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
glutMarkerListT4T2 <- getMarkers(glutMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT4T2, file = "glut_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
glutMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = glut_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
glutMarkerListT4T3 <- getMarkers(glutMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(glutMarkerListT4T3, file = "glut_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)




############################################
############################################
############################################
############################################
############################################
############################################


## Oligo pairwise comparisons

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- oligo_ArchRSubset$Sample

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
oligo_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(oligo_ArchRSubset$treatment)

############################################
############################################
############################################
############################################
############################################

## Get marker features
oligoMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
oligoMarkerListT3T1 <- getMarkers(oligoMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT3T1, file = "oligo_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
oligoMarkerListT3T2 <- getMarkers(oligoMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT3T2, file = "oligo_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
oligoMarkerListT3T4 <- getMarkers(oligoMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT3T4, file = "oligo_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################

## Get marker features
oligoMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
oligoMarkerListT1T3 <- getMarkers(oligoMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT1T3, file = "oligo_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
oligoMarkerListT1T2 <- getMarkers(oligoMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT1T2, file = "oligo_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
oligoMarkerListT1T4 <- getMarkers(oligoMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT1T4, file = "oligo_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################

## Get marker features
oligoMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
oligoMarkerListT2T1 <- getMarkers(oligoMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT2T1, file = "oligo_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
oligoMarkerListT2T3 <- getMarkers(oligoMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT2T3, file = "oligo_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
oligoMarkerListT2T4 <- getMarkers(oligoMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT2T4, file = "oligo_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################

## Get marker features
oligoMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
oligoMarkerListT4T1 <- getMarkers(oligoMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT4T1, file = "oligo_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
oligoMarkerListT4T2 <- getMarkers(oligoMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT4T2, file = "oligo_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
oligoMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = oligo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
oligoMarkerListT4T3 <- getMarkers(oligoMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(oligoMarkerListT4T3, file = "oligo_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)



############################################
############################################
############################################
############################################
############################################
############################################



## Microglia pairwise comparisons


# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- microglia_ArchRSubset$Sample

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
microglia_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(microglia_ArchRSubset$treatment)


############################################
############################################
############################################


## Get marker features
microgliaMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
microgliaMarkerListT3T1 <- getMarkers(microgliaMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT3T1, file = "microglia_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
microgliaMarkerListT3T2 <- getMarkers(microgliaMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT3T2, file = "microglia_T3vsT2Comp_FDR-0-1_Log2FC-0-5_202e-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
microgliaMarkerListT3T4 <- getMarkers(microgliaMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT3T4, file = "microglia_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
microgliaMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
microgliaMarkerListT1T3 <- getMarkers(microgliaMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT1T3, file = "microglia_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
microgliaMarkerListT1T2 <- getMarkers(microgliaMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT1T2, file = "microglia_T1vsT2Comp_FDR-0-1_Log2FC-0-5_202e-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
microgliaMarkerListT1T4 <- getMarkers(microgliaMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT1T4, file = "microglia_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################
############################################
############################################


## Get marker features
microgliaMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
microgliaMarkerListT2T1 <- getMarkers(microgliaMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT2T1, file = "microglia_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
microgliaMarkerListT2T3 <- getMarkers(microgliaMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT2T3, file = "microglia_T2vsT3Comp_FDR-0-1_Log2FC-0-5_202e-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
microgliaMarkerListT2T4 <- getMarkers(microgliaMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT2T4, file = "microglia_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################
############################################
############################################


## Get marker features
microgliaMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
microgliaMarkerListT4T1 <- getMarkers(microgliaMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT4T1, file = "microglia_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
microgliaMarkerListT4T2 <- getMarkers(microgliaMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT4T2, file = "microglia_T4vsT2Comp_FDR-0-1_Log2FC-0-5_202e-01-18.csv", row.names = FALSE)

############################################

## Get marker features
microgliaMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
microgliaMarkerListT4T3 <- getMarkers(microgliaMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(microgliaMarkerListT4T3, file = "microglia_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)



############################################
############################################
############################################
############################################
############################################
############################################


## GABA pairwise comparisons


# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- gaba_ArchRSubset$Sample

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
gaba_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(gaba_ArchRSubset$treatment)


############################################
############################################
############################################
############################################
############################################


## Get marker features
gabaMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
gabaMarkerListT3T1 <- getMarkers(gabaMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT3T1, file = "gaba_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
gabaMarkerListT3T2 <- getMarkers(gabaMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT3T2, file = "gaba_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
gabaMarkerListT3T4 <- getMarkers(gabaMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT3T4, file = "gaba_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
gabaMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
gabaMarkerListT1T3 <- getMarkers(gabaMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT1T3, file = "gaba_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
gabaMarkerListT1T2 <- getMarkers(gabaMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT1T2, file = "gaba_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
gabaMarkerListT1T4 <- getMarkers(gabaMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT1T4, file = "gaba_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
gabaMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
gabaMarkerListT2T3 <- getMarkers(gabaMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT2T3, file = "gaba_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
gabaMarkerListT2T1 <- getMarkers(gabaMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT2T1, file = "gaba_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
gabaMarkerListT2T4 <- getMarkers(gabaMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT2T4, file = "gaba_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
gabaMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
gabaMarkerListT4T1 <- getMarkers(gabaMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT4T1, file = "gaba_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
gabaMarkerListT4T2 <- getMarkers(gabaMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT4T2, file = "gaba_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
gabaMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = gaba_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
gabaMarkerListT4T3 <- getMarkers(gabaMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(gabaMarkerListT4T3, file = "gaba_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)



############################################
############################################
############################################
############################################
############################################
############################################


## Endo-Vasc pairwise comparisons


# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- endo_ArchRSubset$Sample

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
endo_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(endo_ArchRSubset$treatment)


############################################
############################################
############################################
############################################
############################################


## Get marker features
endoMarkersT3T1 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t1"
)

#Get the list of markers
endoMarkerListT3T1 <- getMarkers(endoMarkersT3T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT3T1, file = "endo_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT3T2 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t2"
)

#Get the list of markers
endoMarkerListT3T2 <- getMarkers(endoMarkersT3T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT3T2, file = "endo_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT3T4 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t3",
  bgdGroups = "t4"
)

#Get the list of markers
endoMarkerListT3T4 <- getMarkers(endoMarkersT3T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT3T4, file = "endo_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
endoMarkersT1T3 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t3"
)

#Get the list of markers
endoMarkerListT1T3 <- getMarkers(endoMarkersT1T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT1T3, file = "endo_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT1T2 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t2"
)

#Get the list of markers
endoMarkerListT1T2 <- getMarkers(endoMarkersT1T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT1T2, file = "endo_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT1T4 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t1",
  bgdGroups = "t4"
)

#Get the list of markers
endoMarkerListT1T4 <- getMarkers(endoMarkersT1T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT1T4, file = "endo_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
endoMarkersT2T1 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t1"
)

#Get the list of markers
endoMarkerListT2T1 <- getMarkers(endoMarkersT2T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT2T1, file = "endo_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT2T3 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t3"
)

#Get the list of markers
endoMarkerListT2T3 <- getMarkers(endoMarkersT2T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT2T3, file = "endo_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT2T4 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t2",
  bgdGroups = "t4"
)

#Get the list of markers
endoMarkerListT2T4 <- getMarkers(endoMarkersT2T4, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT2T4, file = "endo_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)


############################################
############################################
############################################


## Get marker features
endoMarkersT4T1 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t1"
)

#Get the list of markers
endoMarkerListT4T1 <- getMarkers(endoMarkersT4T1, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT4T1, file = "endo_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT4T2 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t2"
)

#Get the list of markers
endoMarkerListT4T2 <- getMarkers(endoMarkersT4T2, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT4T2, file = "endo_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)

############################################

## Get marker features
endoMarkersT4T3 <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "t4",
  bgdGroups = "t3"
)

#Get the list of markers
endoMarkerListT4T3 <- getMarkers(endoMarkersT4T3, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

# Write markerList to a CSV file: 
write.csv(endoMarkerListT4T3, file = "endo_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv", row.names = FALSE)



############################################
############################################
############################################
############################################
############################################
############################################
