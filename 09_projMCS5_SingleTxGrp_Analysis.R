#Setup an interactive session
salloc --account=eon -t 0-16:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

############################################
############################################

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

############################################
############################################

# To compare samples amongst tx groups in a particular cluster:
t1ArchRSubset <-projMCS5[
	projMCS5$Sample %in% c("C302_", "C306_", "C309_", "C318_", "C328_", "C332_", "C337_", "C339_",
	"C346_", "C351_", "C353_", "C360_"),]

t1ArchRSubset

############################################
############################################

# Specify which treatment group each sample is in:

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- projMCS5$Sample

treatment <- gsub("C302_", "t1", treatment)
treatment <- gsub("C306_", "t1", treatment)
treatment <- gsub("C309_", "t1", treatment)
treatment <- gsub("C318_", "t1", treatment)
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
treatment <- gsub("C323_", "t3", treatment)
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
projMCS5$treatment <- treatment

# Check that this worked - if not, make sure the previous line was run successfully
head(projMCS5$treatment)

##########################################
##########################################

# To select only a single treatment to run things on:
projT1 <- projMCS5$treatment[which(projMCS5$treatment=="t1",)]
projt2 <- projMCS5$treatment[which(projMCS5$treatment=="t2",)]
projt3 <- projMCS5$treatment[which(projMCS5$treatment=="t3",)]
projt4 <- projMCS5$treatment[which(projMCS5$treatment=="t4",)]
## DOES NOT create a workable project

## OR
proj2N <- projMCS5$treatment[which(projMCS5$treatment %in% c("t1", "t2")),]
# Error in projMCS5$treatment[which(projMCS5$treatment %in% c("t1", "t2")),  : 
#  incorrect number of dimensions


##########################################
##########################################
##########################################
##########################################

## Started running code here 1-31-2024:

# To compare samples amongst tx groups in a particular cluster:
t1ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C302_", "C306_", "C309_", "C318_", "C328_", "C332_", "C337_", "C339_",
        "C346_", "C351_", "C353_", "C360_"),]

t1ArchRSubset


#To make comparisons and plots among specific treatments; use getMarkerFeatures to identify differentially accessible peaks in two different treatments using "treatment" for the "groupby" argument and then specify the two treatments to compare in the "useGroups" and "bgdGroups" arguments

t1markerTest <- getMarkerFeatures(
   ArchRProj = t1ArchRSubset,
   useMatrix = "GeneScoreMatrix",
   groupBy = "Sample", 
   testMethod = "wilcoxon",
   bias = c("TSSEnrichment",
   "log10(nFrags)"), 
)

t1markerTest

--------------

table(t1ArchRSubset$Clusters)

markersGenes <- getMarkerFeatures(
    ArchRProj = t1ArchRSubset, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersGenes

markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", #returnGR = TRUE)
)
markerList

markerList$C1

#markerList$t1ArchRSubset
#NULL

heatmapGenes <- markerHeatmap(
  seMarker = markersGenes, 
  cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25",
  transpose = TRUE
)

draw(heatmapGenes, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGenes, name = "T1_GeneMarker_Heatmap.pdf", width = 8, height = 6, ArchRProj = t1ArchRSubset, addDOC = FALSE)


###############################
###############################
