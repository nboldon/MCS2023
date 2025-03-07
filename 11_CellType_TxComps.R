
## Code did not produce meaningful results


#Setup an interactive session
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
library(tidyr)
library(cowplot)
library(ggplot2)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ProjMCS6/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", 
                             force = FALSE, showLogo = FALSE)
projMCS6
#NumberOfCells: 104455
#medianTSS: 14.521
#medianFrags: 6653

getAvailableMatrices(projMCS6)
# "GeneScoreMatrix" "PeakMatrix"      "TileMatrix" 

#############################################
#############################################

treatment <- projMCS6$Sample
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
projMCS6$treatment <- treatment
# Check that this worked
head(projMCS6$treatment)

#############################################

celltype <- projMCS6$Clusters
celltype[celltype=="C1"] <- "astrocyte"
celltype[celltype=="C2"] <- "oligodendrocyte"
celltype[celltype=="C3"] <- "oligodendrocyte"
celltype[celltype=="C4"] <- "glutamatergic"
celltype[celltype=="C5"] <- "oligodendrocyte"
celltype[celltype=="C6"] <- "oligodendrocyte_precursor"
celltype[celltype=="C7"] <- "astrocyte"
celltype[celltype=="C8"] <- "astrocyte"
celltype[celltype=="C9"] <- "astrocyte"
celltype[celltype=="C10"] <- "microglia"
celltype[celltype=="C11"] <- "microglia"
celltype[celltype=="C12"] <- "endo-vasc"
celltype[celltype=="C13"] <- "glutamatergic"
celltype[celltype=="C14"] <- "endo_vasc"
celltype[celltype=="C15"] <- "glutamatergic"
celltype[celltype=="C16"] <- "glutamatergic"
celltype[celltype=="C17"] <- "glutamatergic"
celltype[celltype=="C18"] <- "glutamatergic"
celltype[celltype=="C19"] <- "glutamatergic"
celltype[celltype=="C20"] <- "glutamatergic"
celltype[celltype=="C21"] <- "glutamatergic"
celltype[celltype=="C22"] <- "GABAergic"
celltype[celltype=="C23"] <- "GABAergic"
celltype[celltype=="C24"] <- "glutamatergic"
celltype[celltype=="C25"] <- "glutamatergic"

#Check that the only entries in this object are the defined cell types now
unique(celltype)

list(celltype)

projMCS6$celltype <- celltype

##########################################

projMCS6$type_treat <- paste0(projMCS6$celltype, "_", projMCS6$treatment)
projMCS6$type_treat

##########################################

##########################################

## Identify marker peaks and account for differences in data quality amongst cell groups
markersGenes_type_treat <- getMarkerFeatures(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "type_treat",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
 
markersPeaks_type_treat

#To retrieve particular slices of the SummarizedExperiment; returns a list of DataFrame objects, 
#one for each cell group
markerList_type_treat <- getMarkers(markersPeaks_type_treat, cutOff = "FDR <= 0.01 & Log2FC >= 1")  
markerList_type_treat

#The GRangesList object can be subset to a GRanges object for a particular cell group
markerList_type_treat$microglia_t1
markerList_type_treat$microglia_t2
markerList_type_treat$microglia_t3
markerList_type_treat$microglia_t4
markerList_type_treat$oligodendrocyte_t1
markerList_type_treat$oligodendrocyte_t2
markerList_type_treat$oligodendrocyte_t3
markerList_type_treat$oligodendrocyte_t4
markerList_type_treat$oligodendrocyte_precursor_t1
markerList_type_treat$oligodendrocyte_precursor_t2
markerList_type_treat$oligodendrocyte_precursor_t3
markerList_type_treat$oligodendrocyte_precursor_t4
markerList_type_treat$endo_vasc_t1
markerList_type_treat$endo_vasc_t2
markerList_type_treat$endo_vasc_t3
markerList_type_treat$endo_vasc_t4
markerList_type_treat$glutamatergic_t1
markerList_type_treat$glutamatergic_t2
markerList_type_treat$glutamatergic_t3
markerList_type_treat$glutamatergic_t4
markerList_type_treat$interneuron_t1
markerList_type_treat$interneuron_t2
markerList_type_treat$interneuron_t3
markerList_type_treat$interneuron_t4
markerList_type_treat$GABAergic_t1
markerList_type_treat$GABAergic_t2
markerList_type_treat$GABAergic_t3
markerList_type_treat$GABAergic_t4

###############################################
###############################################
###############################################
