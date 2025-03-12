

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


###############################################################
###############################################################
###############################################################


## Overlay cells from specific samples on the UMAP of all cells by using the plotEmbedding function


treatment <- projMCS7$Sample

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
projMCS7$treatment <- treatment

# Check that this worked - if not, make sure the previous line was run successfully
head(projMCS7$treatment)

##############################################################


## Colors cells by treatment group

plotEmbedding(
  ArchRProj = projMCS7,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "treatment",
  plotAs = "points"
)


##############################################################


## Colors cells of interest, leaving all other cells grey in the background


# Extract all cell names in the ArchR project
cellNames <- getCellNames(projMCS7)

# Extract cell names for t1 and t2 treatments
highlightCells <- cellNames[projMCS7$treatment %in% c("t3", "t4")]

# Plot with only t1 and t2 cells highlighted
plotEmbedding(
  ArchRProj = projMCS7,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "treatment",
  highlightCells = highlightCells, # Highlight specific cells
  plotAs = "points"
)


##############################################################


## Using custom colors for treatment groups

# Extract all cell names in the ArchR project
cellNames <- getCellNames(projMCS7)

# Extract cell names for t3 and t4 treatments
highlightCells <- cellNames[projMCS7$treatment %in% c("t3")]

# Plot with custom colors for t3 and t4 cells
plotEmbedding(
  ArchRProj = projMCS7,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "treatment",
  highlightCells = highlightCells,
  plotAs = "points",
  # Add custom color scheme
  pal = c(
    "t1" = "lightgrey",  # non-highlighted cells
    "t2" = "lightgrey",  # non-highlighted cells
    "t3" = "purple",  # highlighted cells
    "t4" = "lightgrey"      # highlighted cells
  )
)

