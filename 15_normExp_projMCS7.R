## Experiments for normalization

#######################################
#######################################

#Setup an interactive session
salloc --account=eon -t 0-03:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Updated conda env 12-2023
module load miniconda3/23.1.0
conda activate archr2023_12

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)


#Additional setup
setwd("/project/eon/nboldon/MCS2023/Seurat/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", 
                             force = FALSE, showLogo = FALSE)

projMCS7
getAvailableMatrices(projMCS7)

#######################################

markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS7, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

# Returns a list of dataframe objects
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
markerList

############################################

# Returns a GRangesList object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)
markerList

############################################
############################################
############################################
############################################

## Browser Tracks - Normalization experiments

##################

# Interesting genes

# Normalize by Sample:
p1 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Sample",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Sample",
    #ylim = c(0, 10)
)

plotPDF(p1, name = "normTracksbySample-Interesting-With-PeakFeatures_2024-08-28", width = 5, height = 10, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Interesting genes

# Normalize by Cluster:
p2 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Clusters",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Cluster",
    #ylim = c(0, 10)
)

plotPDF(p2, name = "normTracksbyClusters-Interesting-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

####################

# Chr16 genes

# Normalize by Sample:
p3 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Sample",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("App", "Bach1", "Cldn14", "Cldn17", "Dscam", "Ets2", "Grik1", "Hunk", "Kcnj6", "Olig1", "Olig2", "Runx1", "Tiam1", "Urb1"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Sample",
    #ylim = c(0, 10)
)

plotPDF(p3, name = "normTracksbySample-Chr16-With-PeakFeatures_2024-08-28", width = 5, height = 10, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Normalize by Sample:
p4 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Clusters",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("App", "Bach1", "Cldn14", "Cldn17", "Dscam", "Ets2", "Grik1", "Hunk", "Kcnj6", "Olig1", "Olig2", "Runx1", "Tiam1", "Urb1"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Sample",
    #ylim = c(0, 10)
)

plotPDF(p4, name = "normTracksbyCluster-Chr16-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)





############################################
############################################
############################################
############################################

# Browser Tracks From ArchR Manual 

p <- plotBrowserTrack(
    ArchRProj = projMCS7, 
    groupBy = "Clusters", 
    geneSymbol = c("GATA1", "Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)["C18"],
    upstream = 50000,
    downstream = 50000
)

# Plot browser tracks
grid::grid.newpage()
grid::grid.draw(p$GATA1)
grid::grid.draw(p$Slc17a7)
grid::grid.draw(p$Slc17a6)
grid::grid.draw(p$Slc17a8)
grid::grid.draw(p$Slc32a1)

# Save
plotPDF(p, name = "C18-Tracks-With-PeakFeatures_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

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


# Subset by Cluster
C1ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C1"]
C2ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C2"]
C3ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C3"]
C4ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C4"]
C5ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C5"]
C6ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C6"]
C7ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C7"]
C8ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C8"]
C9ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C9"]
C10ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C10"]
C11ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C11"]
C12ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C12"]
C13ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C13"]
C14ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C14"]
C15ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C15"]
C16ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C16"]
C17ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C17"]
C18ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C18"]
C19ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C19"]
C20ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C20"]
C21ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C21"]
C22ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C22"]
C23ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C23"]
C24ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C24"]
C25ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C25"]

############################################
############################################
############################################

## Cluster 1 

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

# Interesting genes

# Normalize by Sample:
p1 <- plotBrowserTrack(
    ArchRProj = C3ArchRSubset,
    groupBy = "treatment",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Olig1", "Olig2", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "C3 Normalized Fragment Counts by Treatment",
    #ylim = c(0, 10)
)

plotPDF(p1, name = "C3-normTracksbyTx-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Interesting genes

# Normalize by Cluster:
p2 <- plotBrowserTrack(
    ArchRProj = C3ArchRSubset,
    groupBy = "treatment",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Olig1", "Olig2", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "C3 Normalized Fragment Counts by Treatment",
    #ylim = c(0, 10)
)

plotPDF(p2, name = "C3-normTracksbyClusters-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

####################
####################
####################

# Interesting genes

# Normalize by Sample:
p1 <- plotBrowserTrack(
    ArchRProj = C18ArchRSubset,
    groupBy = "treatment",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "C18 Normalized Fragment Counts by Treatment",
    #ylim = c(0, 10)
)

plotPDF(p1, name = "C18-normTracksbyTx-Interesting-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Interesting genes

# Normalize by Cluster:
p2 <- plotBrowserTrack(
    ArchRProj = C18ArchRSubset,
    groupBy = "treatment",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin"),
    upstream = 2000,
    downstream = 2000,
    title = "C18 Normalized Fragment Counts by Treatment",
    #ylim = c(0, 10)
)

plotPDF(p2, name = "C18-normTracksbyTx-Interesting-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

####################


####################
####################
####################
####################
####################
####################
####################
####################
####################

#Below is duplicate code only

# Chr16 genes

# Normalize by Sample:
p3 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Sample",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("App", "Bach1", "Cldn14", "Cldn17", "Dscam", "Ets2", "Grik1", "Hunk", "Kcnj6", "Olig1", "Olig2", "Runx1", "Tiam1", "Urb1"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Sample",
    #ylim = c(0, 10)
)

plotPDF(p3, name = "normTracksbySample-Chr16-With-PeakFeatures_2024-08-28", width = 5, height = 10, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Normalize by Sample:
p4 <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "Clusters",
    #useGroups = c("C301", "C305"),
    geneSymbol = c("App", "Bach1", "Cldn14", "Cldn17", "Dscam", "Ets2", "Grik1", "Hunk", "Kcnj6", "Olig1", "Olig2", "Runx1", "Tiam1", "Urb1"),
    upstream = 2000,
    downstream = 2000,
    title = "Normalized Fragment Counts by Sample",
    #ylim = c(0, 10)
)

plotPDF(p4, name = "normTracksbyCluster-Chr16-With-PeakFeatures_2024-08-28", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

