#Setup an interactive session
salloc --account=eon -t 1-00:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE)

getAvailableMatrices(projMCS7)


############################################
############################################
############################################
############################################

markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS7,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

############################################

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)
markerList


############################################
############################################
############################################
############################################                                            
############################################
############################################


# Motif enrichment in marker peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projMCS7,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
  )

enrichMotifs

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_2024-09-16", width = 8, height = 6, ArchRProj = projMCS7, addDOC = FALSE)


############################################
############################################

## Motif Deviations

# First check to make sure we've added motif annotations to the ArchRProject
if("Motif" %ni% names(projMCS7@peakAnnotation)){c
    projMCS7 <- addMotifAnnotations(ArchRProj = projMCS7, motifSet = "cisbp", name = "Motif")
}

# We also need to add a set of background peaks which are used in computing deviations. 
# Background peaks are chosen using the chromVAR::getBackgroundPeaks() function which samples peaks based on similarity in GC-content and$
projMCS7 <- addBgdPeaks(projMCS7)
## 2024-09-16: Background peaks already exist

# We are now ready to compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. 
# This function has an optional parameter called matrixName that allows us to define the name of the deviations matrix that will be stored in the Arrow files. 
# If we do not provide a value to this parameter, as in the example below, this function creates a matrix name by adding the word “Matrix" to the name of the peakAnnotation. 
# The example below creates a deviations matrix in each of our Arrow files called “MotifMatrix”.
projMCS7 <- addDeviationsMatrix(
  ArchRProj = projMCS7,
  peakAnnotation = "Motif",
  force = TRUE
)

# To access these deviations, we use the getVarDeviations() function. 
# If we want this function to return a ggplot object, we set plot = TRUE otherwise, this function would return the DataFrame object. 
# The head of that DataFrame object is displayed by default when the function is run.
plotVarDev <- getVarDeviations(projMCS7, name = "MotifMatrix", plot = TRUE)

# Plot variable deviations
plotVarDev

# Save plot
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


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
      useMatrix = "MotifMatrix",
      groupBy = "treatment",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = useGroup,
      bgdGroups = bgdGroup
    )

    # Get the list of markers
    markerList <- getMarkers(markerFeatures, cutOff = "FDR <= 0.1")

    # Define file name
    file_name <- paste0(cluster, "_", useGroup, "vs", bgdGroup, "_Motif_FDR-0-1_2024-09-16.csv")

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
############################################
############################################
############################################


# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

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

############################################
############################################
############################################
############################################
############################################
############################################


# Browser Tracks

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
plotPDF(p, name = "C18-Tracks-With-PeakFeatures_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALS$


############################################
############################################


# Browser Tracks

p <- plotBrowserTrack(
    ArchRProj = projMCS7,
    groupBy = "treatment",
    region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
    plotSummary = c("bulkTrack", "scTrack", "featureTrack", "geneTrack"),
    loops = getCoAccessibility(projMCS7),   
    useMatrix = "PeakMatrix",
    upstream = 5000,
    downstream = 5000
)

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-Tracks_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

############################################
############################################
############################################
############################################                                            
############################################
############################################


sessionInfo()
