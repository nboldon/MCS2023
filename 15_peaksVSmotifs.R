library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)


#Additional setup
setwd("/Volumes/DataBox/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE)

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

head(projMCS7$treatment)

############################################
############################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = C3ArchRSubset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks


features <- getMarkers(
  markersPeaks, 
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", 
  returnGR = TRUE
)
features


############################################
############################################
############################################
############################################
############################################
############################################

# Browser Tracks

p <- plotBrowserTrack(
  ArchRProj = C3ArchRSubset,
  groupBy = "treatment",
  region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
  features =  features,
  plotSummary = c("bulkTrack", "scTrack", "featureTrack", "geneTrack"),
  #loops = getCoAccessibility(C3ArchRSubset),
  useMatrix = "PeakMatrix",
  upstream = 5000,
  downstream = 5000
)
p

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-Tracks_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


#################

# Browser Tracks

p <- plotBrowserTrack(
  ArchRProj = C3ArchRSubset,
  groupBy = "treatment",
  region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
  plotSummary = c("scTrack", "featureTrack", "geneTrack"),
  loops = getCoAccessibility(projMCS7),
  useMatrix = "PeakMatrix",
  upstream = 5000,
  downstream = 5000
)
p

#"bulkTrack", removed; these are the usual browser tracks seen.

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-noBulkTrack_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

#################

p <- plotBrowserTrack(
  ArchRProj = projMCS7,
  groupBy = "treatment",
  region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
  plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
  loops = getCoAccessibility(projMCS7),
  useMatrix = "PeakMatrix",
  upstream = 5000,
  downstream = 5000
)
p

#"scTrack", removed; this is the mini duplicate of the bulkTrack

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-noScTrack_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


#################

# Browser Tracks

p <- plotBrowserTrack(
  ArchRProj = projMCS7,
  groupBy = "treatment",
  region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
  plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
  loops = getCoAccessibility(projMCS7),
  useMatrix = "PeakMatrix",
  upstream = 5000,
  downstream = 5000
)
p

#"featureTrack", removed; removed multicolored "peaks" plot that doesn't appear informative

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-noFeatureTrack_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


#################

p <- plotBrowserTrack(
  ArchRProj = projMCS7,
  groupBy = "treatment",
  region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
  plotSummary = c("bulkTrack", "scTrack", "featureTrack"),
  loops = getCoAccessibility(projMCS7),
  useMatrix = "PeakMatrix",
  upstream = 5000,
  downstream = 5000
)
p

# "geneTrack" removed; missing the "gene" section of the plot

# Save
plotPDF(p, name = "Peak-C3-T1vsT3-noGeneTrack_2024-09-16", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


############################################
############################################
############################################
############################################
############################################
############################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = projMCS7,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

########################

markersGenes <- getMarkerFeatures(
  ArchRProj = projMCS7,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersGenes


markerGS <- getMarkers(markersGenes, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")


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

head(projMCS7$treatment)

############################################
############################################


## Loop for all clusters, one region of interest

## T1vsT3_C3


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 97322494, end = 97322994)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T1vsT3-C3_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T4vsT1_C24


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr9", IRanges(start = 43270867, end = 43271367)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T4vsT1_C24_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T4vsT1-1_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 17576434, end = 17576934)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T4vsT1-1_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################

## T4vsT1-2_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr17", IRanges(start = 7169635, end = 7170135)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T4vsT1-2_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T3vsT1-1_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 91044363, end = 91044863)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T3vsT1-1_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T3vsT1-2_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 94411074, end = 94411574)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T3vsT1-2_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T3vsT1-3_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr17", IRanges(start = 7169635, end = 7170135)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T3vsT1-3_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T3vsT1-4_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 90143656, end = 90144156)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T3vsT1-4_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T1vsT3_C20


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr4", IRanges(start = 135831661, end = 135832161)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T1vsT3_C20_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T1vsT4_C18


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 91728462, end = 91728962)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T1vsT4_C18_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T3vsT2_C17


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr16", IRanges(start = 84834699, end = 84835199)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T3vsT2_C17_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}

######################################
######################################
######################################


## T4vsT3_C10


# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each cluster
for (cluster in clusters) {
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% cluster]
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges("chr14", IRanges(start = 27562549, end = 27563049)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE),
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_T4vsT3_C10_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name if there is an error
    message(paste("No result for cluster:", cluster))
  })
}


######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################


## Obtaining nfrag counts for regions of interest


## Loop for all samples & multiple regions of interest

# List of genomic regions of interest
regions <- list(
  c3_T1vsT3 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97322494, end = 97322994)),
  c24_T4vsT1 = GRanges(seqnames = "chr9", ranges = IRanges(start = 43270867, end = 43271367)),
  c20_T4vsT1_1 = GRanges(seqnames = "chr16", ranges = IRanges(start = 17576434, end = 17576934)),
  c20_T4vsT1_2 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7169635, end = 7170135)),
  c20_T3vsT1_1 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91044363, end = 91044863)),
  c20_T3vsT1_2 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94411074, end = 94411574)),
  c20_T3vsT1_3 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7169635, end = 7170135)),
  c20_T3vsT1_4 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90143656, end = 90144156)),
  c20_T1vsT3 = GRanges(seqnames = "chr4", ranges = IRanges(start = 135831661, end = 135832161)),
  c18_T1vsT4 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91728462, end = 91728962)),
  c17_T3vsT2 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834699, end = 84835199)),
  c10_T2vsT3 = GRanges(seqnames = "chr14", ranges = IRanges(start = 27562549, end = 27563049))
)


###########################################
###########################################

# Get the unique sample names from the project
sample_names <- unique(projMCS7$Sample)

# Loop over each genomic region
for (region_name in names(regions)) {
  
  # Initialize an empty data frame to store the results for the current region
  results <- data.frame(Sample = character(), stringsAsFactors = FALSE)
  
  # Get the current genomic region from the list
  region <- regions[[region_name]]
  
  # Loop over each sample
  for (i in sample_names) {
    # Get the indices of the cells corresponding to the current sample
    idxSample <- BiocGenerics::which(projMCS7$Sample %in% i)
    cellsSample <- projMCS7$cellNames[idxSample]
    
    # Subset the project for the current sample
    projSample <- projMCS7[cellsSample, ]
    
    # Extract the peak matrix for the current sample
    peakMatrix <- getMatrixFromProject(
      ArchRProj = projSample,
      useMatrix = "PeakMatrix",
      verbose = FALSE,
      binarize = FALSE,
      threads = getArchRThreads(),
      logFile = createLogFile(paste0("getMatrixFromProject_", i))
    )
    
    # Get the row ranges (genomic coordinates) of the genes
    peakRanges <- rowRanges(peakMatrix)
    
    # Subset the peak ranges by the genomic region of interest
    subset_peaks <- subsetByOverlaps(peakRanges, region)
    
    # Get the indices of the overlapping peaks
    indices <- subset_peaks$idx
    
    # Subset the peak matrix using these indices
    subset_fragments <- peakMatrix[indices, , drop = FALSE]
    
    # Calculate the total fragments for each cell
    total_fragments <- colSums(assay(subset_fragments))
    
    # Sum the total fragments across all cells for the current sample
    sum_fragments <- sum(total_fragments)
    
    # Add the results to the data frame (adding 'TotalFragments' as the column name)
    results <- rbind(results, data.frame(Sample = i, sum_fragments))
  }
  
  # Rename the column 'sum_fragments' to reflect the region name (e.g., "appTotalFrags")
  colnames(results)[2] <- paste0(region_name, "peakFrags")
  
  # Write the results to a CSV file named after the region
  write.csv(results, paste0(region_name, "_peakFrag_counts_2024-09-17.csv"), row.names = FALSE)
}


######################################
######################################
######################################
######################################


# List of CSV files to combine (assuming they're all in the current working directory)
csv_files <- list.files(pattern = "*_peakFrag_counts_2024-09-17.csv")

# Initialize an empty data frame
combined_results <- NULL

# Loop over each file and merge the data
for (csv_file in csv_files) {
  
  # Read the current CSV file
  current_data <- read.csv(csv_file)
  
  # If this is the first file, initialize the combined_results with the current data
  if (is.null(combined_results)) {
    combined_results <- current_data
  } else {
    # Merge the current data with the combined results on the "Sample" column
    combined_results <- merge(combined_results, current_data, by = "Sample", all = TRUE)
  }
}

# Create a treatment column based on the Sample column
combined_results$Treatment <- combined_results$Sample

# Group the samples into treatments using gsub
combined_results$Treatment <- gsub("C302_|C306_|C309_|C318_|C323_|C328_|C332_|C337_|C339_|C346_|C351_|C353_|C360_", "t1", combined_results$Treatment)
combined_results$Treatment <- gsub("C304_|C308_|C312_|C349_|C315_|C321_|C324_|C355_|C327_|C330_|C333_|C358_|C336_|C342_|C348_|C362_", "t2", combined_results$Treatment)
combined_results$Treatment <- gsub("C305_|C307_|C313_|C350_|C316_|C320_|C322_|C352_|C325_|C334_|C359_|C340_|C341_|C345_|C364_", "t3", combined_results$Treatment)
combined_results$Treatment <- gsub("C301_|C303_|C310_|C314_|C319_|C335_|C338_|C344_|C354_|C356_|C361_|C363_", "t4", combined_results$Treatment)

# Reorder the data frame by Treatment
combined_results <- combined_results[order(combined_results$Treatment), ]

# Write the combined and reordered data to a new CSV file
write.csv(combined_results, "combined_peakFrag_counts_byTx_2024-09-17.csv", row.names = FALSE)


######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################


## Combine individual files into one and normalize frag counts
# Normalize by dividing the number of gene frags by total number of frags in sample


# List of fragment count files
fragment_files <- list.files(pattern = "*_peakFrag_counts_2024-09-17.csv")

# Read the sample statistics file
sample_stats <- read.csv("sample_statistics.csv")

# Add treatment information based on the Sample column
sample_stats$Treatment <- sample_stats$Sample

# Group the samples into treatments using gsub
sample_stats$Treatment <- gsub("C302_|C306_|C309_|C318_|C323_|C328_|C332_|C337_|C339_|C346_|C351_|C353_|C360_", "t1", sample_stats$Treatment)
sample_stats$Treatment <- gsub("C304_|C308_|C312_|C349_|C315_|C321_|C324_|C355_|C327_|C330_|C333_|C358_|C336_|C342_|C348_|C362_", "t2", sample_stats$Treatment)
sample_stats$Treatment <- gsub("C305_|C307_|C313_|C350_|C316_|C320_|C322_|C352_|C325_|C334_|C359_|C340_|C341_|C345_|C364_", "t3", sample_stats$Treatment)
sample_stats$Treatment <- gsub("C301_|C303_|C310_|C314_|C319_|C335_|C338_|C344_|C354_|C356_|C361_|C363_", "t4", sample_stats$Treatment)

# Loop through each fragment count file and merge it with sample_stats
for (file in fragment_files) {
  # Read the fragment counts file
  fragment_counts <- read.csv(file)
  
  # Get the column name for fragment counts
  colname <- names(fragment_counts)[2] # Assuming the second column holds the fragment counts
  
  # Merge the fragment counts with the sample statistics
  sample_stats <- merge(sample_stats, fragment_counts, by = "Sample", all.x = TRUE)
  
  # Calculate the ratio for the current fragment counts and add a new column
  ratio_colname <- paste0(colname, "_Ratio")
  sample_stats[[ratio_colname]] <- sample_stats[[colname]] / sample_stats$TotalFragments
}

# Reorder the data frame by Treatment
sample_stats <- sample_stats[order(sample_stats$Treatment), ]

# Save the final merged data to a new CSV file
write.csv(sample_stats, "combined_normPeaks_statistics_with_treatment_2024-09-17.csv", row.names = FALSE)


######################################################
######################################################
######################################################





######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################

## The Ultimate loop for peak regions of interest looped through each cluster
########################
## THIS DOES NOT RUN ##
#Saved files only include C20 T2vsT3 names


# Create a data frame containing your regions
peak_regions <- data.frame(
  Chromosome = c("chr16", "chr9", "chr16", "chr17", "chr16", "chr16", "chr17", "chr16", "chr17", "chr17", "chr17", "chr17", "chr4", "chr16", "chr16", "chr16", "chr16", "chr14", "chr14"),
  Start = c(97322494, 43270867, 17576434, 7169635, 91044363, 94411074, 7169635, 90143656, 7170135, 7169635, 7169635, 7169635, 135831661, 91728462, 90143656, 84834699, 84834699, 27562549, 27562549),
  End = c(97322994, 43271367, 17576934, 7170135, 91044863, 94411574, 7170135, 90144156, 7170135, 7170135, 7170135, 7170135, 135832161, 91728962, 90144156, 84835199, 84835199, 27563049, 27563049),
  Log2FC = c(-2.35, -2.24, -1.13, 1.06, 1.02, 1.81, 0.96, 0.96, -1.02, -1.08, -0.93, -1.15, 1.50, -0.96, -0.95, 0.85, -0.78, -6.22, -2.73),
  FDR = c(9.44e-05, 0.0844, 0.0145, 0.0145, 0.0593, 0.0593, 0.0651, 0.0665, 0.0787, 0.00176, 0.00821, 5.01e-06, 0.0173, 0.00608, 0.0129, 0.0568, 0.0573, 0.0873, 0.0381)
)


#Cluster = c("C3", "C24", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C18", "C18", "C17", "C17", "C10", "C10"),
#Comparison = c("T1vsT3", "T4vsT1", "T4vsT1", "T4vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T2vsT4", "T2vsT3", "T1vsT4", "T1vsT3", "T1vsT3", "T1vsT4", "T1vsT4", "T3vsT2", "T2vsT3", "T4vsT3", "T2vsT3"),

# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each row in the data frame
for (i in 1:nrow(peak_regions)) {
  
  # Extract the region and cluster information
  chrom <- peak_regions$Chromosome[i]
  start <- peak_regions$Start[i]
  end <- peak_regions$End[i]
  
  #cluster <- peak_regions$Cluster[i]
  #comparison <- peak_regions$Comparison[i]
  
  # Subset the ArchR project by the specific cluster for this iteration
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% clusters]
  
  # Flag to check if plot was created
  plot_created <- FALSE
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster and region
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  
      region = GRanges(chrom, IRanges(start = start, end = end)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", clusters, "_", comparison, "_", chrom, "_", start, "_", end, "_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
    # Set flag to TRUE since plot was created
    plot_created <- TRUE
    
  }, error = function(e) {
    # Print the cluster name and comparison if there is an error
    message(paste("Error for cluster:", clusters, "comparison:", comparison, "-", e$message))
  })
  
  # Print message if no plot was created
  if (!plot_created) {
    message(paste("No result for cluster:", clusters, "comparison:", comparison))
  }
}



######################################
######################################
######################################
######################################
######################################
######################################


## Normalizing by nFrags
########################
## THIS DOES NOT RUN ##
#Error in normalizeByNFrags(subsetArchRProj, normFactor) : 
#no slot of name "assays" for this object of class "ArchRProject"

# Create a data frame containing your regions
peak_regions <- data.frame(
  Cluster = c("C3", "C24", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C18", "C18", "C17", "C17", "C10", "C10"),
  Comparison = c("T1vsT3", "T4vsT1", "T4vsT1", "T4vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T2vsT4", "T2vsT3", "T1vsT4", "T1vsT3", "T1vsT3", "T1vsT4", "T1vsT4", "T3vsT2", "T2vsT3", "T4vsT3", "T2vsT3"),
  Chromosome = c("chr16", "chr9", "chr16", "chr17", "chr16", "chr16", "chr17", "chr16", "chr17", "chr17", "chr17", "chr17", "chr4", "chr16", "chr16", "chr16", "chr16", "chr14", "chr14"),
  Start = c(97322494, 43270867, 17576434, 7169635, 91044363, 94411074, 7169635, 90143656, 7170135, 7169635, 7169635, 7169635, 135831661, 91728462, 90143656, 84834699, 84834699, 27562549, 27562549),
  End = c(97322994, 43271367, 17576934, 7170135, 91044863, 94411574, 7170135, 90144156, 7170135, 7170135, 7170135, 7170135, 135832161, 91728962, 90144156, 84835199, 84835199, 27563049, 27563049),
  Log2FC = c(-2.35, -2.24, -1.13, 1.06, 1.02, 1.81, 0.96, 0.96, -1.02, -1.08, -0.93, -1.15, 1.50, -0.96, -0.95, 0.85, -0.78, -6.22, -2.73),
  FDR = c(9.44e-05, 0.0844, 0.0145, 0.0145, 0.0593, 0.0593, 0.0651, 0.0665, 0.0787, 0.00176, 0.00821, 5.01e-06, 0.0173, 0.00608, 0.0129, 0.0568, 0.0573, 0.0873, 0.0381)
)

# Loop through each row in the data frame
for (i in 1:nrow(peak_regions)) {
  
  # Extract the cluster and region information
  cluster <- peak_regions$Cluster[i]
  comparison <- peak_regions$Comparison[i]
  chrom <- peak_regions$Chromosome[i]
  start <- peak_regions$Start[i]
  end <- peak_regions$End[i]
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7@cellColData$Clusters == cluster]
  
  # Extract "nFrags" for normalization
  nFrags <- subsetArchRProj@cellColData$nFrags
  
  # Normalize the matrix by "nFrags"
  normFactor <- nFrags / median(nFrags)  # Normalize to median
  
  # You need to preprocess the data or create a custom matrix with normalization applied before plotting
  # Ensure you are using a preprocessed matrix for the normalization
  # e.g., customMatrix <- yourNormalizationFunction(subsetArchRProj@peakMatrix, normFactor)
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster and region
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges(chrom, IRanges(start = start, end = end)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",  # Adjust if needed based on the matrix type
      upstream = 5000,
      downstream = 5000
      # The normalization needs to be applied before plotting
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_", comparison, "_Tracks_nFragsNorm_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name and comparison if there is an error
    message(paste("No result for cluster:", cluster, "comparison:", comparison))
  })
}



######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################



## The Ultimate loop for peak regions of interest looped through each cluster
########################
## THIS DOES NOT RUN ##
#Saved files only include C20 T2vsT3 names


# Create a data frame containing your regions
peak_regions <- data.frame(
  Chromosome = c("chr16", "chr9", "chr16", "chr17", "chr16", "chr16", "chr17", "chr16", "chr17", "chr17", "chr17", "chr17", "chr4", "chr16", "chr16", "chr16", "chr16", "chr14", "chr14"),
  Start = c(97322494, 43270867, 17576434, 7169635, 91044363, 94411074, 7169635, 90143656, 7170135, 7169635, 7169635, 7169635, 135831661, 91728462, 90143656, 84834699, 84834699, 27562549, 27562549),
  End = c(97322994, 43271367, 17576934, 7170135, 91044863, 94411574, 7170135, 90144156, 7170135, 7170135, 7170135, 7170135, 135832161, 91728962, 90144156, 84835199, 84835199, 27563049, 27563049),
  Log2FC = c(-2.35, -2.24, -1.13, 1.06, 1.02, 1.81, 0.96, 0.96, -1.02, -1.08, -0.93, -1.15, 1.50, -0.96, -0.95, 0.85, -0.78, -6.22, -2.73),
  FDR = c(9.44e-05, 0.0844, 0.0145, 0.0145, 0.0593, 0.0593, 0.0651, 0.0665, 0.0787, 0.00176, 0.00821, 5.01e-06, 0.0173, 0.00608, 0.0129, 0.0568, 0.0573, 0.0873, 0.0381)
)


#Cluster = c("C3", "C24", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C18", "C18", "C17", "C17", "C10", "C10"),
#Comparison = c("T1vsT3", "T4vsT1", "T4vsT1", "T4vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T2vsT4", "T2vsT3", "T1vsT4", "T1vsT3", "T1vsT3", "T1vsT4", "T1vsT4", "T3vsT2", "T2vsT3", "T4vsT3", "T2vsT3"),

# Define the clusters you want to loop through
clusters <- unique(projMCS7$Clusters)

# Loop through each row in the data frame
for (i in 1:nrow(peak_regions)) {
  
  # Extract the region and cluster information
  chrom <- peak_regions$Chromosome[i]
  start <- peak_regions$Start[i]
  end <- peak_regions$End[i]
  
  #cluster <- peak_regions$Cluster[i]
  #comparison <- peak_regions$Comparison[i]
  
  # Subset the ArchR project by the specific cluster for this iteration
  subsetArchRProj <- projMCS7[projMCS7$Clusters %in% clusters]
  
  # Flag to check if plot was created
  plot_created <- FALSE
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster and region
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  
      region = GRanges(chrom, IRanges(start = start, end = end)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",
      upstream = 5000,
      downstream = 5000
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", clusters, "_", comparison, "_", chrom, "_", start, "_", end, "_Tracks_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
    # Set flag to TRUE since plot was created
    plot_created <- TRUE
    
  }, error = function(e) {
    # Print the cluster name and comparison if there is an error
    message(paste("Error for cluster:", clusters, "comparison:", comparison, "-", e$message))
  })
  
  # Print message if no plot was created
  if (!plot_created) {
    message(paste("No result for cluster:", clusters, "comparison:", comparison))
  }
}



######################################
######################################
######################################
######################################
######################################
######################################


## Normalizing by nFrags
########################
## THIS DOES NOT RUN ##
#Error in normalizeByNFrags(subsetArchRProj, normFactor) : 
#no slot of name "assays" for this object of class "ArchRProject"

# Create a data frame containing your regions
peak_regions <- data.frame(
  Cluster = c("C3", "C24", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C20", "C18", "C18", "C17", "C17", "C10", "C10"),
  Comparison = c("T1vsT3", "T4vsT1", "T4vsT1", "T4vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T3vsT1", "T2vsT4", "T2vsT3", "T1vsT4", "T1vsT3", "T1vsT3", "T1vsT4", "T1vsT4", "T3vsT2", "T2vsT3", "T4vsT3", "T2vsT3"),
  Chromosome = c("chr16", "chr9", "chr16", "chr17", "chr16", "chr16", "chr17", "chr16", "chr17", "chr17", "chr17", "chr17", "chr4", "chr16", "chr16", "chr16", "chr16", "chr14", "chr14"),
  Start = c(97322494, 43270867, 17576434, 7169635, 91044363, 94411074, 7169635, 90143656, 7170135, 7169635, 7169635, 7169635, 135831661, 91728462, 90143656, 84834699, 84834699, 27562549, 27562549),
  End = c(97322994, 43271367, 17576934, 7170135, 91044863, 94411574, 7170135, 90144156, 7170135, 7170135, 7170135, 7170135, 135832161, 91728962, 90144156, 84835199, 84835199, 27563049, 27563049),
  Log2FC = c(-2.35, -2.24, -1.13, 1.06, 1.02, 1.81, 0.96, 0.96, -1.02, -1.08, -0.93, -1.15, 1.50, -0.96, -0.95, 0.85, -0.78, -6.22, -2.73),
  FDR = c(9.44e-05, 0.0844, 0.0145, 0.0145, 0.0593, 0.0593, 0.0651, 0.0665, 0.0787, 0.00176, 0.00821, 5.01e-06, 0.0173, 0.00608, 0.0129, 0.0568, 0.0573, 0.0873, 0.0381)
)

# Loop through each row in the data frame
for (i in 1:nrow(peak_regions)) {
  
  # Extract the cluster and region information
  cluster <- peak_regions$Cluster[i]
  comparison <- peak_regions$Comparison[i]
  chrom <- peak_regions$Chromosome[i]
  start <- peak_regions$Start[i]
  end <- peak_regions$End[i]
  
  # Subset the ArchR project by cluster
  subsetArchRProj <- projMCS7[projMCS7@cellColData$Clusters == cluster]
  
  # Extract "nFrags" for normalization
  nFrags <- subsetArchRProj@cellColData$nFrags
  
  # Normalize the matrix by "nFrags"
  normFactor <- nFrags / median(nFrags)  # Normalize to median
  
  # You need to preprocess the data or create a custom matrix with normalization applied before plotting
  # Ensure you are using a preprocessed matrix for the normalization
  # e.g., customMatrix <- yourNormalizationFunction(subsetArchRProj@peakMatrix, normFactor)
  
  # Run plotBrowserTrack inside a try-catch to handle errors
  tryCatch({
    # Generate browser track for the current cluster and region
    p <- plotBrowserTrack(
      ArchRProj = subsetArchRProj,
      groupBy = "treatment",  # Ensure you have the correct metadata field
      region = GRanges(chrom, IRanges(start = start, end = end)),
      plotSummary = c("bulkTrack", "geneTrack"),
      useMatrix = "PeakMatrix",  # Adjust if needed based on the matrix type
      upstream = 5000,
      downstream = 5000
      # The normalization needs to be applied before plotting
    )
    
    # If the plot is successfully created, save it as a PDF
    plotPDF(p, name = paste0("Peak-", cluster, "_", comparison, "_Tracks_nFragsNorm_2024-09-17"), 
            width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)
    
  }, error = function(e) {
    # Print the cluster name and comparison if there is an error
    message(paste("No result for cluster:", cluster, "comparison:", comparison))
  })
}




######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################

######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################



sessionInfo()

R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.1.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

Random number generation:
  RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Denver
tzcode source: internal

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] presto_1.0.0                ggrepel_0.9.6               nabor_0.5.0                 pheatmap_1.0.12            
 [5] enrichplot_1.24.4           clusterProfiler_4.12.6      org.Mm.eg.db_3.19.1         AnnotationDbi_1.66.0       
 [9] rhdf5_2.48.0                SummarizedExperiment_1.34.0 Biobase_2.64.0              MatrixGenerics_1.16.0      
[13] Rcpp_1.0.13                 Matrix_1.7-0                GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
[17] IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.50.0         matrixStats_1.4.1          
[21] data.table_1.16.0           stringr_1.5.1               plyr_1.8.9                  magrittr_2.0.3             
[25] ggplot2_3.5.1               gtable_0.3.5                gtools_3.9.5                gridExtra_2.3              
[29] ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3                 rstudioapi_0.16.0                  jsonlite_1.8.8                    
  [4] shape_1.4.6.1                      farver_2.1.2                       rmarkdown_2.28                    
  [7] BiocIO_1.14.0                      GlobalOptions_0.1.2                fs_1.6.4                          
 [10] zlibbioc_1.50.0                    vctrs_0.6.5                        Rsamtools_2.20.0                  
 [13] Cairo_1.6-2                        memoise_2.0.1                      RCurl_1.98-1.16                   
 [16] ggtree_3.12.0                      S4Arrays_1.4.1                     htmltools_0.5.8.1                 
 [19] curl_5.2.2                         Rhdf5lib_1.26.0                    SparseArray_1.4.8                 
 [22] gridGraphics_0.5-1                 httr2_1.0.3                        cachem_1.1.0                      
 [25] GenomicAlignments_1.40.0           igraph_2.0.3                       lifecycle_1.0.4                   
 [28] iterators_1.0.14                   pkgconfig_2.0.3                    R6_2.5.1                          
 [31] fastmap_1.2.0                      gson_0.1.0                         GenomeInfoDbData_1.2.12           
 [34] clue_0.3-65                        digest_0.6.37                      aplot_0.2.3                       
 [37] colorspace_2.1-1                   patchwork_1.2.0                    RSQLite_2.3.7                     
 [40] labeling_0.4.3                     fansi_1.0.6                        abind_1.4-5                       
 [43] httr_1.4.7                         polyclip_1.10-7                    compiler_4.4.0                    
 [46] bit64_4.0.5                        withr_3.0.1                        doParallel_1.0.17                 
 [49] BiocParallel_1.38.0                viridis_0.6.5                      DBI_1.2.3                         
 [52] ggforce_0.4.2                      R.utils_2.12.3                     MASS_7.3-61                       
 [55] DelayedArray_0.30.1                rappdirs_0.3.3                     rjson_0.2.22                      
 [58] tools_4.4.0                        ape_5.8                            scatterpie_0.2.4                  
 [61] R.oo_1.26.0                        glue_1.7.0                         restfulr_0.0.15                   
 [64] rhdf5filters_1.16.0                nlme_3.1-166                       GOSemSim_2.30.2                   
 [67] shadowtext_0.1.4                   cluster_2.1.6                      reshape2_1.4.4                    
 [70] fgsea_1.30.0                       generics_0.1.3                     BSgenome_1.72.0                   
 [73] R.methodsS3_1.8.2                  tidyr_1.3.1                        tidygraph_1.3.1                   
 [76] utf8_1.2.4                         XVector_0.44.0                     foreach_1.5.2                     
 [79] pillar_1.9.0                       yulab.utils_0.1.7                  circlize_0.4.16                   
 [82] splines_4.4.0                      dplyr_1.1.4                        tweenr_2.0.3                      
 [85] treeio_1.28.0                      lattice_0.22-6                     rtracklayer_1.64.0                
 [88] bit_4.0.5                          tidyselect_1.2.1                   GO.db_3.19.1                      
 [91] ComplexHeatmap_2.20.0              Biostrings_2.72.1                  knitr_1.48                        
 [94] xfun_0.47                          graphlayouts_1.1.1                 stringi_1.8.4                     
 [97] UCSC.utils_1.0.0                   lazyeval_0.2.2                     ggfun_0.1.6                       
[100] yaml_2.3.10                        evaluate_0.24.0                    codetools_0.2-20                  
[103] BSgenome.Mmusculus.UCSC.mm10_1.4.3 ggraph_2.2.1                       tibble_3.2.1                      
[106] qvalue_2.36.0                      BiocManager_1.30.25                ggplotify_0.1.2                   
[109] cli_3.6.3                          munsell_0.5.1                      png_0.1-8                         
[112] XML_3.99-0.17                      blob_1.2.4                         DOSE_3.30.5                       
[115] bitops_1.0-8                       viridisLite_0.4.2                  tidytree_0.4.6                    
[118] scales_1.3.0                       purrr_1.0.2                        crayon_1.5.3                      
[121] GetoptLong_1.0.5                   rlang_1.1.4                        cowplot_1.1.3                     
[124] fastmatch_1.1-4                    KEGGREST_1.44.1 
