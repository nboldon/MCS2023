#Setup an interactive session
salloc --account=eon -t 0-08:00:00 --mem=256G --nodes=1 --ntasks-per-node=16

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
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

############################################
############################################

getPeakSet(projMCS4)

projMCS7 <- addPeakMatrix(projMCS4)

getAvailableMatrices(projMCS7)

saveArchRProject(ArchRProj = projMCS7, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", load = FALSE)

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

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
markerList

############################################

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)
markerList

############################################
############################################

# Heatmaps

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap_2024-08-07", width = 8, height = 6, ArchRProj = projMCS7, addDOC = FALSE)

############################################
############################################

# MA Plots

pma <- plotMarkers(seMarker = markersPeaks, name = "C18", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pma

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markersPeaks, name = "C18", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pma, pv, name = "C18-PeakMarkers-MA-Volcano_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

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
plotPDF(p, name = "C18-Tracks-With-PeakFeatures_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)


############

region <- GRanges(
  seqnames = "chr16",  # Replace with the correct chromosome
  ranges = IRanges(start = 84985335, end = 85173989)
)

p <- plotBrowserTrack(
  ArchRProj = projMCS7, 
  groupBy = "Sample", 
  region = region,
  normMethod = nFrags,
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)["C18"],
  upstream = 2000,
  downstream = 2000
)

plotPDF(p, name = "App-Tracks-With-PeakFeatures_2024-08-29", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

############################################
############################################

# Pairwise testing between groups

markerTest <- getMarkerFeatures(
  ArchRProj = projMCS7, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C18",
  bgdGroups = "C22"
)

# MA Plots

pma <- plotMarkers(seMarker = markerTest, name = "C18", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pma

# Volcano Plots

pv <- plotMarkers(seMarker = markerTest, name = "C18", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

# Print Plots

plotPDF(pma, pv, name = "C18-vs-C22-Markers-MA-Volcano_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

############################################
############################################

# Add motif annotations

projMCS7 <- addMotifAnnotations(ArchRProj = projMCS7, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest, #Completed earlier for C18
    ArchRProj = projMCS7,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

motifsUp

# To prepare this data for plotting with ggplot we can create a simplified data.frame object containing the motif names, the corrected p-values, and the significance rank.
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

# Using ggplot we can plot the rank-sorted TF motifs and color them by the significance of their enrichment. Here we use ggrepel to label each TF motif.
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

# We can perform the same analyses for the peaks that are more accessible in the “C22” cells by using peaks with Log2FC <= -0.5.

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest, #This uses C22
    ArchRProj = projMCS7,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

motifsDo

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

# In this case, the most enriched motifs in the peaks that are more accessible in "C22" cells should correspond to GABA motifs.
head(df)

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

# Save the plots
plotPDF(ggUp, ggDo, name = "C18-vs-C22-Markers-Motifs-Enriched_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

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

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_2024-08-07", width = 8, height = 6, ArchRProj = projMCS7, addDOC = FALSE)


############################################
############################################


## Motif heatmaps by treatment group


# Check treatment assignment
unique(projMCS7$Sample)  # Validate the sample IDs
treatment <- projMCS7$Sample


# Specify which treatment group each sample is in:

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

# Add treatment to ArchR project
projMCS7$treatment <- treatment

# Validate treatment assignment
message("Unique treatments: ", unique(projMCS7$treatment))
table(projMCS7$treatment)  # Check the distribution of cells in each treatment group

# Proceed with splitting and analysis
treatmentGroups <- unique(projMCS7$treatment)
markerFeaturesList <- list()

for (group in treatmentGroups) {
  message("Processing treatment group: ", group)
  cellsForGroup <- projMCS7$cellNames[projMCS7$treatment == group]
  message("Number of cells in group ", group, ": ", length(cellsForGroup))
  
  # Subset ArchR project
  subsetProj <- subsetArchRProject(
    ArchRProj = projMCS7,
    cells = projMCS7$cellNames[projMCS7$treatment == group],
    outputDirectory = paste0("ArchRSubset_", group),
    force = TRUE
  )
  
  # Run getMarkerFeatures
  markerFeatures <- getMarkerFeatures(
    ArchRProj = subsetProj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",  # Use "Clusters" or the relevant column
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  # Store results
  markerFeaturesList[[group]] <- markerFeatures
}

# Create heatmaps
for (group in names(markerFeaturesList)) {
  heatmap <- plotMarkerHeatmap(
    seMarker = markerFeaturesList[[group]],
    cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
    transpose = TRUE
  )
  
  # Save heatmap
  pdf(paste0("Heatmap_", group, ".pdf"))
  print(heatmap)
  dev.off()
}


####Motif heatmaps all look the same after subsetting by treatment group


# Example for a specific treatment group
group <- "t2"  # Change as needed
subsetProj <- subsetArchRProject(
  ArchRProj = projMCS7,
  cells = projMCS7$cellNames[projMCS7$treatment == group],
  force = TRUE
)
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = subsetProj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")




############################################
############################################
############################################
############################################

# ENCODE TFBSs

projMCS7 <- addArchRAnnotations(ArchRProj = projMCS7, collection = "EncodeTFBS")
# Error in addArchRAnnotations(ArchRProj = projMCS7, collection = "EncodeTFBS") : 
#  ArchR mm10 annotations not yet supported! Try LOLA for now!

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LOLA")

############################################
############################################

# Bulk ATAC-seq

projMCS7 <- addArchRAnnotations(ArchRProj = projMCS7, collection = "ATAC")
# Error in addArchRAnnotations(ArchRProj = projMCS7, collection = "ATAC") : 
  ArchR mm10 annotations not yet supported! Try LOLA for now!

############################################
############################################

# Codex TFBS

projMCS7 <- addArchRAnnotations(ArchRProj = projMCS7, collection = "Codex")
# Error in addArchRAnnotations(ArchRProj = projMCS7, collection = "Codex") : 
  ArchR mm10 annotations not yet supported! Try LOLA for now!

############################################
############################################

# Custom Enrichment

# Create custom annotations based on ENCODE ChIP-seq experiments

EncodePeaks <- c(
  Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz",
  Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz",
  Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz",
  Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"
)

# We then add a custom annotation to our ArchRProject using the addPeakAnnotation() function. Here, we call our custom annotation “ChIP”

projMCS7 <- addPeakAnnotations(ArchRProj = projMCS7, regions = EncodePeaks, name = "ChIP")

# Use the custom annotation to perform the peak annotation enrichment

enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projMCS7,
    peakAnnotation = "ChIP",
    cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
  )

enrichRegions

heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
# Error in plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE) : 
#  No enrichments found for your cutoff!

#ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

############################################
############################################

## Motif Deviations

# First check to make sure we've added motif annotations to the ArchRProject
if("Motif" %ni% names(projMCS7@peakAnnotation)){c
    projMCS7 <- addMotifAnnotations(ArchRProj = projMCS7, motifSet = "cisbp", name = "Motif")
}

# We also need to add a set of background peaks which are used in computing deviations. 
# Background peaks are chosen using the chromVAR::getBackgroundPeaks() function which samples peaks based on similarity in GC-content and number of fragments across all samples using the Mahalanobis distance.
projMCS7 <- addBgdPeaks(projMCS7)

# We are now ready to compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. This function has an optional parameter called matrixName that allows us to define the name of the deviations matrix that will be stored in the Arrow files. 
# If we do not provide a value to this parameter, as in the example below, this function creates a matrix name by adding the word “Matrix” to the name of the peakAnnotation. 
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
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

##################################

# Extract a subset of motifs for downstream analysis using the getFeatures() function. 
# The paste(motifs, collapse="|") statement below creates a concatenated or statement that enables selection of all of the motifs.
motifs <- c("GATA1", "Slc17a1", "Slc17a6", "Slc17a8", "Slc32a1")
markerMotifs <- getFeatures(projMCS7, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

# To get just the features corresponding to z-scores, we can use grep.
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

# In the example motifs, in addition to “EBF1”, we also selected “SREBF1” which we do not want to analyze. 
# We can remove it using the %ni% expression which is an ArchR helper function that provides the opposite of %in% from base R.
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

# Now that we have the names of the features that we are interested in, we can plot the distribution of chromVAR deviation scores for each cluster. 
p <- plotGroups(ArchRProj = projMCS7, 
  groupBy = "Clusters", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projMCS7)
)

# Use cowplot to plot the distributions of all these motifs in a single plot

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

# Save the plots
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation_2024-08-07", width = 5, height = 5, ArchRProj = projMCS7, addDOC = FALSE)

##################

# Instead of looking at the distributions of these z-scores, we can overlay the z-scores on a UMAP embedding.
p <- plotEmbedding(
    ArchRProj = projMCS7, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme5)
)

# Plot all of the UMAPs using cowplot
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

####################

# To see how these TF deviation z-scores compare to the inferred gene expression via gene scores of the corresponding TF genes, we can overlay the gene scores for each of these TFs on the UMAP embedding.
markerRNA <- getFeatures(projMCS7, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("Slc17a7","Slc17a6", "Slc17a8", "Slc32a1")]
markerRNA

p <- plotEmbedding(
    ArchRProj = projMCS7, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projMCS7)
)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

