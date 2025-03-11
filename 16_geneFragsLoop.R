library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(GenomicRanges)

#Additional setup
setwd("/Volumes/DataBox/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS1 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS1", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS1)

############################################
############################################

# Step 1: Identify all samples in your project
cellData <- getCellColData(ArchRProj = projMCS1)
sampleNames <- unique(cellData$Sample)

## Store the extracted matrices for each sample

# Prepare a list to store peak matrices for each sample
geneMatrix_list <- list()

# Loop through each sample to extract peak matrices
for (sample in sampleNames) {
  
  # Subset the ArchR project by the current sample
  sample_cells <- rownames(cellData[cellData$Sample == sample, ])
  subsetProj <- subsetArchRProject(
    ArchRProj = projMCS1,
    cells = sample_cells,
    dropCells = TRUE
  )
  
  # Extract the peak matrix for the current sample
  geneMatrix <- getMatrixFromProject(
    ArchRProj = subsetProj,
    useMatrix = "GeneScoreMatrix",
    verbose = FALSE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile(paste0("getMatrixFromProject_", sample))
  )
  
  # Store the peak matrix for later use
  geneMatrix_list[[sample]] <- peakMatrix
}

# After the loop, geneMatrix_list contains peak matrices for all samples

######################################
######################################

## Subset the stored matrices by additional regions of interest

# Define a new region of interest
new_region <- GRanges(
  seqnames = "chr16",
  ranges = IRanges(start = 86000000, end = 86200000)
)

# Prepare a list to store results for the new region
new_region_results <- list()

# Loop through the stored peak matrices and analyze the new region
for (sample in names(peakMatrix_list)) {
  
  peakMatrix <- peakMatrix_list[[sample]]
  
  # Get the peak ranges
  peakRanges <- rowRanges(peakMatrix)
  
  # Subset the peak ranges by the new region of interest
  subset_peaks <- subsetByOverlaps(peakRanges, new_region)
  
  # Extract the indices from the subset_peaks
  indices <- subset_peaks$idx
  
  # Subset the PeakMatrix using the extracted indices
  subset_fragments <- peakMatrix[indices, , drop = FALSE]
  
  # Sum the fragment counts across the subset region for each sample
  total_fragments <- colSums(assay(subset_fragments))
  
  # Store the total fragments result for the current sample
  new_region_results[[sample]] <- total_fragments
}

# View the results for the new region of interest
new_region_results



  
