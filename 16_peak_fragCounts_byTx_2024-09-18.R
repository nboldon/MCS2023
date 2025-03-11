#Setup an interactive session
salloc --account=eon -t 0-16:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Updated conda env 12-2023
module load miniconda3/23.1.0

conda activate archr2023_12

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
############################################
############################################


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
tx_names <- unique(projMCS7$treatment)

# Loop over each genomic region
for (region_name in names(regions)) {
  
  # Initialize an empty data frame to store the results for the current region
  results <- data.frame(Sample = character(), stringsAsFactors = FALSE)
  
  # Get the current genomic region from the list
  region <- regions[[region_name]]
  
  # Loop over each treatment group (assumed to be in tx_names)
  for (i in tx_names) {
    
    # Get the indices of the cells corresponding to the current treatment
    idxTreatment <- BiocGenerics::which(projMCS7$treatment %in% i)
    cellsTreatment <- projMCS7$cellNames[idxTreatment]
    
    # Subset the project for the current treatment group
    projTreatment <- projMCS7[cellsTreatment, ]
    
    # Extract the peak matrix for the current treatment
    peakMatrix <- getMatrixFromProject(
      ArchRProj = projTreatment,
      useMatrix = "PeakMatrix",
      verbose = FALSE,
      binarize = FALSE,
      threads = getArchRThreads(),
      logFile = createLogFile(paste0("getMatrixFromProject_", i))
    )
    
    # Get the row ranges (genomic coordinates) of the peaks
    peakRanges <- rowRanges(peakMatrix)
    
    # Subset the peak ranges by the genomic region of interest
    subset_peaks <- subsetByOverlaps(peakRanges, region)
    
    # Get the indices of the overlapping peaks
    indices <- subset_peaks$idx
    
    # Subset the peak matrix using these indices
    subset_fragments <- peakMatrix[indices, , drop = FALSE]
    
    # Calculate the total fragments for each cell
    total_fragments <- colSums(assay(subset_fragments))
                                     
    # Sum the total fragments across all cells for the current treatment
    sum_fragments <- sum(total_fragments)
                                     
    # Add the results to the data frame
    results <- rbind(results, data.frame(Sample = i, sum_fragments))
  }
  
  # Rename the column 'sum_fragments' to reflect the region name (e.g., "TotalFragments")
  colnames(results)[2] <- paste0(region_name, "_peakFrags")
  
  # Write the results to a CSV file named after the region and treatment
  write.csv(results, paste0(region_name, "_peakFrag_counts_byTx_2024-09-18.csv"), row.names = FALSE)
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
