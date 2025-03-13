

# To preserve both sample assignments and treatment group information while maintaining the integrity of your SummarizedExperiment for downstream analyses, 
# you can add the treatment information as a new column in the colData of your ArchR project. 


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
setwd("/Volumes/DataBox/MCS2023/Stats/Lasso_CellTypes/")
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



# Assign treatments based on sample IDs and add a column to colData
# t1 = 2N, t2 = 2N+, t3 = Ts, t4 = Ts+

# Create a named vector to map sample IDs to treatments
treatment_mapping <- c(
  "C302_" = "t1", "C306_" = "t1", "C309_" = "t1", "C318_" = "t1", 
  "C323_" = "t1", "C328_" = "t1", "C332_" = "t1", "C337_" = "t1", 
  "C339_" = "t1", "C346_" = "t1", "C351_" = "t1", "C353_" = "t1", 
  "C360_" = "t1", "C304_" = "t2", "C308_" = "t2", "C312_" = "t2", 
  "C349_" = "t2", "C315_" = "t2", "C321_" = "t2", "C324_" = "t2", 
  "C355_" = "t2", "C327_" = "t2", "C330_" = "t2", "C333_" = "t2", 
  "C358_" = "t2", "C336_" = "t2", "C342_" = "t2", "C348_" = "t2", 
  "C362_" = "t2", "C305_" = "t3", "C307_" = "t3", "C313_" = "t3", 
  "C350_" = "t3", "C316_" = "t3", "C320_" = "t3", "C322_" = "t3", 
  "C352_" = "t3", "C325_" = "t3", "C334_" = "t3", "C359_" = "t3", 
  "C340_" = "t3", "C341_" = "t3", "C345_" = "t3", "C364_" = "t3", 
  "C301_" = "t4", "C303_" = "t4", "C310_" = "t4", "C314_" = "t4", 
  "C319_" = "t4", "C335_" = "t4", "C338_" = "t4", "C344_" = "t4", 
  "C354_" = "t4", "C356_" = "t4", "C361_" = "t4", "C363_" = "t4"
)

# Use the sample names in the ArchR project to assign treatment groups
microglia_ArchRSubset$treatment <- treatment_mapping[microglia_ArchRSubset$Sample]

# Verify treatment assignments
table(microglia_ArchRSubset$treatment)

# Ensure the column is added to colData (so it is preserved for downstream analysis)
microglia_ArchRSubset <- addCellColData(
  ArchRProj = microglia_ArchRSubset,
  data = microglia_ArchRSubset$treatment,
  name = "TreatmentGroup",
  cells = getCellNames(microglia_ArchRSubset) # Ensure cell barcodes are passed correctly
)

# Check unique samples in the ArchR project
unique(projMCS7$Sample)

# Check how the treatment mapping aligns
table(treatment_mapping[names(treatment_mapping) %in% unique(projMCS7$Sample)])

# Confirm treatment & sample assignments
table(microglia_ArchRSubset$treatment)
table(microglia_ArchRSubset$Sample)


endoMarkers_bySample <- getMarkerFeatures(
  ArchRProj = endo_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)



#########################
#########################
#########################
#########################
#########################



## MICROGLIA ERROR: Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 29 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C334_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C334_>
  Error Found Iteration 51 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C359_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C359_>
  Error Found Iteration 53 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C361_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C361_>
  In addition: Warning message:
  In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  3 function calls resulted in an error
Error Found Iteration 48 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C356_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C356_>
  Error Found Iteration 46 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C354_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C354_>
  

  
## ENDO ERROR:
  
  Error Found Iteration 7 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C307_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C307_>
  Error Found Iteration 22 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C324_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C324_>
  Error Found Iteration 27 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C332_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C332_>
  Error Found Iteration 30 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C335_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C335_>
  Error Found Iteration 32 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C341_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C341_>
  Error Found Iteration 33 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C342_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C342_>
  Error Found Iteration 51 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C363_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C363_>
  Error Found Iteration 52 : 
  [1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C364_\n"
<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C364_>
  
  
  
## Microglia markerFeatures resulted in error for C334, C354, C356, C359, and C361
# Exclude the samples by filtering based on the sample names
# MICROGLIA #excluded_samples <- c("C334_", "C354_", "C356_", "C359_", "C361_")
  
# ENDO-VASC  
excluded_samples <- c("C307_", "C324_", "C332_", "C335_", "C341_", "C342_", "C363_", "C364")
  
# Subset the ArchR project to remove these samples
microglia_ArchRSubset <- microglia_ArchRSubset[!microglia_ArchRSubset$Sample %in% excluded_samples, ]

# Check the updated sample distribution
table(microglia_ArchRSubset$Sample)

# Check the first few rows of the project
head(getCellColData(endo_ArchRSubset))


microgliaMarkers_bySample <- getMarkerFeatures(
  ArchRProj = microglia_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)



#########################
#########################
#########################
#########################
#########################


#Get the list of markers
microglia_markerList_bySample <- getMarkers(microgliaMarkers_bySample, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")


# 1/27/25 cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
# 1/28/25 cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1.15"
#   same as pairwise comp: cutOff = FDR <= 0.1 & abs(Log2FC) >=0.5

##########################


str(microglia_markerList_bySample[[1]])


# Define treatment mapping
treatment_mapping <- c(
  "C302_" = "t1", "C306_" = "t1", "C309_" = "t1", "C318_" = "t1", 
  "C323_" = "t1", "C328_" = "t1", "C332_" = "t1", "C337_" = "t1", 
  "C339_" = "t1", "C346_" = "t1", "C351_" = "t1", "C353_" = "t1", 
  "C360_" = "t1", "C304_" = "t2", "C308_" = "t2", "C312_" = "t2", 
  "C349_" = "t2", "C315_" = "t2", "C321_" = "t2", "C324_" = "t2", 
  "C355_" = "t2", "C327_" = "t2", "C330_" = "t2", "C333_" = "t2", 
  "C358_" = "t2", "C336_" = "t2", "C342_" = "t2", "C348_" = "t2", 
  "C362_" = "t2", "C305_" = "t3", "C307_" = "t3", "C313_" = "t3", 
  "C350_" = "t3", "C316_" = "t3", "C320_" = "t3", "C322_" = "t3", 
  "C352_" = "t3", "C325_" = "t3", "C334_" = "t3", "C359_" = "t3", 
  "C340_" = "t3", "C341_" = "t3", "C345_" = "t3", "C364_" = "t3", 
  "C301_" = "t4", "C303_" = "t4", "C310_" = "t4", "C314_" = "t4", 
  "C319_" = "t4", "C335_" = "t4", "C338_" = "t4", "C344_" = "t4", 
  "C354_" = "t4", "C356_" = "t4", "C361_" = "t4", "C363_" = "t4"
)

# Initialize an empty list
all_markers <- list()

# Loop through each sample's markers
for (sample in names(microglia_markerList_bySample)) {
  # Extract markers for the current sample
  sample_markers <- microglia_markerList_bySample[[sample]]
  
  # Check if the DFrame has any rows
  if (nrow(sample_markers) > 0) {
    # Convert to a data frame (if needed)
    sample_markers <- as.data.frame(sample_markers)
    
    # Add sample name as a new column
    sample_markers$Sample <- sample
    
    # Map the sample to its treatment group and add as a new column
    sample_markers$TreatmentGroup <- treatment_mapping[sample]
    
    # Store the result in the list
    all_markers[[sample]] <- sample_markers
  }
}

# Combine all data frames in the list into one data frame
if (length(all_markers) > 0) {
  combined_markers <- do.call(rbind, all_markers)
  
  # Write the combined data frame to a CSV file
  write.csv(combined_markers, file = "microglia_markers_bySample_2025-01-28.csv", row.names = FALSE)
  message("File successfully saved.")
} else {
  message("No markers passed the filter criteria.")
}

## ENDO-VASC RETURN:
# No markers passed the filter criteria.



##################################


# Print list of Samples with returns
if (exists("combined_markers")) {
  # Extract the unique sample names
  sample_list <- unique(combined_markers$Sample)
  
  # Print the sample list
  print(sample_list)
  
  # Optionally, save the sample list to a CSV file
#  write.csv(sample_list, file = "sample_list.csv", row.names = FALSE, quote = FALSE)
  
  message("Sample list saved to 'sample_list.csv'.")
}

# ASTRO RETURN: 
# t1: C339_, C351_
# t2: C304_, C308_
# t3: C305_, C306_, C307_
# t4: C301_, C303_, C314_, C354_

# GLUT RETURN: 
# t1: C302_, C306_, C309_, C318_, C323_, C328_, C332_, C337_, C339_, C346_, C351_, C353_, C360_
# t2: C304_, C308_, C312_, C315_, C321_, C324_, C327_, C330_, C333_, C336_, C342_, C348_, C349_, C355_, C358_, C362_
# t3: C305_, C307_, C313_, C316_, C320_, C322_, C325_, C334_, C340_, C341_, C345_, C350_, C352_, C359_, C364_
# t4: C301_, C303_, C310_, C314_, C319_, C335_, C338_, C344_, C354_, C356_, C361_, C363_

# GABA RETURN: 
# t1: C306_, C309_, C351_, C360_
# t2: C304_, C308_, C321_, C324_, C342_, C358_
# t3: C307_, C307_, C340_
# t4: C303_, C319_, C338_

# MICROGLIA RETURN: 
# t1: C302_, C306_, C309_, C318_, C323_, C328_, C332_, C337_, C339_, C346_, C351_, C353_, C360_
# t2: C304_, C308_, C312_, C315_, C321_, C324_, C327_, C330_, C333_, C336_, C342_, C348_, C349_, C355_, C358_, C362_
# t3: C305_, C307_, C313_, C316_, C320_, C322_, C325_, C334_, C340_, C341_, C345_, C350_, C352_, C359_, C364_
# t4: C301_, C303_, C310_, C314_, C319_, C335_, C338_, C344_, C354_, C356_, C361_, C363_

# OLIGO RETURN: 
# t1: C302_, C306_, C309_, C318_, C323_, C328_, C332_, C337_, C339_, C346_, C351_, C353_, C360_
# t2: C304_, C308_, C312_, C315_, C321_, C324_, C327_, C330_, C333_, C336_, C342_, C348_, C349_, C355_, C358_, C362_
# t3: C305_, C307_, C313_, C316_, C320_, C322_, C325_, C334_, C340_, C341_, C345_, C350_, C352_, C359_, C364_
# t4: C301_, C303_, C310_, C314_, C319_, C335_, C338_, C344_, C354_, C356_, C361_, C363_

# ENDO-VASC RETURN: 
# t1: C302_, C306_, C309_, C318_, C323_, C328_, C332_, C337_, C339_, C346_, C351_, C353_, C360_
# t2: C304_, C308_, C312_, C315_, C321_, C324_, C327_, C330_, C333_, C336_, C342_, C348_, C349_, C355_, C358_, C362_
# t3: C305_, C307_, C313_, C316_, C320_, C322_, C325_, C334_, C340_, C341_, C345_, C350_, C352_, C359_, C364_
# t4: C301_, C303_, C310_, C314_, C319_, C335_, C338_, C344_, C354_, C356_, C361_, C363_



##################################


## Print total number of returns by sample and treatment group
# Count the total entries for each sample
sample_counts <- table(projMCS7$Sample)

# Map each sample to its treatment group
sample_treatment <- data.frame(
  Sample = names(sample_counts),
  Count = as.vector(sample_counts),
  Treatment = treatment_mapping[names(sample_counts)]
)

# Summarize total counts for each treatment group
treatment_summary <- aggregate(Count ~ Treatment, data = sample_treatment, sum)

# Print the results
print(sample_treatment)
print(treatment_summary)


## Add sig gene count by sample and treatment and print results
# Step 1: Calculate the number of genes for each sample from glut_markerList_bySample
gene_counts <- sapply(names(endo_markerList_bySample), function(sample) {
  # Check if there are markers for this sample
  if (!is.null(endo_markerList_bySample[[sample]])) {
    nrow(endo_markerList_bySample[[sample]]) # Count number of rows (genes) for this sample
  } else {
    0 # If no markers, return 0
  }
})

# Step 2: Create a data frame for gene counts
gene_counts_df <- data.frame(
  Sample = names(gene_counts),
  GeneCount = as.vector(gene_counts) # Convert to a numeric vector
)

# Step 3: Merge the gene counts with sample_treatment
sample_treatment <- merge(
  sample_treatment, 
  gene_counts_df, 
  by = "Sample", 
  all.x = TRUE
)

# Step 4: Replace NA values in GeneCount with 0 (in case some samples have no markers)
sample_treatment$GeneCount[is.na(sample_treatment$GeneCount)] <- 0

# Step 5: Summarize total gene counts for each treatment group
treatment_summary <- aggregate(cbind(Count, GeneCount) ~ Treatment, data = sample_treatment, sum)

# Step 6: Print or save the results
print(sample_treatment)
print(treatment_summary)

# Save the updated sample_treatment to a CSV
write.csv(sample_treatment, file = "endo_bySampleTx_withGene-Cell-Counts_2025-01-28.csv", row.names = FALSE)
message("Sample treatment data with gene counts saved.")

write.csv(treatment_summary, file = "endo_byTx_withGene-Cell-Counts_2025-01-28.csv", row.names = FALSE)
message("Sample treatment data with gene counts saved.")




# Print number of gene markers in each treatment group with returns
if (exists("combined_markers")) {
  # Count the number of entries for each treatment group
  treatment_counts <- table(combined_markers$TreatmentGroup)
  
  # Print the counts
  print(treatment_counts)
}


## ASTRO RETURN:
# t1  t2  t3  t4 
# 139 255  97 100

## GLUT RETURN:
# t1   t2   t3   t4 
# 984 1169 1223  893 

## GABA RETURN:
# t1 t2 t3 t4 
# 8  8  5  8 

## MICROGLIA RETURN:
#   t1   t2   t3   t4 
# 109  922  119 1060 

## OLIGO RETURN:
#   t1   t2   t3   t4 
# 878 1902 1927  888 

## ENDO-VASC RETURN:
# NULL



# Identify samples with no markers
empty_samples <- names(endo_markerList_bySample)[sapply(endo_markerList_bySample, nrow) == 0]
print(empty_samples)

## ASTRO RETURN:
[1] "C302_" "C309_" "C310_" "C312_" "C313_" "C315_" "C316_" "C318_" "C319_" "C320_" "C321_" "C322_" "C323_" "C324_"
[15] "C325_" "C327_" "C328_" "C330_" "C332_" "C333_" "C334_" "C335_" "C336_" "C337_" "C338_" "C340_" "C341_" "C342_"
[29] "C344_" "C345_" "C346_" "C348_" "C349_" "C350_" "C352_" "C353_" "C355_" "C356_" "C358_" "C359_" "C360_" "C361_"
[43] "C362_" "C363_" "C364_"


## ENDO RETURN (49):
[1] "C301_" "C302_" "C303_" "C304_" "C305_" "C306_" "C308_" "C309_" "C310_" "C312_" "C313_" "C314_" "C315_"
[14] "C316_" "C318_" "C319_" "C320_" "C321_" "C322_" "C323_" "C325_" "C327_" "C328_" "C330_" "C333_" "C334_"
[27] "C336_" "C337_" "C338_" "C339_" "C340_" "C344_" "C345_" "C346_" "C348_" "C349_" "C350_" "C351_" "C352_"
[40] "C353_" "C354_" "C355_" "C356_" "C358_" "C359_" "C360_" "C361_" "C362_" "C364_"




#############################################
######################################################################


## ASTRO Returned results for C304, C305, C306
## 1-27-25 cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"


> empty_samples <- names(astro_markerList_bySample)[sapply(astro_markerList_bySample, nrow) == 0]
> print(empty_samples)
[1] "C301_" "C302_" "C303_" "C307_" "C308_" "C309_" "C310_" "C312_" "C313_" "C314_" "C315_" "C316_"
[13] "C318_" "C319_" "C320_" "C321_" "C322_" "C323_" "C324_" "C325_" "C327_" "C328_" "C330_" "C332_"
[25] "C333_" "C334_" "C335_" "C336_" "C337_" "C338_" "C339_" "C340_" "C341_" "C342_" "C344_" "C345_"
[37] "C346_" "C348_" "C349_" "C350_" "C351_" "C352_" "C353_" "C354_" "C355_" "C356_" "C358_" "C359_"
[49] "C360_" "C361_" "C362_" "C363_" "C364_"

> all_markers
$C304_
seqnames     start       end strand          name  idx   Log2FC         FDR  MeanDiff Sample
6652     chr14  57890262  57832702      2       Zdhhc20  538 2.320057 0.002898615 0.4179886  C304_
21648     chr8 105827350 105768308      2       Ranbp10  883 1.758948 0.002898615 0.5725293  C304_
21649     chr8 105827744 105844676      1      Tsnaxip1  884 2.009558 0.002898615 0.6568741  C304_
9475     chr17  43667015  43641900      2      Slc25a27  803 2.541578 0.009599654 0.3003103  C304_
12783     chr2 146542931 146546108      1 4933406D12Rik 1530 1.325559 0.009599654 0.4879197  C304_
17923     chr6  85587531  85702751      1         Alms1  636 1.423373 0.009599654 0.1998046  C304_

$C305_
seqnames    start      end strand  name idx    Log2FC         FDR  MeanDiff Sample
18869     chr7 17062427 17030993      2 Hif3a 273 -1.328039 0.007222315 -1.589377  C305_

$C306_
seqnames     start       end strand   name  idx    Log2FC         FDR   MeanDiff Sample
14537     chr4  24851086  24617273      2 Klhl32   98 -1.737125 0.001760317 -0.6694671  C306_
14435     chr3 158562221 158082895      2  Lrrc7 1137 -2.162724 0.008341922 -0.5505981  C306_
19224     chr7  30090510  30056034      2  Zfp82  628  1.557175 0.008341922  0.6246821  C306_
19225     chr7  30090510  30077337      2 Zfp566  629  1.594460 0.009923369  0.7592240  C306_

>
  
  
  
  ##################################
##################################
##################################



## ASTRO Returned results for C301, C303, C304, C305, C306, C307, C314
## 1-28-25 cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1.15"



> print(empty_samples)
[1] "C302_" "C309_" "C310_" "C312_" "C313_" "C315_" "C316_" "C318_" "C319_" "C320_" "C321_" "C322_"
[13] "C323_" "C324_" "C325_" "C327_" "C328_" "C330_" "C332_" "C333_" "C334_" "C335_" "C336_" "C337_"
[25] "C338_" "C339_" "C340_" "C341_" "C342_" "C344_" "C345_" "C346_" "C348_" "C349_" "C350_" "C351_"
[37] "C352_" "C353_" "C354_" "C355_" "C356_" "C358_" "C359_" "C360_" "C361_" "C362_" "C363_" "C364_"
> all_markers
$C301_
seqnames    start      end strand  name idx   Log2FC       FDR  MeanDiff Sample
8593     chr16 91011308 90936097      2 Synj1 690 1.696372 0.0143197 0.4727632  C301_
10377    chr18 74216212 74221491      1 Cxxc1 521 2.261943 0.0143197 1.7783624  C301_

$C303_
seqnames     start       end strand          name  idx    Log2FC        FDR   MeanDiff Sample
7613     chr15  88314874  88372651      1 B230214G05Rik  623  3.102189 0.01249812  0.2890414  C303_
9729     chr17  71526857  71496100      2         Ndc80 1057  2.903894 0.01249812  0.5361651  C303_
23676     chrX  61709616  61710950      1         Ldoc1  386 -2.855758 0.01249812 -0.9113838  C303_
18550     chr6 145210970 145174834      2         Casc1 1263  2.677164 0.01754320  0.6039052  C303_
11765     chr2  42653598  40596773      2         Lrp1b  512 -1.185300 0.03318108 -1.9839440  C303_

$C304_
seqnames     start       end strand          name  idx   Log2FC         FDR  MeanDiff Sample
6652     chr14  57890262  57832702      2       Zdhhc20  538 2.320057 0.002898615 0.4179886  C304_
21648     chr8 105827350 105768308      2       Ranbp10  883 1.758948 0.002898615 0.5725293  C304_
21649     chr8 105827744 105844676      1      Tsnaxip1  884 2.009558 0.002898615 0.6568741  C304_
9475     chr17  43667015  43641900      2      Slc25a27  803 2.541578 0.009599654 0.3003103  C304_
12783     chr2 146542931 146546108      1 4933406D12Rik 1530 1.325559 0.009599654 0.4879197  C304_
17923     chr6  85587531  85702751      1         Alms1  636 1.423373 0.009599654 0.1998046  C304_
6178     chr14  18271136  18284003      1       Nkiras1   64 1.631326 0.013008612 0.4681380  C304_
6739     chr14  65833969  65819027      2         Esco2  625 2.070564 0.013008612 0.5672992  C304_
4770     chr12  83597147  83546941      2        Zfyve1  427 1.511045 0.016207940 0.4052070  C304_
7950     chr16   5255956   5244155      2       Eef2kmt   47 1.203761 0.016207940 0.4994251  C304_
6179     chr14  18270986  18267823      2         Rpl15   65 2.047163 0.026597284 0.6350114  C304_
12839     chr2 151741310 151716062      2         Psmf1 1586 2.761011 0.026597284 0.2851497  C304_
13166     chr2 172345577 172355749      1       Fam210b 1913 1.212260 0.037016988 0.5323112  C304_
22650     chr9  64022059  63953076      2         Smad6  712 1.188664 0.037016988 0.8586789  C304_
20387     chr7 123500449 123478631      2       Zkscan2 1791 1.689417 0.049081760 0.1199082  C304_

$C305_
seqnames     start       end strand   name idx    Log2FC         FDR   MeanDiff Sample
18869     chr7  17062427  17030993      2  Hif3a 273 -1.328039 0.007222315 -1.5893767  C305_
13466     chr3  38886940  39011985      1   Fat4 168 -1.433273 0.011304926 -0.7078240  C305_
9120     chr17  30612659  30592861      2   Glo1 448  1.163437 0.015733073  1.0619410  C305_
15071     chr4 107253930 107279888      1 Hspb11 632  1.972234 0.044547919  0.6866349  C305_

$C306_
seqnames     start       end strand          name  idx    Log2FC         FDR   MeanDiff Sample
14537     chr4  24851086  24617273      2        Klhl32   98 -1.737125 0.001760317 -0.6694671  C306_
14435     chr3 158562221 158082895      2         Lrrc7 1137 -2.162724 0.008341922 -0.5505981  C306_
19224     chr7  30090510  30056034      2         Zfp82  628  1.557175 0.008341922  0.6246821  C306_
19225     chr7  30090510  30077337      2        Zfp566  629  1.594460 0.009923369  0.7592240  C306_
1722     chr10  62252438  62265990      1         Tacr2  364  1.558557 0.016951801  0.6086250  C306_
18227     chr6 121636173 121679238      1           A2m  940 -2.944084 0.020480799 -0.2078223  C306_
12891     chr2 153345810 153341157      2 2500004C02Rik 1638  1.799348 0.024832438  0.7265233  C306_
292       chr1  58445486  58430993      2         Ppil3  292  2.246174 0.025002919  0.3288973  C306_
8478     chr16  62786716  62734852      2         Nsun3  575 -2.435043 0.025002919 -0.3460746  C306_
9758     chr17  75861023  75889039      1        Gm4710 1086 -3.189276 0.025002919 -0.2379533  C306_
10349    chr18  67800107  67885170      1        Cep192  493  2.394168 0.025002919  0.2480842  C306_
12783     chr2 146542931 146546108      1 4933406D12Rik 1530  1.185961 0.025002919  0.5585985  C306_
13331     chr3  13471655  14182287      1         Ralyl   33 -1.518110 0.025002919 -0.5174862  C306_
17594     chr6  45060061  47301388      1       Cntnap2  307 -1.235152 0.025002919 -1.4103376  C306_
3573     chr11  85235166  85238304      1      Appbp2os 1068  2.636675 0.026063644  0.5644297  C306_
8516     chr16  77646273  77646343      1     Mir125b-2  613  1.410863 0.030299927  1.3072157  C306_
12890     chr2 153346139 153404007      1         Asxl1 1637  1.554815 0.030533250  0.3459531  C306_
19226     chr7  30095076  30107614      1        Zfp260  630  1.364106 0.032103960  0.6868223  C306_
6729     chr14  64950045  65039835      1         Ints9  615  1.752288 0.032543848  0.6865237  C306_
21063     chr8  35589020  35582503      2       Gm16793  298  1.194360 0.032543848  0.7054949  C306_
14172     chr3 114904078 115125764      1         Olfm3  874 -1.789329 0.035047183 -0.4012430  C306_
5258     chr13  20090507  20608353      1         Elmo1   83 -1.419298 0.038289676 -0.6058130  C306_
416       chr1  74588036  74566426      2        Zfp142  416  1.162198 0.038295628  0.6883641  C306_
20512     chr7 128461513 128439777      2         Tial1 1916  1.170578 0.045032049  0.5258602  C306_
21778     chr8 116978959 116969599      2 1700030J22Rik 1013  2.008450 0.045517677  0.4857662  C306_
21674     chr8 106415339 106425895      1         Zfp90  909  1.314224 0.048945385  0.2932243  C306_
21776     chr8 116921740 116941503      1         Cenpn 1011  1.231640 0.048948892  0.3048128  C306_
11133    chr19  46304737  46312385      1         Nfkb2  677  1.717336 0.049453489  0.8578508  C306_

$C307_
seqnames     start       end strand          name  idx    Log2FC        FDR   MeanDiff Sample
18367     chr6 128891124 128848044      2      BC035044 1080 -1.717699 0.01419848 -0.5585953  C307_
10348    chr18  67799943  67776033      2 4930549G23Rik  492  2.050767 0.01638705  0.4772810  C307_
12717     chr2 132816054 132804215      2         Trmt6 1464  3.692899 0.01638705  0.4213816  C307_
12718     chr2 132816141 132844197      1          Mcm8 1465  3.628988 0.01638705  0.3285052  C307_
17841     chr6  77979667  76881637      2        Ctnna2  554 -1.175130 0.01638705 -1.1726072  C307_
18366     chr6 128826315 128813706      2        Klrb1b 1079 -1.256902 0.01638705 -0.4417752  C307_
19866     chr7  89632505  89854359      1           Me3 1270 -1.864900 0.01638705 -0.7198090  C307_
21921     chr8 125995102 126030685      1         Kcnk1 1156  1.209914 0.01638705  1.0002014  C307_
16152     chr5  33898971  33898910      2       Mir7024  250  2.002188 0.02727196  0.5125176  C307_
1523     chr10  28074820  28597397      1         Ptprk  165 -1.498798 0.03103340 -0.6475167  C307_
21713     chr8 109868603 109944671      1        Phlpp2  948  1.958678 0.03568975  0.3302918  C307_
21472     chr8  85432841  85408759      2 4921524J17Rik  707  1.687529 0.04190570  0.3699664  C307_
1221      chr1 179546267 179528056      2         Tfb2m 1221  1.442743 0.04730534  0.3673258  C307_
20759     chr7 144896532 144900953      1         Fgf15 2163  3.475352 0.04730534  0.3857752  C307_

$C308_
seqnames     start       end strand          name  idx    Log2FC        FDR   MeanDiff Sample
4082     chr11 108425266 108860615      1        Cep112 1577  1.693930 0.01201420  0.5907835  C308_
10196    chr18  46741876  46790826      1         Ap3s1  340  1.747448 0.01201420  0.5936993  C308_
11765     chr2  42653598  40596773      2         Lrp1b  512 -1.207934 0.01201420 -1.6050834  C308_
1787     chr10  73099342  74649737      1        Pcdh15  429 -1.438762 0.01560690 -0.8842166  C308_
16153     chr5  33983433  33985009      1        Gm1673  251  1.962896 0.01560690  0.6938442  C308_
21301     chr8  71463657  71456700      2         Abhd8  536  1.270438 0.01560690  1.4051335  C308_
21303     chr8  71469194  71476097      1          Dda1  538  1.257547 0.01560690  1.6447839  C308_
23171     chr9 114640200 114585874      2        Cnot10 1233  2.407862 0.01560690  0.3850242  C308_
21302     chr8  71464926  71465753      1        Mrpl34  537  1.295458 0.03206069  1.4607000  C308_
12608     chr2 127208280 127186355      2 1810024B03Rik 1355  1.797320 0.03669216  0.2735543  C308_
1445     chr10  19591958  19610225      1        Ifngr1   87  1.932474 0.04148162  0.5572413  C308_
4772     chr12  83688203  83735199      1         Psen1  429  2.663753 0.04148162  0.3596870  C308_
4890     chr12  99627974  99626053      2 1700064M15Rik  547  1.203820 0.04148162  0.7651938  C308_
5664     chr13  55693124  55703500      1 B230219D22Rik  489  2.849327 0.04148162  0.5337710  C308_
8523     chr16  81200697  81624290      1         Ncam2  620 -1.269834 0.04148162 -0.6736505  C308_
8590     chr16  90810413  90751527      2          Urb1  687  1.611071 0.04148162  0.3524822  C308_
11397     chr2  19909780  20810535      1          Etl4  144 -2.001582 0.04148162 -0.7932667  C308_
13146     chr2 168601657 168476410      2        Nfatc2 1893  1.821352 0.04148162  0.4671329  C308_
17987     chr6  88842854  88847274      1       Gm15612  700  1.213099 0.04148162  1.4872420  C308_
18097     chr6 110645598 111567230      1          Grm7  810 -2.084004 0.04148162 -1.0327230  C308_
21898     chr8 124433936 124369049      2         Pgbd5 1133  1.424997 0.04148162  0.2696023  C308_
17989     chr6  88898780  88883474      2          Mcm2  702  2.338225 0.04162979  0.4337351  C308_
8649     chr16  95929077  95831123      2 1600002D24Rik  746  1.906483 0.04184924  0.1094939  C308_
12758     chr2 144189290 144173150      2        Gm5535 1505 -1.503361 0.04184924 -0.3377108  C308_
16498     chr5  92278181  92257660      2          Naaa  596  1.262598 0.04184924  0.6312415  C308_
17986     chr6  88841935  88835914      2         Abtb1  699  1.172591 0.04184924  1.4238203  C308_

$C314_
seqnames   start     end strand   name idx   Log2FC        FDR MeanDiff Sample
8725    chr17 8526801 8986648      1 Pde10a  53 1.896657 0.03759833 3.043807  C314_

>
  
  
  
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################



## Lasso for all samples/treatment groups:

## Run for each cell type
#gaba, glut, oligo, astro, microglia
## No markers passed the filter criteria for endo-vasc above.

## Run for each behavioral test
# TaskA-45_AllStats, TaskA-45_WrongStats



setwd("/Volumes/DataBox/MCS2023/Stats/Lasso_CellTypes")

# Load necessary libraries
library(glmnet)
library(readr)
library(dplyr)

# Load your data
microglia_data <- read_csv("microglia_markers_bySample_2025-01-28.csv")
behavior_data <- read_csv("TaskA-45_AllStats.csv")

# Clean column names
colnames(behavior_data) <- make.names(colnames(behavior_data))

# Clean up sample IDs
microglia_data <- microglia_data %>% mutate(Sample = gsub("_$", "", Sample))

# Merge the datasets by 'Sample'
merged_data <- merge(microglia_data, behavior_data, by = "Sample")

# Select relevant features and responses
# Use multiple features for Lasso regression
X <- merged_data[, c("Log2FC", "FDR", "MeanDiff")]  # Include additional columns if needed
y <- merged_data[, c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")]

# Standardize the data (standardize both X and y)
X_scaled <- scale(X)
y_scaled <- scale(y)

# Set up Lasso regression model (use cross-validation to select the best lambda)
lasso_model <- cv.glmnet(x = X_scaled,  # Use the matrix X_scaled as input
                         y = y_scaled, 
                         alpha = 1, # alpha = 1 for Lasso (Ridge would be 0)
                         family = "mgaussian")  # Change to "mgaussian" for multivariate Gaussian

# Get the best lambda (regularization parameter)
best_lambda <- lasso_model$lambda.min
cat("Best lambda:", best_lambda, "\n")


# Get the coefficients for the best model
coefficients <- coef(lasso_model, s = "lambda.min")
cat("Coefficients:\n")
print(coefficients)


###################################


# Plot the Lasso path
plot(lasso_model)


# Get predictions for each behavioral response
predictions <- predict(lasso_model, newx = X_scaled, s = "lambda.min")

# Check the structure of predictions
str(predictions)  # Should return an array with dimensions

# Ensure predictions has the same dimensions as y_scaled
if (dim(predictions)[1] == nrow(y_scaled) && dim(predictions)[2] == ncol(y_scaled)) {
  # Compute correlations for each response variable
  correlations <- sapply(1:ncol(y_scaled), function(i) cor(predictions[, i, 1], y_scaled[, i]))
  
  # Print correlations
  names(correlations) <- colnames(y_scaled)
  print(correlations)
} else {
  cat("Error: Dimension mismatch between predictions and y_scaled.\n")
  print(dim(predictions))
  print(dim(y_scaled))
}


##################################
##################################
##################################
##################################


## Loop to run the above through all behavioral tests and cell types


setwd("/Volumes/DataBox/MCS2023/Stats/Lasso_CellTypes")


# Load necessary libraries
library(glmnet)
library(readr)
library(dplyr)

# Define file lists
behavior_files <- c("TaskA-45_AllStats.csv", "TaskA-45_WrongStats.csv", "TaskB-59_AllStats.csv", "TaskB-59_WrongStats.csv",
                    "TaskC-48_AllStats.csv", "TaskC-48_WrongStats.csv", "TaskD-149_AllStats.csv", "TaskD-149_WrongStats.csv")
gene_files <- c("gaba_markers_bySample_2025-01-28.csv", "microglia_markers_bySample_2025-01-28.csv",
                "glut_markers_bySample_2025-01-28.csv", "oligo_markers_bySample_2025-01-28.csv", 
                "astro_markers_bySample_2025-01-28.csv")

# Output directories
dir.create("Lasso_Coefficients", showWarnings = FALSE)
dir.create("Lasso_Correlations", showWarnings = FALSE)

# Loop through all combinations of behavioral and gene score files
for (gene_file in gene_files) {
  for (behavior_file in behavior_files) {
    
    cat("\nProcessing:", gene_file, "with", behavior_file, "\n")
    
    tryCatch({
      # Load data
      gene_data <- read_csv(gene_file)
      behavior_data <- read_csv(behavior_file)
      
      # Clean column names
      colnames(behavior_data) <- make.names(colnames(behavior_data))
      
      # Clean up sample IDs
      gene_data <- gene_data %>% mutate(Sample = gsub("_$", "", Sample))
      
      # Merge datasets
      merged_data <- merge(gene_data, behavior_data, by = "Sample")
      colnames(merged_data) <- make.names(colnames(merged_data))
      
      # Select features and responses
      if (!all(c("Log2FC", "FDR", "MeanDiff") %in% colnames(merged_data))) {
        stop("Missing necessary predictor columns in merged_data")
      }
      
      if (!all(c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice") %in% colnames(merged_data))) {
        stop("Missing necessary response columns in merged_data")
      }
      
      X <- merged_data[, c("Log2FC", "FDR", "MeanDiff"), drop = FALSE]
      y <- merged_data[, c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice"), drop = FALSE]
      
      # Standardize data
      X_scaled <- scale(X)
      y_scaled <- scale(y)
      
      # Fit Lasso regression
      lasso_model <- cv.glmnet(x = X_scaled, y = y_scaled, alpha = 1, family = "mgaussian")
      best_lambda <- lasso_model$lambda.min
      
      # Save coefficients
      coefficients <- coef(lasso_model, s = "lambda.min")
      
      # Convert coefficients to a data frame
      coef_df <- data.frame(
        Predictor = rownames(as.matrix(coefficients[[1]])),
        Coefficient = as.matrix(coefficients[[1]])[, 1]
      )
      
      # Save coefficients as CSV
      coef_filename <- paste0("Lasso_Coefficients/Coefficients_", sub(".csv", "", gene_file), "_", sub(".csv", "", behavior_file), ".csv")
      write.csv(coef_df, coef_filename, row.names = FALSE)
      
      # Save plot
      plot_filename <- paste0("Lasso_Coefficients/LassoPlot_", sub(".csv", "", gene_file), "_", sub(".csv", "", behavior_file), ".png")
      png(plot_filename)
      plot(lasso_model, main = paste("Lasso Path:", gene_file, "vs", behavior_file))
      dev.off()
      
      # Get predictions
      predictions <- predict(lasso_model, newx = X_scaled, s = "lambda.min")
      
      # Compute and save correlations if dimensions match
      if (dim(predictions)[1] == nrow(y_scaled) && dim(predictions)[2] == ncol(y_scaled)) {
        correlations <- sapply(1:ncol(y_scaled), function(i) cor(predictions[, i, 1], y_scaled[, i]))
        names(correlations) <- colnames(y_scaled)
        
        # Convert to data frame
        correlation_df <- data.frame(
          Behavior = names(correlations),
          Correlation = correlations
        )
        
        # Save correlations as CSV
        corr_filename <- paste0("Lasso_Correlations/Correlations_", sub(".csv", "", gene_file), "_", sub(".csv", "", behavior_file), ".csv")
        write.csv(correlation_df, corr_filename, row.names = FALSE)
        
      } else {
        cat("Error: Dimension mismatch between predictions and y_scaled.\n")
      }
      
    }, error = function(e) {
      cat("Error in processing", gene_file, "with", behavior_file, ":\n", e$message, "\n")
    })
  }
}








##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################



## Lasso to compare differences between treatment groups:
# Include if you want a global model that adjusts for treatment.

### Does NOT run by treatment group 


library(glmnet)
library(readr)
library(dplyr)
library(caret)  # For one-hot encoding

# Load your data
gene_data <- read_csv("gaba_markers_bySample_2025-01-28.csv")
behavior_data <- read_csv("TaskA-45_AllStats.csv")

# Clean column names
colnames(behavior_data) <- make.names(colnames(behavior_data))

# Clean up sample IDs
gene_data <- gene_data %>% mutate(Sample = gsub("_$", "", Sample))

# Merge datasets
merged_data <- merge(gene_data, behavior_data, by = "Sample")

# Ensure TreatmentGroup is a factor
merged_data$TreatmentGroup <- as.factor(merged_data$TreatmentGroup)

# One-hot encode treatment groups
treatment_dummies <- model.matrix(~ TreatmentGroup - 1, data = merged_data)  # Removes intercept

# Select relevant features (add treatment group variables)
X <- cbind(merged_data[, c("Log2FC", "FDR", "MeanDiff")], treatment_dummies)  
y <- merged_data[, c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")]

# Standardize
X_scaled <- scale(X)
y_scaled <- scale(y)

# Fit Lasso regression
lasso_model <- cv.glmnet(x = X_scaled, y = y_scaled, alpha = 1, family = "mgaussian")

# Extract best lambda and coefficients
best_lambda <- lasso_model$lambda.min
coefficients <- coef(lasso_model, s = "lambda.min")

# Print results
cat("Best lambda:", best_lambda, "\n")
# gaba: Best lambda: 0.1145194 

print(coefficients)
# gaba results:
$Correct.Response
8 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)      -2.965302e-15
Log2FC           -1.825300e-01
FDR               9.748654e-02
MeanDiff          4.998982e-02
TreatmentGroupt1  3.354944e-02
TreatmentGroupt2  1.561580e-01
TreatmentGroupt3 -4.788061e-01
TreatmentGroupt4  7.073663e-02

$Premature.Response
8 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)      -4.418122e-16
Log2FC            2.139707e-01
FDR              -1.600932e-01
MeanDiff          3.800334e-02
TreatmentGroupt1  3.043640e-01
TreatmentGroupt2 -1.698173e-01
TreatmentGroupt3 -5.103080e-02
TreatmentGroupt4 -3.132136e-01

$Missed.Response.Window
8 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)      -8.963889e-16
Log2FC            5.700172e-02
FDR               2.699748e-02
MeanDiff         -1.102980e-01
TreatmentGroupt1 -3.194993e-02
TreatmentGroupt2 -2.888073e-01
TreatmentGroupt3  1.153566e-01
TreatmentGroupt4  4.908145e-01

$Wrong.Choice
8 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)      -4.321222e-16
Log2FC            1.062891e-01
FDR              -7.189855e-02
MeanDiff          4.147376e-02
TreatmentGroupt1 -1.194457e-01
TreatmentGroupt2  1.456608e-01
TreatmentGroupt3  4.677798e-01
TreatmentGroupt4 -4.202170e-01




## Run Lasso seperately for each treatment group
# Run if you think gene-behavior relationships differ across treatments.


### Does NOT run by treatment group 




for (treatment in unique(merged_data$TreatmentGroup)) {
  subset_data <- merged_data %>% filter(TreatmentGroup == treatment)
  
  X_subset <- subset_data[, c("Log2FC", "FDR", "MeanDiff")]
  y_subset <- subset_data[, c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")]
  
  X_scaled <- scale(X_subset)
  y_scaled <- scale(y_subset)
  
  lasso_model <- cv.glmnet(x = X_scaled, y = y_scaled, alpha = 1, family = "mgaussian")
  
  cat("\nTreatment:", treatment, "\n")
  print(coef(lasso_model, s = "lambda.min"))
}


## Gaba return:
Treatment: t4 
$Correct.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -4.600944e-16
Log2FC       .           
FDR          .           
MeanDiff     7.664611e-02

$Premature.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)  5.013578e-16
Log2FC       .           
FDR          .           
MeanDiff    -6.331808e-02

$Missed.Response.Window
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -1.644129e-16
Log2FC       .           
FDR          .           
MeanDiff    -7.640030e-02

$Wrong.Choice
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) 1.176530e-17
Log2FC      .           
FDR         .           
MeanDiff    7.611049e-02


Treatment: t2 
$Correct.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) 8.743006e-16
Log2FC      .           
FDR         .           
MeanDiff    .           

$Premature.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) .
Log2FC      .
FDR         .
MeanDiff    .

$Missed.Response.Window
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -1.665335e-16
Log2FC       .           
FDR          .           
MeanDiff     .           

$Wrong.Choice
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -3.191891e-16
Log2FC       .           
FDR          .           
MeanDiff     .           


Treatment: t1 
$Correct.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept)  5.570840e-16
Log2FC      -7.669679e-02
FDR          7.887455e-01
MeanDiff    -1.048946e-01

$Premature.Response
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -9.447020e-17
Log2FC       3.957904e-01
FDR         -4.001443e-02
MeanDiff     2.289685e-01

$Missed.Response.Window
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -7.082624e-17
Log2FC      -2.532513e-02
FDR         -8.295281e-01
MeanDiff     5.686158e-02

$Wrong.Choice
4 x 1 sparse Matrix of class "dgCMatrix"
1
(Intercept) -1.004164e-16
Log2FC      -2.771981e-01
FDR         -7.973440e-01
MeanDiff    -1.100663e-01

Error in matrix(fit$a0[seq(lmu * nc)], nc, lmu, dimnames = list(classnames,  : 
                                                                  length of 'dimnames' [2] not equal to array extent
                                                                In addition: Warning messages:
                                                                  1: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold 
                                                                2: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold 
                                                                3: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold 
                                                                4: from glmnet C++ code (error code -1); Convergence for 1th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned 
                                                                


## Compare prediction accuracy across treatments
# Run if you want to see if the model works better for some treatments.

for (treatment in unique(merged_data$TreatmentGroup)) {
  subset_data <- merged_data %>% filter(TreatmentGroup == treatment)
  
  X_subset <- subset_data[, c("Log2FC", "FDR", "MeanDiff")]
  y_subset <- subset_data[, c("Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")]
  
  X_scaled <- scale(X_subset)
  y_scaled <- scale(y_subset)
  
  lasso_model <- cv.glmnet(x = X_scaled, y = y_scaled, alpha = 1, family = "mgaussian")
  
  predictions <- predict(lasso_model, newx = X_scaled, s = "lambda.min")
  predictions_df <- do.call(cbind, predictions)
  
  correlations <- sapply(1:ncol(y_scaled), function(i) cor(predictions_df[, i], y_scaled[, i]))
  
  cat("\nTreatment:", treatment, "\n")
  print(correlations)
}

                                                                
                                                                
### Does NOT run by treatment group                                                               


cor_matrix_pearson <- cor(merged_data[, c("Log2FC", "FDR", "MeanDiff", 
                                          "Correct.Response", "Premature.Response", 
                                          "Missed.Response.Window", "Wrong.Choice")], 
                          use = "pairwise.complete.obs", method = "pearson")

print(cor_matrix_pearson)
                                                                
                                                                
cor_matrix_spearman <- cor(merged_data[, c("Log2FC", "FDR", "MeanDiff", 
                                           "Correct.Response", "Premature.Response", 
                                           "Missed.Response.Window", "Wrong.Choice")], 
                           use = "pairwise.complete.obs", method = "spearman")

print(cor_matrix_spearman)

                                                                
library(ggcorrplot)

ggcorrplot(cor_matrix_pearson, method = "circle", type = "lower", lab = TRUE)

                                                                
                                                                
###################


## Compute correlations seperately for each treatment group

treatment_groups <- unique(merged_data$TreatmentGroup)

cor_list <- list()

for (treatment in treatment_groups) {
  subset_data <- merged_data %>% filter(TreatmentGroup == treatment)
  
  cor_matrix <- cor(subset_data[, c("Log2FC", "FDR", "MeanDiff", 
                                    "Correct.Response", "Premature.Response", 
                                    "Missed.Response.Window", "Wrong.Choice")], 
                    use = "pairwise.complete.obs", method = "pearson")
  
  cor_list[[treatment]] <- cor_matrix
  
  cat("\nCorrelation Matrix for Treatment:", treatment, "\n")
  print(cor_matrix)
}


# Gaba return:
Correlation Matrix for Treatment: t4 
Log2FC         FDR   MeanDiff Correct.Response Premature.Response
Log2FC                  1.00000000 -0.06292701 -0.3472047       -0.1420054          0.2317134
FDR                    -0.06292701  1.00000000 -0.5706564       -0.6023859          0.1075453
MeanDiff               -0.34720472 -0.57065642  1.0000000        0.8627700         -0.7127425
Correct.Response       -0.14200541 -0.60238588  0.8627700        1.0000000         -0.8582247
Premature.Response      0.23171344  0.10754525 -0.7127425       -0.8582247          1.0000000
Missed.Response.Window  0.10515049  0.72457078 -0.8600031       -0.9865748          0.7628799
Wrong.Choice           -0.09756260 -0.74633564  0.8567408        0.9808142         -0.7416987
Missed.Response.Window Wrong.Choice
Log2FC                              0.1051505   -0.0975626
FDR                                 0.7245708   -0.7463356
MeanDiff                           -0.8600031    0.8567408
Correct.Response                   -0.9865748    0.9808142
Premature.Response                  0.7628799   -0.7416987
Missed.Response.Window              1.0000000   -0.9994830
Wrong.Choice                       -0.9994830    1.0000000

Correlation Matrix for Treatment: t2 
Log2FC        FDR   MeanDiff Correct.Response Premature.Response
Log2FC                  1.0000000 -0.4941320  0.8717175        0.1266738         0.23847919
FDR                    -0.4941320  1.0000000 -0.1990174        0.3681821        -0.55686313
MeanDiff                0.8717175 -0.1990174  1.0000000        0.2966368         0.02825100
Correct.Response        0.1266738  0.3681821  0.2966368        1.0000000        -0.49352669
Premature.Response      0.2384792 -0.5568631  0.0282510       -0.4935267         1.00000000
Missed.Response.Window -0.4449876 -0.1113784 -0.6488368       -0.8496714         0.50202739
Wrong.Choice            0.1613423 -0.3164743  0.1669027       -0.6785632        -0.06544347
Missed.Response.Window Wrong.Choice
Log2FC                             -0.4449876   0.16134232
FDR                                -0.1113784  -0.31647432
MeanDiff                           -0.6488368   0.16690272
Correct.Response                   -0.8496714  -0.67856320
Premature.Response                  0.5020274  -0.06544347
Missed.Response.Window              1.0000000   0.23595759
Wrong.Choice                        0.2359576   1.00000000

Correlation Matrix for Treatment: t1 
Log2FC        FDR    MeanDiff Correct.Response Premature.Response
Log2FC                  1.0000000 -0.3588356  0.87406928       -0.4727483          0.7205502
FDR                    -0.3588356  1.0000000 -0.45746686        0.9414830         -0.2906248
MeanDiff                0.8740693 -0.4574669  1.00000000       -0.5838467          0.7047408
Correct.Response       -0.4727483  0.9414830 -0.58384667        1.0000000         -0.4170755
Premature.Response      0.7205502 -0.2906248  0.70474079       -0.4170755          1.0000000
Missed.Response.Window  0.3149475 -0.9276952  0.44190161       -0.9538535          0.1299634
Wrong.Choice           -0.1645413 -0.7256690 -0.04120578       -0.6769218         -0.3788958
Missed.Response.Window Wrong.Choice
Log2FC                              0.3149475  -0.16454127
FDR                                -0.9276952  -0.72566899
MeanDiff                            0.4419016  -0.04120578
Correct.Response                   -0.9538535  -0.67692180
Premature.Response                  0.1299634  -0.37889582
Missed.Response.Window              1.0000000   0.85367964
Wrong.Choice                        0.8536796   1.00000000

Correlation Matrix for Treatment: t3 
Log2FC        FDR   MeanDiff Correct.Response Premature.Response
Log2FC                  1.0000000  0.8426683  0.4925862       -0.7203504          0.7203504
FDR                     0.8426683  1.0000000  0.8576563       -0.9803231          0.9803231
MeanDiff                0.4925862  0.8576563  1.0000000       -0.9259175          0.9259175
Correct.Response       -0.7203504 -0.9803231 -0.9259175        1.0000000         -1.0000000
Premature.Response      0.7203504  0.9803231  0.9259175       -1.0000000          1.0000000
Missed.Response.Window -0.7203504 -0.9803231 -0.9259175        1.0000000         -1.0000000
Wrong.Choice            0.7203504  0.9803231  0.9259175       -1.0000000          1.0000000
Missed.Response.Window Wrong.Choice
Log2FC                             -0.7203504    0.7203504
FDR                                -0.9803231    0.9803231
MeanDiff                           -0.9259175    0.9259175
Correct.Response                    1.0000000   -1.0000000
Premature.Response                 -1.0000000    1.0000000
Missed.Response.Window              1.0000000   -1.0000000
Wrong.Choice                       -1.0000000    1.0000000









##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################





## Scatterplots, Pearson's R, and Spearman correlations by treatment group


## Run for each cell type



setwd("/Volumes/DataBox/MCS2023/Stats/Lasso_CellTypes")



# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)
library(gridExtra)


# Load data
behavior_data <- read_csv("TaskA-45_AllStats.csv")
microglia_data <- read_csv("microglia_markers_bySample_2025-01-28.csv")

# Clean column names (remove spaces and special characters)
colnames(behavior_data) <- make.names(colnames(behavior_data))

# Ensure Sample IDs match (remove trailing underscores in gene data)
microglia_data <- microglia_data %>% mutate(Sample = gsub("_$", "", Sample))

# Merge the datasets
merged_data <- merge(behavior_data, microglia_data, by="Sample", all=FALSE)

# Create a long format for behavioral tests
behavior_long <- merged_data %>%
  pivot_longer(
    cols = c("Correct.Response", "Premature.Response", 
             "Missed.Response.Window", "Wrong.Choice"),
    names_to = "behavioral_test",
    values_to = "behavior_score"
  )

# Calculate overall correlation results (across all treatments)
overall_correlation_results <- behavior_long %>%
  group_by(behavioral_test) %>%
  summarise(
    Treatment = "Overall",
    n = n(),
    pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
    pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
    spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
    spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
    .groups = "drop"
  )

# Calculate correlation results by treatment group
treatment_correlation_results <- behavior_long %>%
  group_by(behavioral_test, Treatment) %>%
  summarise(
    n = n(),
    pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
    pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
    spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
    spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
    .groups = "drop"
  )

# Combine overall and treatment-specific results
correlation_results <- bind_rows(overall_correlation_results, treatment_correlation_results) %>%
  mutate(across(c(pearson_r, pearson_p, spearman_rho, spearman_p), 
                ~round(., digits = 3)))


# Create clean plot without text annotations
plot <- ggplot(behavior_long, aes(x = Log2FC, y = behavior_score, color = Treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  facet_wrap(~behavioral_test, scales = "free") +
  scale_color_viridis_d() +
  theme_bw() +
  labs(
    x = "Gene Accessibility Log2FC Score",
    y = "Behavioral Test Score",
    color = "Treatment Group",
    title = "Task 45 Microglia Gene Accessibility and Behavioral Scores"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgray"),
    panel.grid.minor = element_blank()
  )

# Create summary table for display
summary_table <- correlation_results %>%
  arrange(behavioral_test, Treatment) %>%
  mutate(
    # Format correlation and p-value together
    correlation = sprintf("r = %.3f (p = %.3f)", pearson_r, pearson_p),
    # Add asterisks for significance
    correlation = ifelse(pearson_p < 0.05, paste0(correlation, "*"), correlation)
  ) %>%
  select(behavioral_test, Treatment, n, correlation) %>%
  # Reshape to wide format for display
  pivot_wider(
    names_from = Treatment,
    values_from = correlation,
    id_cols = behavioral_test
  )

# Convert summary table to grob for plotting
table_grob <- tableGrob(
  summary_table,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 8,
    padding = unit(c(2, 2), "mm"),
    core = list(fg_params = list(hjust = 0, x = 0.1))
  )
)

# Arrange plot and table together
combined_plot <- grid.arrange(
  plot, 
  table_grob,
  heights = c(0.7, 0.3),  # Adjust these values to change relative sizes
  nrow = 2
)

# Save the combined figure
#ggsave("gaba_Task45_correlationPlot_2025-01-30.png", combined_plot, width = 12, height = 10, dpi = 300)
ggsave("microglia_Task45_correlationPlot_2025-01-30.pdf", combined_plot, width = 12, height = 10)

# Add note about approximation to the correlation results
correlation_results <- correlation_results %>%
  mutate(
    Note = case_when(
      !is.na(spearman_p) ~ "Spearman p-values are approximate due to ties in the data",
      TRUE ~ NA_character_
    )
  )

# Save detailed correlation results with note
write.csv(
  correlation_results %>%
    arrange(behavioral_test, Treatment) %>%
    select(
      behavioral_test,
      Treatment,
      n,
      pearson_r,
      pearson_p,
      spearman_rho,
      spearman_p,
      Note
    ),
  "microglia_Task45_correlationResults_2025-01-30.csv",
  row.names = FALSE
)

# Print summary to console with note about approximation
cat("Note: Spearman correlation p-values are approximate due to ties in the data\n\n")
print(correlation_results %>%
        arrange(behavioral_test, Treatment) %>%
        select(behavioral_test, Treatment, n,
               pearson_r, pearson_p,
               spearman_rho, spearman_p))



##################################
##################################
##################################
##################################


## Loop for the above code using AllStats files


# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)
library(gridExtra)
library(readr)


# Define file lists
behavior_files <- list.files(pattern = "Task.*_AllStats\\.csv$")
cell_type_files <- list.files(pattern = ".*_markers_bySample_2025-01-28\\.csv$")

# Loop through all behavioral test files and cell type marker files
for (behavior_file in behavior_files) {
  for (cell_file in cell_type_files) {
    
    # Extract test and cell type names from file names
    test_name <- gsub("_AllStats\\.csv$", "", behavior_file)
    cell_type <- gsub("_markers_bySample_2025-01-28\\.csv$", "", cell_file)
    
    # Load data
    behavior_data <- read_csv(behavior_file)
    cell_data <- read_csv(cell_file)
    
    # Clean column names
    colnames(behavior_data) <- make.names(colnames(behavior_data))
    cell_data <- cell_data %>% mutate(Sample = gsub("_$", "", Sample))  # Remove trailing underscores
    
    # Merge datasets
    merged_data <- merge(behavior_data, cell_data, by="Sample", all=FALSE)
    
    # Reshape behavioral data
    behavior_long <- merged_data %>%
      pivot_longer(
        cols = c("Correct.Response", "Premature.Response", 
                 "Missed.Response.Window", "Wrong.Choice"),
        names_to = "behavioral_test",
        values_to = "behavior_score"
      )
    
    # Compute correlations (overall)
    overall_correlation_results <- behavior_long %>%
      group_by(behavioral_test) %>%
      summarise(
        Treatment = "Overall",
        n = n(),
        pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
        pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
        spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
        spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
        .groups = "drop"
      )
    
    # Compute correlations by treatment
    treatment_correlation_results <- behavior_long %>%
      group_by(behavioral_test, Treatment) %>%
      summarise(
        n = n(),
        pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
        pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
        spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
        spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
        .groups = "drop"
      )
    
    # Combine results
    correlation_results <- bind_rows(overall_correlation_results, treatment_correlation_results) %>%
      mutate(across(c(pearson_r, pearson_p, spearman_rho, spearman_p), ~round(., digits = 3)))
    
    # Create scatter plots with regression lines
    plot <- ggplot(behavior_long, aes(x = Log2FC, y = behavior_score, color = Treatment)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      facet_wrap(~behavioral_test, scales = "free") +
      scale_color_viridis_d() +
      theme_bw() +
      labs(
        x = "Gene Accessibility Log2FC Score",
        y = "Behavioral Test Score",
        color = "Treatment Group",
        title = paste(test_name, "-", cell_type, "Gene Accessibility and Behavioral Scores")
      ) +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),
        panel.grid.minor = element_blank()
      )
    
    # Create summary table for display
    summary_table <- correlation_results %>%
      arrange(behavioral_test, Treatment) %>%
      mutate(
        correlation = sprintf("r = %.3f (p = %.3f)", pearson_r, pearson_p),
        correlation = ifelse(pearson_p < 0.05, paste0(correlation, "*"), correlation)
      ) %>%
      select(behavioral_test, Treatment, n, correlation) %>%
      pivot_wider(
        names_from = Treatment,
        values_from = correlation,
        id_cols = behavioral_test
      )
    
    # Convert summary table to grob for plotting
    table_grob <- tableGrob(
      summary_table,
      rows = NULL,
      theme = ttheme_minimal(
        base_size = 8,
        padding = unit(c(2, 2), "mm"),
        core = list(fg_params = list(hjust = 0, x = 0.1))
      )
    )
    
    # Arrange plot and table together
    combined_plot <- grid.arrange(
      plot, 
      table_grob,
      heights = c(0.7, 0.3),
      nrow = 2
    )
    
    # Save plot
    plot_filename <- paste0(cell_type, "_", test_name, "_AllScores_correlationPlot_2025-01-30.pdf")
    ggsave(plot_filename, combined_plot, width = 12, height = 10)
    
    # Add note about approximation to correlation results
    correlation_results <- correlation_results %>%
      mutate(
        Note = case_when(
          !is.na(spearman_p) ~ "Spearman p-values are approximate due to ties in the data",
          TRUE ~ NA_character_
        )
      )
    
    # Save correlation results
    csv_filename <- paste0(cell_type, "_", test_name, "_AllScores_correlationResults_2025-01-30.csv")
    write.csv(
      correlation_results %>%
        arrange(behavioral_test, Treatment) %>%
        select(behavioral_test, Treatment, n,
               pearson_r, pearson_p,
               spearman_rho, spearman_p, Note),
      csv_filename,
      row.names = FALSE
    )
    
    # Print note to console
    cat("\n", test_name, "-", cell_type, "Results Saved to", csv_filename, "\n")
  }
}




##################################
##################################
##################################




## Loop for the above code using WrongStats files


# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(viridis)
library(gridExtra)
library(readr)


# Define file lists
behavior_files <- list.files(pattern = "Task.*_WrongStats\\.csv$")
cell_type_files <- list.files(pattern = ".*_markers_bySample_2025-01-28\\.csv$")

# Loop through all behavioral test files and cell type marker files
for (behavior_file in behavior_files) {
  for (cell_file in cell_type_files) {
    
    # Extract test and cell type names from file names
    test_name <- gsub("_WrongStats\\.csv$", "", behavior_file)
    cell_type <- gsub("_markers_bySample_2025-01-28\\.csv$", "", cell_file)
    
    # Load data
    behavior_data <- read_csv(behavior_file)
    cell_data <- read_csv(cell_file)
    
    # Clean column names
    colnames(behavior_data) <- make.names(colnames(behavior_data))
    cell_data <- cell_data %>% mutate(Sample = gsub("_$", "", Sample))  # Remove trailing underscores
    
    # Merge datasets
    merged_data <- merge(behavior_data, cell_data, by="Sample", all=FALSE)
    
    # Reshape behavioral data
    behavior_long <- merged_data %>%
      pivot_longer(
        cols = c("Premature.Response", "Missed.Response.Window", "Wrong.Choice"),
        names_to = "behavioral_test",
        values_to = "behavior_score"
      )
    
    # Compute correlations (overall)
    overall_correlation_results <- behavior_long %>%
      group_by(behavioral_test) %>%
      summarise(
        Treatment = "Overall",
        n = n(),
        pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
        pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
        spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
        spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
        .groups = "drop"
      )
    
    # Compute correlations by treatment
    treatment_correlation_results <- behavior_long %>%
      group_by(behavioral_test, Treatment) %>%
      summarise(
        n = n(),
        pearson_r = cor.test(Log2FC, behavior_score, method = "pearson")$estimate,
        pearson_p = cor.test(Log2FC, behavior_score, method = "pearson")$p.value,
        spearman_rho = cor.test(Log2FC, behavior_score, method = "spearman")$estimate,
        spearman_p = cor.test(Log2FC, behavior_score, method = "spearman")$p.value,
        .groups = "drop"
      )
    
    # Combine results
    correlation_results <- bind_rows(overall_correlation_results, treatment_correlation_results) %>%
      mutate(across(c(pearson_r, pearson_p, spearman_rho, spearman_p), ~round(., digits = 3)))
    
    # Create scatter plots with regression lines
    plot <- ggplot(behavior_long, aes(x = Log2FC, y = behavior_score, color = Treatment)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      facet_wrap(~behavioral_test, scales = "free") +
      scale_color_viridis_d() +
      theme_bw() +
      labs(
        x = "Gene Accessibility Log2FC Score",
        y = "Behavioral Test Score",
        color = "Treatment Group",
        title = paste(test_name, "-", cell_type, "Gene Accessibility and Wrong Behavioral Scores")
      ) +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),
        panel.grid.minor = element_blank()
      )
    
    # Create summary table for display
    summary_table <- correlation_results %>%
      arrange(behavioral_test, Treatment) %>%
      mutate(
        correlation = sprintf("r = %.3f (p = %.3f)", pearson_r, pearson_p),
        correlation = ifelse(pearson_p < 0.05, paste0(correlation, "*"), correlation)
      ) %>%
      select(behavioral_test, Treatment, n, correlation) %>%
      pivot_wider(
        names_from = Treatment,
        values_from = correlation,
        id_cols = behavioral_test
      )
    
    # Convert summary table to grob for plotting
    table_grob <- tableGrob(
      summary_table,
      rows = NULL,
      theme = ttheme_minimal(
        base_size = 8,
        padding = unit(c(2, 2), "mm"),
        core = list(fg_params = list(hjust = 0, x = 0.1))
      )
    )
    
    # Arrange plot and table together
    combined_plot <- grid.arrange(
      plot, 
      table_grob,
      heights = c(0.7, 0.3),
      nrow = 2
    )
    
    # Save plot
    plot_filename <- paste0(cell_type, "_", test_name, "_WrongScores_correlationPlot_2025-01-30.pdf")
    ggsave(plot_filename, combined_plot, width = 12, height = 10)
    
    # Add note about approximation to correlation results
    correlation_results <- correlation_results %>%
      mutate(
        Note = case_when(
          !is.na(spearman_p) ~ "Spearman p-values are approximate due to ties in the data",
          TRUE ~ NA_character_
        )
      )
    
    # Save correlation results
    csv_filename <- paste0(cell_type, "_", test_name, "_WrongScores_correlationResults_2025-01-30.csv")
    write.csv(
      correlation_results %>%
        arrange(behavioral_test, Treatment) %>%
        select(behavioral_test, Treatment, n,
               pearson_r, pearson_p,
               spearman_rho, spearman_p, Note),
      csv_filename,
      row.names = FALSE
    )
    
    # Print note to console
    cat("\n", test_name, "-", cell_type, "Results Saved to", csv_filename, "\n")
  }
}





##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################


## Additional scatterplot code...






setwd("/Volumes/DataBox/MCS2023/Stats/Lasso_CellTypes")


# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)

# Load data
behavior_data <- read_csv("TaskA-45_AllStats.csv")
gene_data <- read_csv("gaba_markers_bySample_2025-01-28.csv")

# Clean column names
colnames(behavior_data) <- make.names(colnames(behavior_data))

# Clean up sample IDs
gene_data <- gene_data %>% mutate(Sample = gsub("_$", "", Sample))

# Merge the datasets by Sample column
merged_data <- merge(gene_data, behavior_data, by.x = "Sample", by.y = "Sample")

# Convert columns to numeric (if they aren't already)
behavioral_score <- as.numeric(merged_data$`Correct.Response`)  # Adjusting the column name for proper access
gene_accessibility_score <- as.numeric(merged_data$Log2FC)

# Check if conversion worked
str(behavioral_score)
str(gene_accessibility_score)

# Pearson correlation test
pearson_result <- cor.test(behavioral_score, gene_accessibility_score, method = "pearson")
print(pearson_result)

# Spearman correlation test (if you think the relationship is non-linear)
spearman_result <- cor.test(behavioral_score, gene_accessibility_score, method = "spearman")
print(spearman_result)


# Scatter plot with regression line for Pearson's correlation
ggplot(merged_data, aes(x = `Correct.Response`, y = Log2FC)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Correlation between Correct Response and Gene Accessibility",
       x = "Correct Response",
       y = "Gene Accessibility (Log2FC)") +
  theme_minimal()



#################

## Updates the above code to put behavioral result on y-axis, gene score on x-axis,
# point color by gene name, and point shape by treatment group: 
















#################



## Subsetting by treatment group results in: 
# Error in cor.test.default(behavioral_score_group1, gene_accessibility_score_group1,  : 
# not enough finite observations


# Subset merged data by treatment group
treatment_group <- merged_data[merged_data$Treatment == "Group1", ]
# Then perform correlation on the subset
behavioral_score_group1 <- as.numeric(treatment_group$`Correct.Response`)
gene_accessibility_score_group1 <- as.numeric(treatment_group$Log2FC)

# Pearson correlation for this group
cor.test(behavioral_score_group1, gene_accessibility_score_group1, method = "pearson")

















## Code results in Pearson and Spearman correlations, but only for 1 gene; plots look terrible


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Load data
behavior_data <- read_csv("TaskA-45_AllStats.csv")
gene_data <- read_csv("gaba_markers_bySample_2025-01-28.csv")

# Clean column names
colnames(behavior_data) <- make.names(colnames(behavior_data))

# Clean up sample IDs
gene_data <- gene_data %>% mutate(Sample = gsub("_$", "", Sample))

# Merge the datasets by Sample column
merged_data <- merge(gene_data, behavior_data, by.x = "Sample", by.y = "Sample")

# List of behavioral response columns to correlate with
behavior_columns <- c("Response.Trials", "Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")

# Create a list to store the correlation results
cor_results <- data.frame(Gene = character(),
                          TreatmentGroup = character(),
                          Behavior = character(),
                          Pearson_r = numeric(),
                          Spearman_rho = numeric())

# For each treatment group and each gene, calculate correlations with Log2FC and FDR
for (treatment in unique(merged_data$TreatmentGroup)) {
  treatment_data <- merged_data %>% filter(TreatmentGroup == treatment)
  
  for (gene in unique(treatment_data$name)) {  # Loop over each gene
    gene_data <- treatment_data %>% filter(name == gene)  # Subset for the current gene
    
    # For each behavioral column, calculate Pearson and Spearman correlations with Log2FC and FDR
    for (behavior in behavior_columns) {
      
      # Remove rows where Log2FC, FDR, or the behavioral column has NA or non-numeric values
      valid_data <- gene_data %>%
        filter(!is.na(Log2FC) & !is.na(FDR) & !is.na(gene_data[[behavior]]))
      
      if (nrow(valid_data) > 1) {  # Ensure there is more than 1 valid data point for correlation
        # Pearson correlation with Log2FC
        pearson_log2fc <- cor(valid_data$Log2FC, valid_data[[behavior]], method = "pearson", use = "complete.obs")
        
        # Spearman correlation with Log2FC
        spearman_log2fc <- cor(valid_data$Log2FC, valid_data[[behavior]], method = "spearman", use = "complete.obs")
        
        # Pearson correlation with FDR
        pearson_fdr <- cor(valid_data$FDR, valid_data[[behavior]], method = "pearson", use = "complete.obs")
        
        # Spearman correlation with FDR
        spearman_fdr <- cor(valid_data$FDR, valid_data[[behavior]], method = "spearman", use = "complete.obs")
        
        # Append the results for each gene and behavior
        cor_results <- rbind(cor_results, 
                             data.frame(Gene = gene,
                                        TreatmentGroup = treatment, 
                                        Behavior = paste(behavior, "- Log2FC"),
                                        Pearson_r = pearson_log2fc,
                                        Spearman_rho = spearman_log2fc))
        cor_results <- rbind(cor_results, 
                             data.frame(Gene = gene,
                                        TreatmentGroup = treatment, 
                                        Behavior = paste(behavior, "- FDR"),
                                        Pearson_r = pearson_fdr,
                                        Spearman_rho = spearman_fdr))
      }
    }
  }
}

# Save correlation results to CSV
write.csv(cor_results, "correlation_results_byGene_Treatment.csv", row.names = FALSE)

# Create scatterplots for each treatment group
scatter_plots <- lapply(behavior_columns, function(behavior) {
  ggplot(merged_data, aes(x = Log2FC, y = get(behavior))) +
    geom_point(aes(color = TreatmentGroup)) +
    facet_wrap(~ TreatmentGroup) +
    labs(title = paste("Log2FC vs", behavior)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# Combine all scatterplots into a grid
facet_grid_plots <- ggarrange(plotlist = scatter_plots, ncol = 2, nrow = 3)

# Save the grid of scatterplots
ggsave("scatter_plots_byTreatment.pdf", facet_grid_plots)








####################################################################






# There is issues in running correlation stats because there is only 
# 1 treatment group response per gene score. 


################



###########################################################


## Prints gene vs behavioral test with regression lines for each tx group


# Load required packages
library(ggplot2)
library(dplyr)
library(readr)
library(viridis)

# Load data
behavior <- read_csv("TaskA-45_AllStats.csv")
genes <- read_csv("gaba_markers_bySample_2025-01-28.csv")

# Clean column names (remove spaces and special characters)
colnames(behavior) <- make.names(colnames(behavior))

# Ensure Sample IDs match (remove trailing underscores in gene data)
genes <- genes %>% mutate(Sample = gsub("_$", "", Sample))

# Merge datasets (keeping only samples that have both behavior and gene data)
merged_data <- behavior %>%
  inner_join(genes, by = "Sample")

# Convert Treatment to a factor for consistency
merged_data$Treatment <- as.factor(merged_data$Treatment)
merged_data$TreatmentGroup <- as.factor(merged_data$TreatmentGroup)  # Ensure TreatmentGroup is a factor

# Define behavior tests to plot
behavior_tests <- c("Correct.Response", "Premature.Response", 
                    "Missed.Response.Window", "Wrong.Choice")

# Loop through behavior tests to create scatter plots
for (test in behavior_tests) {
  # Create the plot for the current behavioral test
  p <- ggplot(merged_data, aes(x = Log2FC, y = .data[[test]])) +
    # Plot behavioral points in black with TreatmentGroup shape
    geom_point(aes(shape = TreatmentGroup), color = "black", size = 4, alpha = 0.8) +
    
    # Plot gene points with Viridis color palette for gene scores and TreatmentGroup shape
    geom_point(aes(color = name, shape = TreatmentGroup), size = 4, alpha = 0.8) +
    
    # Linear regression line for gene scores, using color based on treatment group
    geom_smooth(method = "lm", se = FALSE, aes(color = TreatmentGroup)) +
    
    # Apply Viridis color scale for genes
    scale_color_viridis_d(option = "plasma") +  # Color palette for genes
    
    # Shape manual for TreatmentGroup differentiation
    scale_shape_manual(values = c(16, 17, 18, 15, 8)) +  # Custom shape palette (modify as needed)
    
    theme_minimal() +
    labs(x = "Gene Score (Log2FC)", y = test,
         title = paste("Gene Expression vs", test),
         color = "Gene", shape = "Treatment Group") +  # Updated legend labels
    
    theme(legend.position = "right")
  
  # Save the plot with a unique filename
  filename <- paste0("Gene_vs_", gsub("\\.", "_", test), "_with_regression.png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  
  # Display the plot
  print(p)  # Optional: If you want to display the plot in R
  
  
  # Optionally, you can also compute correlations for each treatment group per behavior test
  cor_results_pearson <- data.frame()
  cor_results_spearman <- data.frame()
  
  for (group in unique(merged_data$TreatmentGroup)) {
    # Subset data for the current treatment group
    group_data <- merged_data %>%
      filter(TreatmentGroup == group) %>%
      select(Sample, name, Log2FC, .data[[test]]) %>%
      drop_na()  # Remove rows with missing values
    
    # Ensure there are at least two data points for correlation
    if (nrow(group_data) >= 2) {
      # Pearson correlation
      cor_pearson <- cor.test(group_data$Log2FC, group_data[[test]], method = "pearson")
      
      # Spearman correlation
      cor_spearman <- cor.test(group_data$Log2FC, group_data[[test]], method = "spearman")
      
      # Store the correlation results
      cor_results_pearson <- rbind(cor_results_pearson,
                                   data.frame(TreatmentGroup = group,
                                              BehaviorTest = test, 
                                              Correlation = cor_pearson$estimate, 
                                              P_value = cor_pearson$p.value))
      
      cor_results_spearman <- rbind(cor_results_spearman,
                                    data.frame(TreatmentGroup = group,
                                               BehaviorTest = test, 
                                               Correlation = cor_spearman$estimate, 
                                               P_value = cor_spearman$p.value))
    }
  }
  
  # Optionally, save correlation results to CSV
  write.csv(cor_results_pearson, "Pearson_Correlation_Results_by_Treatment_and_Gene.csv", row.names = FALSE)
  write.csv(cor_results_spearman, "Spearman_Correlation_Results_by_Treatment_and_Gene.csv", row.names = FALSE)
}




#######################################################################







## has some helpful hints to look at the data
## Does not create plots or run stats


# Load libraries
library(tidyverse)
library(viridis)
library(glmnet)



# Load the behavioral and gene data
behavior_data <- read.csv("TaskA-45_AllStats.csv", check.names = FALSE)
gene_data <- read.csv("glut_markers_bySample_2025-01-28.csv", check.names = FALSE)

# Check the structure of the datasets
glimpse(behavior_data)
glimpse(gene_data)

# Rename columns for consistency (replace spaces with underscores)
colnames(behavior_data) <- make.names(colnames(behavior_data))
colnames(gene_data) <- make.names(colnames(gene_data))


# Perform the join
merged_data <- inner_join(behavior_data, gene_data, by = c("Sample", "Treatment" = "TreatmentGroup"))

ggplot(merged_data, aes(x = Log2FC, y = Correct.Response, color = Treatment, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE) +
  labs(
    title = "Scatterplot of Gene Score vs Correct Responses",
    x = "Gene Score (Log2FC)",
    y = "Correct Response"
  ) +
  theme_minimal() +
  theme(legend.position = "right")



# Function to calculate Pearson's correlation for each metric
# Check for mismatched sample names
unique(behavior_data$Sample)
unique(gene_data$Sample)

# If necessary, remove trailing underscores from gene_data$Sample
gene_data <- gene_data %>%
  mutate(Sample = str_replace(Sample, "_$", ""))

# Perform the join again
merged_data <- inner_join(behavior_data, gene_data, by = c("Sample", "Treatment" = "TreatmentGroup"))

# Check if the merged data has rows
print(nrow(merged_data))

# Check if Log2FC and behavior metrics contain NA or NaN
print(summary(merged_data$Log2FC))
print(summary(merged_data$Correct.Response))

# Re-run correlation calculations
cor_results <- behavior_metrics %>%
  set_names() %>%
  map(~ cor.test(merged_data[[.x]], merged_data$Log2FC))



## Run the Lasso
# Check for variability
table(merged_data$Correct.Response)

# Fit Lasso regression if variability exists
if (length(unique(merged_data$Correct.Response)) > 1) {
  X <- model.matrix(~ Log2FC + Response.Trials + Premature.Response + Missed.Response.Window + Wrong.Choice - 1, data = merged_data)
  y <- merged_data$Correct.Response
  
  # Fit Lasso model
  lasso_model <- cv.glmnet(X, y, alpha = 1)
} else {
  message("Correct.Response has no variability, regression cannot be performed.")
}



# Save correlation results
if (exists("cor_results")) {
  cor_results_df <- map_df(cor_results, ~data.frame(
    metric = .x$data.name,
    estimate = .x$estimate,
    p.value = .x$p.value
  ))
} else {
  message("Correlation results are not available.")
}













##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################
##################################


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################


## Create scatterplots of gene vs behavioral scores, perform Pearson's R correlations, and Lasso



## Prints out a merged list of gene scores and behavioral scores



# Load libraries
library(tidyverse)
library(viridis)
library(glmnet)
library(dplyr)
library(readr)
library(stringr)


# Read the data file
gene_data <- read_csv('glut_markers_bySample_2025-01-28.csv')
gene_data$Sample <- gsub("_$", "", gene_data$Sample)

# Check for duplicate matches in gene data
gene_data_counts <- gene_data %>%
  group_by(Sample, TreatmentGroup) %>%
  summarise(count = n(), .groups = "drop")

print(gene_data_counts %>% filter(count > 1))  # Show duplicates


# Read the behavioral statistics file
behavior_data <- read_csv('TaskA-45_AllStats.csv')


# Rename columns for consistency (replace spaces with underscores)
colnames(behavior_data) <- make.names(colnames(behavior_data))
colnames(gene_data) <- make.names(colnames(gene_data))

# Check column names in both dataframes
print(colnames(behavior_data))
print(colnames(gene_data))



####################


# Ensure only one row per sample in the behavior data (if there are duplicates)
behavior_data_unique <- behavior_data %>%
  distinct(Sample, .keep_all = TRUE)  # Keep only one row per Sample

# The behavioral data will be repeated correctly for each gene data row
combined_data <- gene_data %>%
  left_join(behavior_data_unique, by = "Sample")

# Replace NA behavioral data with "No Data"
combined_data[is.na(combined_data)] <- "No Data"

# Save the combined dataset to a CSV
write_csv(combined_data, "glut_TaskA-AllStats_2025-01-28.csv")

####################
####################
####################

































# Generate scatterplots
ggplot(merged_data, aes(x = Log2FC, y = Correct.Response, color = Treatment, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE) +
  labs(
    title = "Scatterplot of Gene Score vs Correct Responses",
    x = "Gene Score (Log2FC)",
    y = "Correct Response"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


# Function to calculate Pearson's correlation for each metric
behavior_metrics <- c("Response.Trials", "Correct.Response", "Premature.Response", "Missed.Response.Window", "Wrong.Choice")

cor_results <- behavior_metrics %>%
  set_names() %>%
  map(~cor.test(merged_data[[.x]], merged_data$Log2FC))

# Print correlation results
cor_results


# Prepare the data for Lasso regression
X <- model.matrix(~ Log2FC + Response.Trials + Premature.Response + Missed.Response.Window + Wrong.Choice - 1, data = merged_data)
y <- merged_data$Correct.Response

# Fit Lasso model
lasso_model <- cv.glmnet(X, y, alpha = 1)

# Plot cross-validation results
plot(lasso_model)

# Extract the best lambda and coefficients
best_lambda <- lasso_model$lambda.min
lasso_coefficients <- coef(lasso_model, s = best_lambda)

# View significant coefficients
print(lasso_coefficients)


# Save correlation results
cor_results_df <- map_df(cor_results, ~data.frame(metric = .x$data.name, estimate = .x$estimate, p.value = .x$p.value))
write.csv(cor_results_df, "correlation_results.csv", row.names = FALSE)

# Save Lasso coefficients
lasso_coef_df <- as.data.frame(as.matrix(lasso_coefficients))
write.csv(lasso_coef_df, "lasso_coefficients.csv", row.names = TRUE)















######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################






















































######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################




## Scatterplots & Pearson's correlation stats

## CANNOT RUN WITH ONLY 1 DATAPOINT (ex: 1 gene score per treatment group)

## Individual correlation plots between each gene score and the behavioral response for all genes
## Uses z-scores (created for sorted heatmaps) 


## Prepare data for Lasso analysis


# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Lasso")

# Read and clean behavioral data
behavioral_file <- "TaskA-45_AllStats-Ave.csv"
behavioral_data <- read_csv(behavioral_file, col_names = TRUE) %>%
  rename_with(~ gsub("[ \n]", "_", .x))  # Replace spaces and newlines with underscores

# Verify the column names
print(colnames(behavioral_data))


# Read the genomic data (with treatment scores already as columns)
genomic_file <- "astro_byTx_zscores_2025-01-16.csv"
genomic_data <- read_csv(genomic_file, col_names = TRUE)

# Rename the first column to 'Gene' to reflect its content
genomic_data <- genomic_data %>% rename(Gene = ...1)

# Verify the updated column names and structure
print(colnames(genomic_data))
print(head(genomic_data))

# Reshape genomic data to long format for merging
genomic_data_long <- genomic_data %>%
  pivot_longer(cols = -Gene, names_to = "Treatment", values_to = "Score")

# Verify reshaped genomic data
print(head(genomic_data_long))



# Merge the behavioral data with genomic data by 'Treatment'
merged_data <- behavioral_data %>%
  left_join(genomic_data_long, by = "Treatment")

# Reshape merged data so that gene scores become columns
final_data <- merged_data %>%
  pivot_wider(names_from = Gene, values_from = Score)

# Check the reshaped final data
print(head(final_data))



# Generate scatter plots
# Plot each behavioral metric against the gene scores
behavioral_metrics <- colnames(behavioral_data)[-1]  # Exclude 'Treatment' column
genes <- unique(merged_data$Gene)

for (metric in behavioral_metrics) {
  for (gene in genes) {
    plot <- ggplot(merged_data %>% filter(Gene == gene), aes(x = Score, y = .data[[metric]])) +
      geom_point() +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      labs(
        title = paste("Gene Score for", gene, "vs", metric, "by Treatment"),
        x = "Gene Score",
        y = metric
      ) +
      theme_minimal()
    
    # Print or save the plot
    print(plot)
    # Uncomment below to save the plots if needed
    ggsave(filename = paste0(metric, "_vs_", gene, "_GeneScore.png"), plot = plot, width = 6, height = 4)
  }
}



# Calculate correlation between behavioral scores and gene scores
# Example for 'Correct_Response' and 'Score'
correlation <- cor(
  merged_data$Correct_Response,
  merged_data$Score,
  use = "complete.obs"
)
print(correlation)


## To save correlation results for all comparisons:

# Create an empty data frame to store the correlation results
correlation_results <- data.frame(
  Metric = character(),
  Gene = character(),
  Treatment = character(),
  Pearson_R = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the behavioral metrics, genes, and treatment groups to calculate correlation
for (metric in behavioral_metrics) {
  for (gene in genes) {
    for (treatment in unique(merged_data$Treatment)) {
      # Filter data for the specific gene and treatment
      gene_data <- merged_data %>% filter(Gene == gene & Treatment == treatment)
      
      # Calculate the Pearson R correlation
      correlation <- cor(gene_data[[metric]], gene_data$Score, use = "complete.obs")
      
      # Store the results in the data frame
      correlation_results <- correlation_results %>%
        add_row(Metric = metric, Gene = gene, Treatment = treatment, Pearson_R = correlation)
    }
  }
}

# Save the correlation results to a .csv file
write.csv(correlation_results, "astro_correlation_results_byTreatment_2025-01-27.csv", row.names = FALSE)

## CANNOT RUN WITH ONLY 1 DATAPOINT



##################################################
##################################################
##################################################



# Load necessary libraries
library(dplyr)

# Read the behavioral data
behavioral_data <- read.csv("TaskA-45_AllStats.csv")

# Read the genomic data (assuming the first column contains gene names and rest are treatment data)
genomic_data <- read.csv("astro_byTx_zscores_2025-01-16.csv")

# Reshape the genomic data: you want the treatment columns to align with the behavioral data treatment groups
genomic_data_long <- genomic_data %>%
  pivot_longer(cols = starts_with("t"), names_to = "Treatment", values_to = "GenomicData") %>%
  select(Treatment, gene_name = X, GenomicData)

# Now merge the behavioral data with the reshaped genomic data
# Note: We assume that you have columns "Sample" and "Treatment" in both datasets, and you need to merge based on Treatment
merged_data <- behavioral_data %>%
  left_join(genomic_data_long, by = "Treatment")

# Check the merged data
head(merged_data)






##################################################
##################################################
##################################################
##################################################
##################################################











# Step 5: Prepare data for lasso regression
merged_data <- na.omit(merged_data)  # Remove rows with missing values

# Gene score as the predictor
x <- as.matrix(merged_data$Score)  # Predictor (gene score)

# Loop through behavioral outcomes as dependent variables
behavioral_outcomes <- c("Response_Trials", "Correct_Response", "Premature_Response", 
                         "Missed_Response_Window", "Wrong_Choice")

# Initialize a list to store results
lasso_results <- list()

# Fit a lasso model for each behavioral outcome
for (outcome in behavioral_outcomes) {
  y <- merged_data[[outcome]]  # Dependent variable (behavioral outcome)
  
  # Check if y has variability
  if (length(unique(y)) <= 1) {
    warning(paste("Skipping", outcome, "- no variability in the data."))
    next
  }
  
  # Fit lasso regression
  set.seed(42)
  lasso_model <- cv.glmnet(x, y, alpha = 1)  # Lasso with cross-validation
  
  # Store results
  lasso_results[[outcome]] <- list(
    model = lasso_model,
    coefficients = coef(lasso_model, s = "lambda.min")
  )
  
  # Plot cross-validation results
  print(paste("Lasso regression results for", outcome))
  plot(lasso_model)
}

# Example: Access coefficients for a specific outcome
print(lasso_results$Correct_Response$coefficients)




















library(readr)
library(dplyr)
library(ggplot2)
library(viridis)

# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Lasso")



# Step 1: Read and clean behavioral data
behavioral_file <- "TaskA-45_AllStats.csv"
behavioral_data <- read_csv(behavioral_file, col_names = TRUE) %>%
  rename_with(~ gsub("\\n", " ", .x))  # Replace newline characters in column names

# Check column names
print(colnames(behavioral_data))

# Step 2: Read genomic data
genomic_file <- "astro_byTx_zscores_2025-01-16.csv"
genomic_data <- read_csv(genomic_file)

# Step 3: Inspect column names for consistency
print(colnames(genomic_data))

# Step 4: Ensure that column names match in `behavioral_data` and `genomic_data`
# Confirm that "Treatment" exists in both datasets
if (!("Treatment" %in% colnames(behavioral_data)) || !("Treatment" %in% colnames(genomic_data))) {
  stop("Column 'Treatment' is missing in one or both datasets.")
}

# Step 5: Merge behavioral and genomic data
# Do not pre-aggregate; keep data granular for correlations
merged_data <- genomic_data %>%
  left_join(behavioral_data, by = "Treatment")



## Plot the Data: Before calculating the correlation, visually inspect the 
# relationship between gene scores (Score) and the behavioral outcome (Correct Response) using scatter plots. 
## This will help you see if the relationship looks linear or if there are patterns suggesting a non-linear relationship.


# Scatterplot with behavioral outcome on the y-axis, gene score on the x-axis
# and different shapes for each treatment group

# Scatterplot for Gene Score vs Correct Response (Behavioral Outcome)
ggplot(merged_data, aes(x = Score, y = `Correct Response`, shape = Treatment, color = Treatment)) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE) +  # Apply viridis color palette
  labs(title = "Gene Score vs Correct Response (Behavioral Outcome)",
       x = "Gene Score", y = "Correct Response") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_line(color = "lightgray"))

# Scatterplot for Gene Score vs Correct Response (Behavioral Outcome) with a different theme
ggplot(merged_data, aes(x = Score, y = `Correct Response`, shape = Treatment, color = Treatment)) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE) +  # Apply viridis color palette
  labs(title = "Gene Score vs Correct Response (Behavioral Outcome)",
       x = "Gene Score", y = "Correct Response") +
  theme_dark() +  # A different theme for genomic scores
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_line(color = "darkgray"))




ggplot(merged_data, aes(x = Score, y = `Correct Response`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red") +
  labs(title = "Linear vs. Polynomial Fit")


## Calculate and Compare Correlations for Linear and Non-Linear Models: 
## You can calculate both the Pearson correlation and a Spearman rank correlation, which does not assume linearity. 
## If the Pearson correlation is significantly different from the Spearman correlation, it might indicate that the relationship is not linear.

correlation_pearson <- cor(merged_data$Score, merged_data$`Correct Response`, method = "pearson", use = "complete.obs")
correlation_spearman <- cor(merged_data$Score, merged_data$`Correct Response`, method = "spearman", use = "complete.obs")

print(correlation_pearson)
print(correlation_spearman)


# Step 6: Choose one behavioral outcome and correlate with each gene score

# The code calculates the correlation between gene scores (Score) and a behavioral outcome (Correct Response) for each gene. 
# The cor function is used, with the use = "complete.obs" argument to handle missing values by excluding them from the correlation calculation. 
# The result is the Pearson correlation coefficient for each gene.

# Example: Correct_Response
correlation_results <- merged_data %>%
  group_by(Gene) %>%
  summarize(
    correlation = cor(Score, `Correct Response`, use = "complete.obs"),
    .groups = "drop"
  )


# In Step 7, the top 10 genes with the highest absolute correlations are plotted. 
# This visualization helps highlight the strongest relationships between gene scores and the behavioral outcome.
# Thus, Pearson's correlation is the statistical method used here to examine the linear relationship between gene scores and behavioral outcomes.


# Step 7: Visualize top correlations (optional)
# Example: Plot top 10 genes with the highest absolute correlations
top_genes <- correlation_results %>%
  arrange(desc(abs(correlation))) %>%
  slice_head(n = 10)

ggplot(top_genes, aes(x = reorder(Gene, -abs(correlation)), y = correlation)) +
  geom_col(fill = "skyblue") +
  coord_flip() +
  labs(
    title = "Top 10 Genes by Correlation with Correct Response",
    x = "Gene",
    y = "Correlation"
  ) +
  theme_minimal()












##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################


## Uses the gene score matrix for the Lasso 


## Converts ArchR Project into gene matrix 


#Setup an interactive session
salloc --account=eon -t 01:00:00 --mem=256G --nodes=4 --ntasks-per-node=16

module load miniconda3/24.3.0
conda activate /apps/u/opt/linux/miniconda3/24.3.0

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
setwd("/project/eon/nboldon/MCS2023/Lasso")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", 
                             force = FALSE, showLogo = FALSE)

projMCS6
getAvailableMatrices(projMCS5)


################################
################################
################################


#peakMat <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
#proj <- filterPeaks(proj, min.cells = 10)


################################
################################
################################



geneMat <- getMatrixFromProject(projMCS6, useMatrix = "GeneScoreMatrix")



# Convert geneMat to a data frame
geneMat_df <- as.data.frame(assay(geneMat))

# Save to a CSV file
write.csv(geneMat_df, "gene_score_matrix.csv", row.names = TRUE)

# Save other metadata if necessary (e.g., row or column data)
write.csv(as.data.frame(rowData(geneMat)), "gene_score_rowData.csv", row.names = TRUE)
write.csv(as.data.frame(colData(geneMat)), "gene_score_colData.csv", row.names = TRUE)



## Load the matrix back in later


# Load the matrix
geneMat_df <- read.csv("gene_score_matrix.csv", row.names = 1)

# Convert back to a SummarizedExperiment object if needed
library(SummarizedExperiment)

# Load metadata (if needed)
rowData_df <- read.csv("gene_score_rowData.csv", row.names = 1)
colData_df <- read.csv("gene_score_colData.csv", row.names = 1)

# Recreate SummarizedExperiment object
geneMat <- SummarizedExperiment(
  assays = list(counts = as.matrix(geneMat_df)),
  rowData = rowData_df,
  colData = colData_df
)




################################
################################
################################



library(glmnet)

# Load data
mat <- read.csv("normalized_matrix.csv", row.names = 1)
response <- read.csv("phenotype_data.csv")$phenotype

# Fit Lasso model
lasso_model <- cv.glmnet(as.matrix(mat), response, alpha = 1)

# Check selected features
selected_features <- coef(lasso_model, s = "lambda.min")
print(selected_features)




##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################



