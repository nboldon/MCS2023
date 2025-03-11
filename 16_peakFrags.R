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
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

############################################
############################################
############################################

idxC301 <- BiocGenerics::which(projMCS7$Sample %in% "C301_")
cellsC301 <- projMCS7$cellNames[idxC301]
projMCS7[cellsC301, ]

############

# Define your region of interest
appGene <- GRanges(
  seqnames = "chr16",
  ranges = IRanges(start = 84985335, end = 85173989)
)

# Extract the peak matrix
peakMatrixC301 <- getMatrixFromProject(
  ArchRProj = projMCS7[cellsC301, ],
  useMatrix = "PeakMatrix",
  verbose = FALSE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

# Get the peak ranges
peakRanges <- rowRanges(peakMatrixC301)

# Subset the peak ranges by your region of interest
subset_peaks <- subsetByOverlaps(peakRanges, appGene)


#######

# Extract the indices from the subset_peaks
indices <- subset_peaks$idx

# Subset the PeakMatrix using the extracted indices
subset_fragments <- peakMatrix[indices, , drop = FALSE]

# Sum the fragment counts across the subset region for each sample
total_fragments <- colSums(assay(subset_fragments))

# View the total fragments for the sample of interest
total_fragments

total_fragments

######################################
######################################
######################################


## Loop for all sample statistics

# Get the unique sample names from the project
sample_names <- unique(projMCS7$Sample)

# Initialize a data frame to store the results for each sample
results <- data.frame(Sample = character(), 
                      MedianFrags = numeric(), 
                      TSSEnrichment = numeric(), 
                      NumCells = numeric(), 
                      TotalFragments = numeric(), 
                      stringsAsFactors = FALSE)

# Loop over each sample
for (i in sample_names) {
  # Get the indices of the cells corresponding to the current sample
  idxSample <- BiocGenerics::which(projMCS7$Sample %in% i)
  cellsSample <- projMCS7$cellNames[idxSample]
  
  # Subset the project for the current sample
  projSample <- projMCS7[cellsSample, ]
  
  # Calculate median fragments for the sample
  median_frags <- median(projSample$nFrags)
  
  # Extract TSS Enrichment values for the sample
  tss_enrichment <- median(projSample$TSSEnrichment)
  
  # Calculate the number of cells in the sample
  num_cells <- length(cellsSample)
  
  # Calculate the total fragment counts (sum of all fragments across cells)
  total_fragments <- sum(projSample$nFrags)
  
  # Store the results in the data frame
  results <- rbind(results, data.frame(Sample = i, 
                                       MedianFrags = median_frags, 
                                       TSSEnrichment = tss_enrichment, 
                                       NumCells = num_cells, 
                                       TotalFragments = total_fragments))
}

# Write the results to a CSV file
write.csv(results, "sample_statistics.csv", row.names = FALSE)

# Print the results to the console
print(results)


######################################
######################################
######################################


## Loop for all samples & 1 region of interest

# Initialize an empty data frame to store the results
results <- data.frame(Sample = character(), TotalFragments = numeric(), stringsAsFactors = FALSE)

# Get the unique sample names from the project
sample_names <- unique(projMCS7$Sample)

# Loop over each sample
for (i in sample_names) {
  # Get the indices of the cells corresponding to the current sample
  idxSample <- BiocGenerics::which(projMCS7$Sample %in% i)
  cellsSample <- projMCS7$cellNames[idxSample]
  
  # Subset the project for the current sample
  projSample <- projMCS7[cellsSample, ]
  
  # Define the genomic region of interest
  app <- GRanges(
    seqnames = "chr16",
    ranges = IRanges(start = 84985335, end = 85173989)
  )
  
  # Extract the peak matrix for the current sample
  peakMatrix <- getMatrixFromProject(
    ArchRProj = projSample,
    useMatrix = "PeakMatrix",
    verbose = FALSE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile(paste0("getMatrixFromProject_", i))
  )
  
  # Get the row ranges (genomic coordinates) of the peaks
  peakRanges <- rowRanges(peakMatrix)
  
  # Subset the peak ranges by the genomic region of interest
  subset_peaks <- subsetByOverlaps(peakRanges, appGene)
  
  # Get the indices of the overlapping peaks
  indices <- subset_peaks$idx
  
  # Subset the peak matrix using these indices
  subset_fragments <- peakMatrix[indices, , drop = FALSE]
  
  # Calculate the total fragments for each cell
  total_fragments <- colSums(assay(subset_fragments))
  
  # Sum the total fragments across all cells for the current sample
  sum_fragments <- sum(total_fragments)
  
  # Store the result in the data frame
  results <- rbind(results, data.frame(Sample = i, TotalFragments = sum_fragments))
}

# Write the results to a CSV file
write.csv(results, "app_fragment_counts.csv", row.names = FALSE)

######################################
######################################
######################################


## Loop for all samples & multiple regions of interest

# List of genomic regions of interest
regions <- list(
  app_84985335_85173989 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84985335, end = 85173989)),
  olig1_91269768_91271939 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91269768, end = 91271939))
)

regions <- list(
g1600002D24Rik_95609515_95677106 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95609515, end = 95677106)),
g1600002D24Rik_95620278_95707470 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95620278, end = 95707470)),
g1700010I14Rik_8988332_9008319 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8988332, end = 9008319)),
g1700029J03Rik_93396569_93458627 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93396569, end = 93458627)),
g1700102H20Rik_3557823_3559863 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3557823, end = 3559863)),
g1700122H20Rik_6895939_6901718 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6895939, end = 6901718)),
g1810053B23Rik_93343470_93353826 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93343470, end = 93353826)),
g1810053B23Rik_93343470_93358467 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93343470, end = 93358467)),
g1810053B23Rik_93343478_93352944 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93343478, end = 93352944)),
g1810053B23Rik_93344498_93353826 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93344498, end = 93353826)),
g1810053B23Rik_93345811_93359298 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93345811, end = 93359298)),
g1810053B23Rik_93346551_93359266 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93346551, end = 93359266)),
g2310034C09Rik_88758601_88759544 = GRanges(seqnames = "chr16", ranges = IRanges(start = 88758601, end = 88759544)),
g2310043M15Rik_93791890_93794656 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93791890, end = 93794656)),
g2310061N02Rik_88706925_88707717 = GRanges(seqnames = "chr16", ranges = IRanges(start = 88706925, end = 88707717)),
g2410124H12Rik_92478496_92497125 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92478496, end = 92497125)),
g3300005D01Rik_5798656_5803240 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5798656, end = 5803240)),
g3300005D01Rik_5799489_5803240 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5799489, end = 5803240)),
g3300005D01Rik_5799489_5803242 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5799489, end = 5803242)),
g4930404I05Rik_91011003_91016258 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91011003, end = 91016258)),
g4930506C21Rik_8293365_8311118 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8293365, end = 8311118)),
g4930517M08Rik_4634979_4637477 = GRanges(seqnames = "chr17", ranges = IRanges(start = 4634979, end = 4637477)),
g4930548J01Rik_4119445_4122102 = GRanges(seqnames = "chr17", ranges = IRanges(start = 4119445, end = 4122102)),
g4930556C24Rik_85818965_85826527 = GRanges(seqnames = "chr16", ranges = IRanges(start = 85818965, end = 85826527)),
g4930590A17Rik_87909025_87952009 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87909025, end = 87952009)),
g4932438H23Rik_91054631_91056364 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91054631, end = 91056364)),
g4932438H23Rik_91054631_91069378 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91054631, end = 91069378)),
g4932438H23Rik_91054631_91129858 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91054631, end = 91129858)),
g4932438H23Rik_91055274_91068907 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91055274, end = 91068907)),
g6530411M01Rik_9147718_9168022 = GRanges(seqnames = "chr17", ranges = IRanges(start = 9147718, end = 9168022)),
A930006K02Rik_91464858_91469878 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91464858, end = 91469878)),
Adamts1_85793581_85802868 = GRanges(seqnames = "chr16", ranges = IRanges(start = 85793581, end = 85802868)),
Adamts5_85857911_85900880 = GRanges(seqnames = "chr16", ranges = IRanges(start = 85857911, end = 85900880)),
App_84952420_85175010 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84952420, end = 85175010)),
App_84954190_85173462 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84954190, end = 85173462)),
App_84985090_85173744 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84985090, end = 85173744)),
Arid1b_4995073_5347656 = GRanges(seqnames = "chr17", ranges = IRanges(start = 4995073, end = 5347656)),
Atp5j_84827625_84834552 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834552)),
Atp5j_84827625_84834731 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834731)),
Atp5j_84827625_84834733 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834733)),
Atp5j_84827625_84834753 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834753)),
Atp5j_84827625_84834756 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834756)),
Atp5j_84827625_84834777 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84834777)),
Atp5j_84827625_84835380 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84835380)),
Atp5j_84827625_84835462 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84835462)),
Atp5j_84827625_84835612 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84827625, end = 84835612)),
Atp5o_91924977_91931385 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91924977, end = 91931385)),
B130034C11Rik_87495827_87503793 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87495827, end = 87503793)),
B3galt5_96014193_96098251 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96014193, end = 96098251)),
B3galt5_96037000_96098252 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96037000, end = 96098252)),
B3galt5_96051905_96098251 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96051905, end = 96098251)),
B3galt5_96059424_96098251 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96059424, end = 96098251)),
Bace2_97135120_97221329 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97135120, end = 97221329)),
Bace2_97135409_97215496 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97135409, end = 97215496)),
Bach1_87698708_87733101 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87698708, end = 87733101)),
Bach1_87698718_87733101 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87698718, end = 87733101)),
Bach1_87713884_87733101 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87713884, end = 87733101)),
Brwd1_95770481_95832319 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95770481, end = 95832319)),
Brwd1_95770481_95844987 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95770481, end = 95844987)),
Brwd1_95770481_95861096 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95770481, end = 95861096)),
Brwd1_95770484_95860821 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95770484, end = 95860821)),
Brwd1_95779940_95860821 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95779940, end = 95860821)),
Brwd1_95803107_95861096 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95803107, end = 95861096)),
Brwd1_95845060_95860821 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95845060, end = 95860821)),
C2cd2_97633601_97701026 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97633601, end = 97701026)),
C2cd2_97633606_97701003 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97633606, end = 97701003)),
C2cd2_97633606_97704188 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97633606, end = 97704188)),
C2cd2_97643146_97704188 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97643146, end = 97704188)),
C2cd2_97648679_97704188 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97648679, end = 97704188)),
Cbr3_93682973_93690746 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93682973, end = 93690746)),
Ccr6_8236042_8257129 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8236042, end = 8257129)),
Ccr6_8242525_8257129 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8242525, end = 8257129)),
Ccr6_8243806_8257129 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8243806, end = 8257129)),
Ccr6_8245057_8257129 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8245057, end = 8257129)),
Cct8_87483079_87495624 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87483079, end = 87495624)),
Cfap298_90925565_90934604 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90925565, end = 90934604)),
Chaf1b_93883654_93905870 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93883654, end = 93905870)),
Chaf1b_93883655_93905861 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93883655, end = 93905861)),
Chaf1b_93884300_93905870 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93884300, end = 93905870)),
Cldn14_93918785_93929322 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93918785, end = 93929322)),
Cldn14_93918785_93929572 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93918785, end = 93929572)),
Cldn14_93918785_93929601 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93918785, end = 93929601)),
Cldn14_93918785_94008592 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93918785, end = 94008592)),
Cldn17_88505561_88506733 = GRanges(seqnames = "chr16", ranges = IRanges(start = 88505561, end = 88506733)),
Cldn8_88560580_88562938 = GRanges(seqnames = "chr16", ranges = IRanges(start = 88560580, end = 88562938)),
Clic6_92497901_92540996 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92497901, end = 92540996)),
Cryzl1_91689054_91728960 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91689054, end = 91728960)),
Cryzl1_91689054_91728961 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91689054, end = 91728961)),
Cryzl1_91689056_91728955 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91689056, end = 91728955)),
Cryzl1_91689056_91728956 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91689056, end = 91728956)),
Cryzl1_91689076_91728557 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91689076, end = 91728557)),
Cryzl1_91694096_91728955 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91694096, end = 91728955)),
Cryzl1_91694212_91728954 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91694212, end = 91728954)),
Cyyr1_85452612_85550660 = GRanges(seqnames = "chr16", ranges = IRanges(start = 85452612, end = 85550660)),
Cyyr1_85455997_85550128 = GRanges(seqnames = "chr16", ranges = IRanges(start = 85455997, end = 85550128)),
Donson_91666610_91688687 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91666610, end = 91688687)),
Donson_91679000_91684146 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91679000, end = 91684146)),
Donson_91679019_91688483 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91679019, end = 91688483)),
Dop1b_93711661_93740033 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711661, end = 93740033)),
Dop1b_93711661_93810343 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711661, end = 93810343)),
Dop1b_93711663_93783392 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711663, end = 93783392)),
Dop1b_93711663_93810335 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711663, end = 93810335)),
Dop1b_93711663_93810338 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711663, end = 93810338)),
Dop1b_93711663_93810343 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711663, end = 93810343)),
Dop1b_93711663_93792905 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711663, end = 93792905)),
Dop1b_93711666_93810335 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711666, end = 93810335)),
Dop1b_93711669_93810335 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711669, end = 93810335)),
Dop1b_93711810_93810335 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93711810, end = 93810335)),
Dop1b_93725546_93810335 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93725546, end = 93810335)),
Dop1b_93729853_93810343 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93729853, end = 93810343)),
Dscam_96365844_96949311 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96365844, end = 96949311)),
Dscam_96370471_96949128 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96370471, end = 96949128)),
Dscam_96371307_96493821 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96371307, end = 96493821)),
Dscam_96371307_96611886 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96371307, end = 96611886)),
Dscam_96407760_96949347 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96407760, end = 96949347)),
Dscam_96452104_96949506 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96452104, end = 96949506)),
Dyrk1a_94348339_94473457 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94348339, end = 94473457)),
Dyrk1a_94348598_94473912 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94348598, end = 94473912)),
Dyrk1a_94348843_94473457 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94348843, end = 94473457)),
Dyrk1a_94348853_94473457 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94348853, end = 94473457)),
Dyrk1a_94349299_94473912 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94349299, end = 94473912)),
Dyrk1a_94384592_94473460 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94384592, end = 94473460)),
Dyrk1a_94389811_94473912 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94389811, end = 94473912)),
Erg_95137561_95197064 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95197064)),
Erg_95137561_95237162 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95237162)),
Erg_95137561_95237741 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95237741)),
Erg_95137561_95237766 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95237766)),
Erg_95137561_95308758 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95308758)),
Erg_95137561_95363914 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95363914)),
Erg_95137561_95363924 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95363924)),
Erg_95137561_95364981 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95364981)),
Erg_95137561_95364986 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95364986)),
Erg_95137561_95364991 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95364991)),
Erg_95137561_95364992 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95364992)),
Erg_95137561_95364993 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95137561, end = 95364993)),
Ets2_95480799_95499442 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95480799, end = 95499442)),
Eva1c_90826473_90904864 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90826473, end = 90904864)),
Eva1c_90826596_90904450 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90826596, end = 90904450)),
Eva1c_90829831_90904450 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90829831, end = 90904450)),
Eva1c_90829833_90904450 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90829833, end = 90904450)),
Eva1c_90830613_90904864 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90830613, end = 90904864)),
Eva1c_90831220_90904450 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90831220, end = 90904450)),
#ArchR logging to : ArchRLogs/ArchR-getMatrixFromProject_C338_-33030a38d96-Date-2024-09-07_Time-21-32-04.253929.log
#If there is an issue, please report to github with logFile!
#  Error in h5ls(ArrowFiles[x]) : HDF5. Object header. Can't get value.


regions2 <- list(
Ezr_6738130_6782780 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6738130, end = 6782780)),
Fam243_92318517_92321196 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92318517, end = 92321196)),
Fndc1_7738568_7827302 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7738568, end = 7827302)),
Gabpa_84834645_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834645, end = 84863538)),
Gabpa_84834655_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834655, end = 84863538)),
Gabpa_84834679_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834679, end = 84863538)),
Gabpa_84834817_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834817, end = 84863538)),
Gabpa_84834818_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834818, end = 84863538)),
Gabpa_84834820_84863538 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834820, end = 84863538)),
Gabpa_84834878_84863534 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84834878, end = 84863534)),
Gart_91621149_91646727 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91621149, end = 91646727)),
Gart_91621161_91646686 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91621161, end = 91646686)),
Gart_91621161_91646696 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91621161, end = 91646696)),
Gm10785_91688652_91715510 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91688652, end = 91715510)),
Gm10789_90141991_90146582 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90141991, end = 90146582)),
Gm10791_84971965_84979206 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84971965, end = 84979206)),
Gm29805_92376825_92382438 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92376825, end = 92382438)),
Gm31641_94962929_95038095 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94962929, end = 95038095)),
Gm32865_87124790_87128194 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87124790, end = 87128194)),
Gm8363_5480988_5483069 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5480988, end = 5483069)),
Grik1_87895655_88289971 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87895655, end = 88289971)),
Grik1_87895656_88289457 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87895656, end = 88289457)),
Grik1_87895656_88290618 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87895656, end = 88290618)),
Grik1_87896013_88290618 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87896013, end = 88290618)),
Grik1_87911808_88290618 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87911808, end = 88290618)),
Grik1_87912170_88289971 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87912170, end = 88289971)),
Grik1_87945344_88290618 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87945344, end = 88290618)),
Gtf2h5_6079827_6085485 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6079827, end = 6085485)),
Gtf2h5_6079920_6085485 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6079920, end = 6085485)),
Hmgn1_95899980_95906118 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95899980, end = 95906118)),
Hunk_90385038_90499310 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90385038, end = 90499310)),
Hunk_90386151_90499308 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90386151, end = 90499308)),
Hunk_90401429_90499310 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90401429, end = 90499310)),
Ifnar1_91484786_91507193 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91484786, end = 91507193)),
Ifnar1_91484969_91510278 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91484969, end = 91510278)),
Ifnar2_91372537_91393798 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91372537, end = 91393798)),
Ifnar2_91372537_91405342 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91372537, end = 91405342)),
Ifnar2_91372537_91405344 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91372537, end = 91405344)),
Ifnar2_91372548_91405337 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91372548, end = 91405337)),
Ifngr2_91546813_91563762 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91546813, end = 91563762)),
Igsf5_96139871_96200517 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96139871, end = 96200517)),
Igsf5_96139895_96243965 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96139895, end = 96243965)),
Igsf5_96139898_96200517 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96139898, end = 96200517)),
Igsf5_96139920_96200517 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96139920, end = 96200517)),
Igsf5_96140149_96182424 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96140149, end = 96182424)),
Igsf5_96140153_96200514 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96140153, end = 96200514)),
Igsf5_96140156_96182424 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96140156, end = 96182424)),
Il10rb_91405911_91425589 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91405911, end = 91425589)),
Il10rb_91405989_91425589 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91405989, end = 91425589)),
Itgb2l_96200680_96217352 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96200680, end = 96217352)),
Itgb2l_96200680_96217658 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96200680, end = 96217658)),
Itgb2l_96200680_96222057 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96200680, end = 96222057)),
Itgb2l_96200690_96222007 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96200690, end = 96222007)),
Itgb2l_96205372_96217664 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96205372, end = 96217664)),
Itsn1_91729064_91920346 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729064, end = 91920346)),
Itsn1_91729125_91832108 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729125, end = 91832108)),
Itsn1_91729125_91871410 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729125, end = 91871410)),
Itsn1_91729125_91920334 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729125, end = 91920334)),
Itsn1_91729362_91920346 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729362, end = 91920346)),
Itsn1_91729863_91853015 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729863, end = 91853015)),
Itsn1_91729863_91863073 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729863, end = 91863073)),
Itsn1_91729863_91871410 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729863, end = 91871410)),
Itsn1_91729863_91920346 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91729863, end = 91920346)),
Itsn1_91739910_91920346 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91739910, end = 91920346)),
Itsn1_91803028_91920346 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91803028, end = 91920346)),
Jam2_84773877_84826135 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84773877, end = 84826135)),
Jam2_84793099_84825687 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84793099, end = 84825687)),
Kcne1_92345755_92358629 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92345755, end = 92358629)),
Kcne1_92345755_92359223 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92345755, end = 92359223)),
Kcne2_92292143_92297888 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92292143, end = 92297888)),
Kcne2_92295370_92297888 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92295370, end = 92297888)),
Kcnj15_95035950_95078653 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95035950, end = 95078653)),
Kcnj15_95036002_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95036002, end = 95075804)),
Kcnj15_95036002_95078606 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95036002, end = 95078606)),
Kcnj15_95036003_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95036003, end = 95075804)),
Kcnj15_95036428_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95036428, end = 95075804)),
Kcnj15_95036511_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95036511, end = 95075804)),
Kcnj15_95067263_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95067263, end = 95075804)),
Kcnj15_95069204_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95069204, end = 95075804)),
Kcnj15_95069621_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95069621, end = 95075804)),
Kcnj15_95071515_95075804 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95071515, end = 95075804)),
Kcnj6_94523372_94635100 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94523372, end = 94635100)),
Kcnj6_94527075_94776089 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94527075, end = 94776089)),
Kcnj6_94539472_94776089 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94539472, end = 94776089)),
Kcnj6_94600924_94635083 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94600924, end = 94635083)),
Kcnj6_94610122_94776089 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94610122, end = 94776089)),
Krtap11_1_89569930_89570938 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89569930, end = 89570938)),
Krtap16_3_88962058_88962628 = GRanges(seqnames = "chr16", ranges = IRanges(start = 88962058, end = 88962628)),
Krtap20_2_89205615_89206149 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89205615, end = 89206149)),
Krtap21_1_89402781_89403529 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89402781, end = 89403529)),
Krtap6_2_89419077_89419866 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89419077, end = 89419866)),
Krtap7_1_89507457_89508078 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89507457, end = 89508078)),
Krtap8_1_89487128_89487707 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89487128, end = 89487707)),
Lca5l_95936798_95970650 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95936798, end = 95970650)),
Lca5l_95936799_95970628 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95936799, end = 95970628)),
Lca5l_95936799_95970650 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95936799, end = 95970650)),
Ldhal6b_5417322_5418767 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5417322, end = 5418767)),
Ltn1_87376405_87432361 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87376405, end = 87432361)),
Ltn1_87376405_87432310 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87376405, end = 87432310)),
Ltn1_87377008_87412982 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87377008, end = 87412982)),
Ltn1_87400132_87432374 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87400132, end = 87432374)),
Map3k7cl_87553084_87595091 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87553084, end = 87595091)),
Map3k7cl_87553357_87594492 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87553357, end = 87594492)),
Map3k7cl_87575347_87594492 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87575347, end = 87594492)),
Mir6964_97655625_97655684 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97655625, end = 97655684)),
Mir802_93369474_93369571 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93369474, end = 93369571)),
Mis18a_90719066_90727126 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90719066, end = 90727126)),
Morc3_93831875_93875828 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93831875, end = 93875828)),
Morc3_93832017_93875832 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93832017, end = 93875832)),
Morc3_93844087_93875832 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93844087, end = 93875832)),
Mpc1_8283786_8297667 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8283786, end = 8297667)),
Mpc1_8284453_8297667 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8284453, end = 8297667)),
Mrap_90737961_90749540 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90737961, end = 90749540)),
Mrpl39_84718035_84734999 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84718035, end = 84734999)),
Mrpl39_84718601_84734992 = GRanges(seqnames = "chr16", ranges = IRanges(start = 84718601, end = 84734992)),
Mrps6_92058090_92111982 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92058090, end = 92111982)),
Mx1_97225427_97241299 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97225427, end = 97241299)),
Mx2_97314473_97339294 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97314473, end = 97339294)),
N6amt1_87353939_87368404 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87353939, end = 87368404)),
Nox3_3635239_3696261 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3635239, end = 3696261)),
Olig1_91269523_91271694 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91269523, end = 91271694)),
Olig2_91225304_91228432 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91225304, end = 91228432)),
Paxbp1_91013791_91044134 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91013791, end = 91044134)),
Pcp4_96245998_96304186 = GRanges(seqnames = "chr16", ranges = IRanges(start = 96245998, end = 96304186)),
Pde10a_8526800_8986648 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8526800, end = 8986648)),
Pde10a_8801744_8986648 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8801744, end = 8986648)),
Pde10a_8849984_8986648 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8849984, end = 8986648)),
Pigp_94137156_94148723 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94137156, end = 94148723)),
Pigp_94142845_94148723 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94142845, end = 94148723)),
Pigp_94142845_94148943 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94142845, end = 94148943)),
Pigp_94142845_94149103 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94142845, end = 94149103)),
Pigp_94142845_94149408 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94142845, end = 94149408)),
Pisd_ps2 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3064317, end = 3084183)),
Prdm15_97569859_97629620 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569859, end = 97629620)),
Prdm15_97569859_97630243 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569859, end = 97630243)),
Prdm15_97569912_97616051 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97616051)),
Prdm15_97569912_97630057 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630057)),
Prdm15_97569912_97630214 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630214)),
Prdm15_97569912_97630230 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630230)),
Prdm15_97569912_97630236 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630236)),
Prdm15_97569912_97630259 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630259)),
Prdm15_97569912_97630267 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97569912, end = 97630267)),
Prdm15_97572917_97630268 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97572917, end = 97630268)),
Prr18_8340405_8344113 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8340405, end = 8344113)),
Prr18_8340738_8344113 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8340738, end = 8344113)),
Psmg1_95758325_95769353 = GRanges(seqnames = "chr16", ranges = IRanges(start = 95758325, end = 95769353)),
Rcan1_92391705_92399833 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92391705, end = 92399833)),
Rcan1_92391705_92465924 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92391705, end = 92465924)),
Ripk4_97520325_97541877 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97520325, end = 97541877)),
Ripk4_97520325_97542172 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97520325, end = 97542172)),
Rnaset2a_8128579_8148246 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8128579, end = 8148246)),
Rnaset2a_8128579_8147787 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8128579, end = 8147787)),
Rnaset2b_6978859_6998193 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6978859, end = 6998193)),
Rnaset2b_8128590_8147832 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8128590, end = 8147832)),
Rps6ka2_7170114_7303315 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7170114, end = 7303315)),
Rsph3a_7945613_7979556 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7945613, end = 7979556)),
Runx1_92601220_92697084 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92601220, end = 92697084)),
Runx1_92601220_92825829 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92601220, end = 92825829)),
Rwdd2b_87433085_87440347 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87433085, end = 87440347)),
Scaf4_90228704_90256991 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228704, end = 90256991)),
Scaf4_90228704_90284241 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228704, end = 90284241)),
Scaf4_90228704_90284251 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228704, end = 90284251)),
Scaf4_90228704_90285895 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228704, end = 90285895)),
Scaf4_90228705_90284241 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228705, end = 90284241)),
Scaf4_90228705_90284258 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228705, end = 90284258)),
Scaf4_90228905_90284472 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90228905, end = 90284472)),
Scaf8_3114971_3198859 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3114971, end = 3198859)),
Serac1_6040570_6079739 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6040570, end = 6079739)),
Setd4_93583211_93599036 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93583211, end = 93599036)),
Setd4_93583211_93604541 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93583211, end = 93604541)),
Setd4_93583211_93604542 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93583211, end = 93604542)),
Setd4_93583211_93604545 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93583211, end = 93604545)),
Setd4_93583214_93603570 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93583214, end = 93603570)),
Setd4_93584674_93604544 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93584674, end = 93604544)),
Setd4_93584976_93604543 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93584976, end = 93604543)),
Setd4_93589526_93604544 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93589526, end = 93604544)),
Setd4_93589747_93604546 = GRanges(seqnames = "chr16", ranges = IRanges(start = 93589747, end = 93604546)),
Sft2d1_8311087_8327442 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8311087, end = 8327442)),
Sft2d1_8311102_8327442 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8311102, end = 8327442)),
Slc5a3_92058076_92087228 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92058076, end = 92087228)),
Smim11_92300995_92312796 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92300995, end = 92312796)),
Smim11_92301057_92312796 = GRanges(seqnames = "chr16", ranges = IRanges(start = 92301057, end = 92312796)),
Snx9_5841323_5931956 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5841323, end = 5931956)),
Sod1_90220516_90226088 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90220516, end = 90226088)),
Son_91647578_91663073 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91647578, end = 91663073)),
Son_91647578_91678947 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91647578, end = 91678947)),
Synj1_90935846_90971188 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935846, end = 90971188)),
Synj1_90935846_90976348 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935846, end = 90976348)),
Synj1_90935846_91010850 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935846, end = 91010850)),
Synj1_90935846_91013482 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935846, end = 91013482)),
Synj1_90935851_91011063 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935851, end = 91011063)),
Synj1_90935855_91010850 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90935855, end = 91010850)),
Synj1_90939386_91013482 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90939386, end = 91013482)),
Synj2_5941279_6039390 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5941279, end = 6039390)),
Synj2_5941279_6040191 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5941279, end = 6040191)),
Synj2_5975579_6040191 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5975579, end = 6040191)),
Synj2_5975585_6044290 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5975585, end = 6044290)),
Synj2_6007579_6039390 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6007579, end = 6039390)),
T_8434422_8442496 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8434422, end = 8442496)),
T2_8372395_8422726 = GRanges(seqnames = "chr17", ranges = IRanges(start = 8372395, end = 8422726)),
Tagap_7925999_7934897 = GRanges(seqnames = "chr17", ranges = IRanges(start = 7925999, end = 7934897)),
Tfb1m_3519250_3557755 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3519250, end = 3557755)),
Tiam1_89786865_89818107 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89818107)),
Tiam1_89786865_89897144 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89897144)),
Tiam1_89786865_89960374 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89960374)),
Tiam1_89786865_89960582 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89960582)),
Tiam1_89786865_89960587 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89960587)),
Tiam1_89786865_89974454 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 89974454)),
Tiam1_89786865_90143148 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 90143148)),
Tiam1_89786865_90143285 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 90143285)),
Tiam1_89786865_90145355 = GRanges(seqnames = "chr16", ranges = IRanges(start = 89786865, end = 90145355)),
Tiam2_3326572_3519397 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3326572, end = 3519397)),
Tiam2_3397206_3519397 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3397206, end = 3519397)),
Tiam2_3482314_3519397 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3482314, end = 3519397)),
Tiam2_3483489_3519397 = GRanges(seqnames = "chr17", ranges = IRanges(start = 3483489, end = 3519397)),
Tmem242_5410863_5440260 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5410863, end = 5440260)),
Tmem50b_91574262_91597435 = GRanges(seqnames = "chr16", ranges = IRanges(start = 91574262, end = 91597435)),
Tmprss2_97343074_97389052 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97343074, end = 97389052)),
Tmprss2_97343074_97389588 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97343074, end = 97389588)),
Tmprss2_97343074_97437934 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97343074, end = 97437934)),
Ttc3_94149017_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149017, end = 94247615)),
Ttc3_94149037_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149037, end = 94247615)),
Ttc3_94149083_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149083, end = 94247615)),
Ttc3_94149098_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149098, end = 94247615)),
Ttc3_94149101_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149101, end = 94247615)),
Ttc3_94149112_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149112, end = 94247615)),
Ttc3_94149117_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149117, end = 94247615)),
Ttc3_94149118_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149118, end = 94247615)),
Ttc3_94149124_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149124, end = 94247615)),
Ttc3_94149131_94247614 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149131, end = 94247614)),
Ttc3_94149131_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149131, end = 94247615)),
Ttc3_94149138_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149138, end = 94247615)),
Ttc3_94149141_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149141, end = 94247615)),
Ttc3_94149143_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149143, end = 94247615)),
Ttc3_94149150_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149150, end = 94247615)),
Ttc3_94149159_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149159, end = 94247615)),
Ttc3_94149165_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149165, end = 94247615)),
Ttc3_94149177_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149177, end = 94247615)),
Ttc3_94149184_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149184, end = 94247615)),
Ttc3_94149190_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149190, end = 94247615)),
Ttc3_94149275_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149275, end = 94247615)),
Ttc3_94149284_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94149284, end = 94247615)),
Ttc3_94150254_94247615 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94150254, end = 94247615)),
Tulp4_6106829_6240637 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6106829, end = 6240637)),
Tulp4_6138147_6240637 = GRanges(seqnames = "chr17", ranges = IRanges(start = 6138147, end = 6240637)),
Urb1_90751281_90810168 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90751281, end = 90810168)),
Urb1_90751703_90810222 = GRanges(seqnames = "chr16", ranges = IRanges(start = 90751703, end = 90810222)),
Usp16_87454739_87483270 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87454739, end = 87483270)),
Usp16_87454744_87483270 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87454744, end = 87483270)),
Usp16_87454803_87483270 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87454803, end = 87483270)),
Usp16_87454891_87483270 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87454891, end = 87483270)),
Usp16_87454921_87483270 = GRanges(seqnames = "chr16", ranges = IRanges(start = 87454921, end = 87483270)),
Vps26c_94276116_94305022 = GRanges(seqnames = "chr16", ranges = IRanges(start = 94276116, end = 94305022)),
Zbtb21_97724187_97739639 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97724187, end = 97739639)),
Zbtb21_97724187_97740645 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97724187, end = 97740645)),
Zbtb21_97724187_97741019 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97724187, end = 97741019)),
Zbtb21_97724235_97741019 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97724235, end = 97741019)),
Zbtb21_97725827_97741015 = GRanges(seqnames = "chr16", ranges = IRanges(start = 97725827, end = 97741015)),
Zdhhc14_5492599_5753891 = GRanges(seqnames = "chr17", ranges = IRanges(start = 5492599, end = 5753891))
)


###########################################
###########################################

# Get the unique sample names from the project
sample_names <- unique(projMCS7$Sample)

# Loop over each genomic region
for (region_name in names(regions2)) {
  
  # Initialize an empty data frame to store the results for the current region
  results <- data.frame(Sample = character(), stringsAsFactors = FALSE)
  
  # Get the current genomic region from the list
  region <- regions2[[region_name]]
  
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
    
    # Sum the total fragments across all cells for the current sample
    sum_fragments <- sum(total_fragments)
    
    # Add the results to the data frame (adding 'TotalFragments' as the column name)
    results <- rbind(results, data.frame(Sample = i, sum_fragments))
  }
  
  # Rename the column 'sum_fragments' to reflect the region name (e.g., "appTotalFrags")
  colnames(results)[2] <- paste0(region_name, "TotalFrags")
  
  # Write the results to a CSV file named after the region
  write.csv(results, paste0(region_name, "_fragment_counts.csv"), row.names = FALSE)
}
