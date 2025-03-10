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
library(RColorBrewer)

#Additional setup
setwd("/Volumes/DataBox/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

projMCS7

getAvailableMatrices(projMCS7)


############################################
############################################
############################################
############################################
############################################


# Define the sample renaming (refactor to avoid repetitive gsub calls)
rename_treatment <- function(sample) {
  if (sample %in% c("C302_", "C306_", "C309_", "C318_", "C323_", "C328_", "C332_", "C337_", "C339_", "C346_", "C351_", "C353_", "C360_")) {
    return("t1")
  } else if (sample %in% c("C304_", "C308_", "C312_", "C349_", "C315_", "C321_", "C324_", "C355_", "C327_", "C330_", "C333_", "C358_", "C336_", "C342_", "C348_", "C362_")) {
    return("t2")
  } else if (sample %in% c("C305_", "C307_", "C313_", "C350_", "C316_", "C320_", "C322_", "C352_", "C325_", "C334_", "C359_", "C340_", "C341_", "C345_", "C364_")) {
    return("t3")
  } else if (sample %in% c("C301_", "C303_", "C310_", "C314_", "C319_", "C335_", "C338_", "C344_", "C354_", "C356_", "C361_", "C363_")) {
    return("t4")
  }
}



# Apply renaming function to the samples
projMCS7$treatment <- sapply(projMCS7$Sample, rename_treatment)

# Check if renaming worked correctly
unique(projMCS7$treatment)

# Define dynamic clusters
clusters <- unique(projMCS7$Clusters)

# Create a custom dynamic color palette for treatment groups
n_colors <- length(unique(projMCS7$treatment))
myColors <- brewer.pal(n = min(n_colors, 9), name = "Set1")  # Use at most 9 colors from 'Set1'

# Ensure coaccessibility has been added to the project
projMCS7 <- addCoAccessibility(
  ArchRProj = projMCS7,
  reducedDims = "Harmony"  # Ensure 'Harmony' reduced dimensions are used
)

############################################


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

############################################


# Loop through each region and cluster and generate browser tracks
for (region_name in names(regions)) {
  
  # Extract the specific region from the list
  dynamicRegion <- regions[[region_name]]
  
  for (cluster in clusters) {
    
    # Subset the project by the current cluster
    cells_in_cluster <- projMCS7$Clusters == cluster
    
    # Try-catch block to handle errors during plotting
    tryCatch({
      # Generate the browser track plot with co-accessibility
      p <- plotBrowserTrack(
        ArchRProj = projMCS7,
        groupBy = "treatment",  # Make sure 'treatment' is in the metadata
        region = dynamicRegion,
        plotSummary = c("bulkTrack", "geneTrack", "loopTrack"),  # Adding loopTrack for co-accessibility
        useMatrix = "PeakMatrix",
        features = getMarkers(
          markersPeaks, 
          cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", 
          returnGR = TRUE
        ),
        upstream = 5000,
        downstream = 5000,
        pal = myColors,  # Apply the dynamic color palette
        loops = getCoAccessibility(
          ArchRProj = projMCS7,
          corCutOff = 0.5,  # Co-accessibility correlation cutoff
          returnLoops = TRUE
        )
      )
      
      # Save the plot if successful
      plotPDF(
        p, 
        name = paste0("Cluster_", cluster, "_", region_name, "_Coaccessibility_Tracks_", Sys.Date(), ".pdf"), 
        width = 6, height = 6,  # Adjust plot size for better readability
        ArchRProj = projMCS7, 
        addDOC = FALSE
      )
      
    }, error = function(e) {
      # Log error messages with cluster and region info
      message(paste("Error in cluster:", cluster, "for region:", region_name, " - ", e$message))
    })
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
############################################
############################################


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

############################################


#https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A56694976%2D56714605&hgsid=2352825847_bMvtNckSXYjKfv8T5fehlaLw8yn0

chr16:97322494-97322994
chr9:43270867-43271367
chr16:17576434-17576934 #Slc7a4
chr17:7169635-7170135
chr16:91044363-91044863 #Paxbp1 #4931406F06Rik
chr16:94411074-94411574 #Ttc3
chr17:7169635-7170135
chr16:90143656-90144156
chr4:135831661-135832161
chr16:91728462-91728962 #Atp5o #Cryzl1
chr16:84834699-84835199 #Atp5j
chr14:27562549-27563049


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

#####################

# Specify which treatment group each sample is in:
# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C8ArchRSubset$Sample

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
C8ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C8ArchRSubset$treatment)


############################################
############################################


ArchRBrowser(
  ArchRProj = C8ArchRSubset,
  features = getPeakSet(C8ArchRSubset),
  loops = getCoAccessibility(C8ArchRSubset),
  minCells = 25,
  baseSize = 10,
  borderWidth = 0.5,
  tickWidth = 0.5,
  facetbaseSize = 12,
  geneAnnotation = getGeneAnnotation(C8ArchRSubset),
  browserTheme = "cosmo",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("ArchRBrowser")
)


## Function does not accept specific region arguments.


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


## To check the coverage depth across different regions in ArchR, you can use a few approaches, 
# including generating coverage plots or directly extracting fragment information from specific genomic regions.

## 1. Using plotBrowserTrack to Visualize Coverage:
# ArchR's plotBrowserTrack can be used to visualize the fragment coverage across different genomic regions. 
# The bulk track in this plot reflects the coverage, which shows how many fragments align to each genomic position.


# You can loop through your regions to plot the coverage and inspect them visually:

# Loop through regions to check coverage
for (region_name in names(regions)) {
  region <- regions[[region_name]]
  
  # Try to plot the browser track for each region
  tryCatch({
    p <- plotBrowserTrack(
      ArchRProj = projMCS7,
      groupBy = "treatment", # Grouping by treatment for comparison
      region = region,
      plotSummary = "bulkTrack", # This will display the coverage (bulkTrack)
      useMatrix = "PeakMatrix"   # Use the peak matrix for plotting accessibility
    )
    
    # Save the plot as a PDF for inspection
    plotPDF(p, name = paste0("Coverage_", region_name, "_", Sys.Date(), ".pdf"), ArchRProj = projMCS7)
    
  }, error = function(e) {
    message(paste("Error in plotting for region:", region_name, " - ", e$message))
  })
}


# The plotBrowserTrack will create coverage plots where you can visually inspect the depth of coverage (as shown in the "bulkTrack").
# This method is useful for getting an overall sense of how coverage differs across regions or between treatments.



############################################
############################################
############################################


## 2. Extracting Fragment Counts for Regions:
# You can directly extract the fragment counts for the regions of interest using getCounts() or getGroupBW(). 
# This method gives you a more quantitative view of the coverage in specific regions.

## Extracting Coverage Data:

# Define the regions you want to check coverage for
regions_list <- regions

# Loop through regions and extract fragment counts
for (region_name in names(regions_list)) {
  region <- regions_list[[region_name]]
  
  # Extract the coverage for each region
  coverage_data <- getGroupBW(
    ArchRProj = projMCS7,
    groupBy = "treatment",
    useMatrix = "PeakMatrix", # Choose the appropriate matrix (e.g., "PeakMatrix")
    region = region,
    tileSize = 500,           # Adjust tile size as needed for better resolution
    method = "Simple"         # "Simple" to directly count fragments in each tile
  )
  
  # Output or save the coverage data for each region
  print(paste0("Coverage for region ", region_name, ":"))
  print(coverage_data)
  
  # You can also save this data for further analysis or plotting
}


# Tile Size: Adjust the tileSize parameter to change the resolution of the coverage. Smaller tiles give more fine-grained coverage but can be slower to process.
# Method: You can choose from different methods, such as "Simple" for basic counts or "Summit" if you're interested in peak summits.


############################################
############################################
############################################


## 3. Using the addGroupCoverages Function:
# You can also precompute the group coverages for different groups (e.g., treatments) in your project, 
# which can then be used in visualizations or exported as BigWig tracks for further analysis.


projMCS7 <- addGroupCoverages(
  ArchRProj = projMCS7,
  groupBy = "treatment",   # Group by treatment or other metadata
  force = TRUE             # Set to TRUE to overwrite any existing coverage data
)

# After adding group coverages, visualize them using plotBrowserTrack
p <- plotBrowserTrack(
  ArchRProj = projMCS7,
  groupBy = "treatment",   # Compare treatments
  region = GRanges("chr16", IRanges(84985335, 85173989)),
  plotSummary = "bulkTrack" # To visualize coverage
)



############################################
############################################
############################################


## 4. Exporting BigWig Files for Detailed Inspection:
# If you want to inspect the coverage in detail, you can export the BigWig files 
# for the different treatment groups and inspect them in genome browsers like IGV:


# Export BigWig files for each treatment group
projMCS7 <- getGroupBW(
  ArchRProj = projMCS7,
  groupBy = "treatment", # Group by treatment for comparison
  normMethod = "nFrags"  # Normalize by fragment counts
)

# You can load the BigWig files into a genome browser like IGV to view the coverage depth
# This will generate BigWig tracks, allowing you to visualize the coverage in a genome browser with high resolution.



############################################
############################################
############################################


## 5. Checking Raw Fragment Data in ArchR:
# You can extract the actual fragment counts directly from your Arrow files, 
# allowing you to analyze the raw fragment data:


# Extract raw fragment data for the region of interest
fragments <- getFragments(
  ArchRProj = projMCS7,
  region = GRanges("chr16", IRanges(84985335, 85173989))
)

# Summarize the fragment data
summary(fragments)


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



sessionInfo()


