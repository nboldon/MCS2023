#Setup an interactive session
salloc --account=eon -t 18:00:00 --mem=64G --nodes=1 --ntasks-per-node=16

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
library(BSgenome.Mmusculus.UCSC.mm10)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ProjMCS9")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

############################################
############################################ 

#Load project
projMCS2 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS2", force = FALSE, showLogo = FALSE)

############################################
############################################

#Pseudo-bulk replicates in ArchR

projMCS9 <- addGroupCoverages(ArchRProj = projMCS2, groupBy = "Clusters") 

# 10 Calling Peaks w/ Macs2

#ProjMCS4 code:
#pathToMacs2 <- findMacs2()
#projMCS4 <- addReproduciblePeakSet(
#    ArchRProj = projMCS4, 
#    groupBy = "Clusters", 
#    pathToMacs2 = pathToMacs2
#)
#getPeakSet(projMCS4)
#saveArchRProject(ArchRProj = projMCS4, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS4", load = FALSE)


pathToMacs2 <- findMacs2()

projMCS9 <- addReproduciblePeakSet(
  ArchRProj = projMCS9, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  extsize = 30, # 20% of 150
  extendSummits = 50 # 20% of 250
)


extsize	
#The number of basepairs to extend the MACS2 fragment after shift has been applied. 
#When combined with extsize this allows you to create proper fragments, 
#centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).

extendSummits	
#The number of basepairs to extend peak summits (in both directions) to obtain final fixed-width peaks. 
#For example, extendSummits = 250 will create 501-bp fixed-width peaks from the 1-bp summits.

getPeakSet(projMCS9)

projMCS9 <- addPeakMatrix(projMCS9)

getAvailableMatrices(projMCS9)

saveArchRProject(ArchRProj = projMCS9, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS9", load = FALSE)


############################################
############################################

#Load project
projMCS9 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS9", force = FALSE, showLogo = FALSE)

############################################
############################################

markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS9,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
markerList

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
  archRSubset <- projMCS9[projMCS9$Clusters %in% cluster]
  
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
      useMatrix = "PeakMatrix",
      groupBy = "treatment",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = useGroup,
      bgdGroups = bgdGroup
    )
    
    # Get the list of markers
    markerList <- getMarkers(markerFeatures, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")
    
    # Define file name
    file_name <- paste0(cluster, "_", useGroup, "vs", bgdGroup, "_Peak9_FDR-0-1_Log2FC-0-5_2024-09-12.csv")
    
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
############################################
############################################
############################################
############################################


#Pseudo-bulk replicates in ArchR

projMCS9_5 <- addGroupCoverages(ArchRProj = projMCS2, groupBy = "Clusters")

################


pathToMacs2 <- findMacs2()

projMCS9_5 <- addReproduciblePeakSet(
  ArchRProj = projMCS9_5,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  extsize = 270, # 80% more than 150
  extendSummits = 450 # 80% more than 250
)


getPeakSet(projMCS9_5)

projMCS9_5 <- addPeakMatrix(projMCS9_5)

getAvailableMatrices(projMCS9_5)

saveArchRProject(ArchRProj = projMCS9_5, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS9_5", load = FALSE)


############################################
############################################

#Load project
projMCS9_5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS9_5", force = FALSE, showLogo = FALSE)

############################################
############################################

markersPeaks <- getMarkerFeatures(
    ArchRProj = projMCS9_5,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
markerList

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
  archRSubset <- projMCS9_5[projMCS9_5$Clusters %in% cluster]

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
      useMatrix = "PeakMatrix",
      groupBy = "treatment",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = useGroup,
      bgdGroups = bgdGroup
    )

    # Get the list of markers
    markerList <- getMarkers(markerFeatures, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

    # Define file name
    file_name <- paste0(cluster, "_", useGroup, "vs", bgdGroup, "_Peak9-5_FDR-0-1_Log2FC-0-5_2024-09-13.csv")

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
############################################
############################################
############################################
############################################


#Pseudo-bulk replicates in ArchR

projMCS9_7 <- addGroupCoverages(ArchRProj = projMCS2, groupBy = "Clusters")

################


pathToMacs2 <- findMacs2()

projMCS9_7 <- addReproduciblePeakSet(
  ArchRProj = projMCS9_7,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  extsize = 75,  # Intermediate between default and reduced; default is 150
  extendSummits = 125  # Intermediate size; default is 250
)

getPeakSet(projMCS9_7)

projMCS9_7 <- addPeakMatrix(projMCS9_7)

getAvailableMatrices(projMCS9_7)

saveArchRProject(ArchRProj = projMCS9_7, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS9_7", load = FALSE)


############################################
############################################

#Load project
projMCS9_7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS9_7", force = FALSE, showLogo = FALSE)

############################################
############################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = projMCS9_7,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
markerList

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
  archRSubset <- projMCS9_7[projMCS9_7$Clusters %in% cluster]
  
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
      useMatrix = "PeakMatrix",
      groupBy = "treatment",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = useGroup,
      bgdGroups = bgdGroup
    )
    
    # Get the list of markers
    markerList <- getMarkers(markerFeatures, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")
    
    # Define file name
    file_name <- paste0(cluster, "_", useGroup, "vs", bgdGroup, "_Peak9-7_FDR-0-1_Log2FC-0-5_2024-09-14.csv")
    
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
############################################
############################################


sessionInfo()

R version 4.3.2 (2023-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.9 (Ootpa)

Matrix products: default
BLAS/LAPACK: /pfs/tc1/home/nboldon/.conda/envs/archr2023_12/lib/libopenblasp-r0.3.25.so;  LAPACK version 3.11.0

Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Denver
tzcode source: system (glibc)

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] presto_1.0.0                nabor_0.5.0                
 [3] pheatmap_1.0.12             enrichplot_1.20.3          
 [5] clusterProfiler_4.8.3       org.Mm.eg.db_3.17.0        
 [7] AnnotationDbi_1.64.1        rhdf5_2.44.0               
 [9] SummarizedExperiment_1.32.0 Biobase_2.62.0             
[11] MatrixGenerics_1.14.0       Rcpp_1.0.12                
[13] Matrix_1.6-4                GenomicRanges_1.54.1       
[15] GenomeInfoDb_1.38.5         IRanges_2.36.0             
[17] S4Vectors_0.40.2            BiocGenerics_0.48.1        
[19] matrixStats_1.2.0           data.table_1.14.10         
[21] stringr_1.5.1               plyr_1.8.9                 
[23] magrittr_2.0.3              ggplot2_3.4.4              
[25] gtable_0.3.4                gtools_3.9.5               
[27] gridExtra_2.3               ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3                 jsonlite_1.8.8                    
  [3] farver_2.1.1                       fs_1.6.3                          
  [5] BiocIO_1.12.0                      zlibbioc_1.48.0                   
  [7] vctrs_0.6.5                        Rsamtools_2.18.0                  
  [9] memoise_2.0.1                      Cairo_1.6-2                       
 [11] RCurl_1.98-1.14                    ggtree_3.8.2                      
 [13] S4Arrays_1.2.0                     Rhdf5lib_1.22.1                   
 [15] SparseArray_1.2.3                  gridGraphics_0.5-1                
 [17] cachem_1.0.8                       GenomicAlignments_1.38.2          
 [19] igraph_1.5.1                       lifecycle_1.0.4                   
 [21] pkgconfig_2.0.3                    R6_2.5.1                          
 [23] fastmap_1.1.1                      gson_0.1.0                        
 [25] GenomeInfoDbData_1.2.11            digest_0.6.33                     
 [27] aplot_0.2.2                        colorspace_2.1-0                  
 [29] patchwork_1.1.3                    RSQLite_2.3.5                     
 [31] fansi_1.0.6                        httr_1.4.7                        
 [33] polyclip_1.10-6                    abind_1.4-5                       
 [35] compiler_4.3.2                     bit64_4.0.5                       
 [37] withr_3.0.0                        downloader_0.4                    
 [39] BiocParallel_1.36.0                viridis_0.6.4                     
 [41] DBI_1.2.1                          ggforce_0.4.1                     
 [43] MASS_7.3-60                        DelayedArray_0.28.0               
 [45] rjson_0.2.21                       HDO.db_0.99.1                     
 [47] tools_4.3.2                        ape_5.7-1                         
 [49] scatterpie_0.2.1                   glue_1.7.0                        
 [51] restfulr_0.0.15                    nlme_3.1-164                      
 [53] GOSemSim_2.26.1                    rhdf5filters_1.12.1               
 [55] shadowtext_0.1.2                   reshape2_1.4.4                    
 [57] fgsea_1.26.0                       generics_0.1.3                    
 [59] BSgenome_1.70.1                    tidyr_1.3.1                       
 [61] tidygraph_1.2.3                    utf8_1.2.4                        
 [63] XVector_0.42.0                     ggrepel_0.9.4                     
 [65] pillar_1.9.0                       yulab.utils_0.1.0                 
 [67] splines_4.3.2                      dplyr_1.1.4                       
 [69] tweenr_2.0.2                       treeio_1.24.3                     
 [71] lattice_0.22-5                     rtracklayer_1.62.0                
 [73] bit_4.0.5                          tidyselect_1.2.0                  
 [75] GO.db_3.18.0                       Biostrings_2.70.1                 
 [77] graphlayouts_1.0.2                 stringi_1.8.3                     
 [79] yaml_2.3.8                         lazyeval_0.2.2                    
 [81] ggfun_0.1.3                        codetools_0.2-19                  
 [83] BSgenome.Mmusculus.UCSC.mm10_1.4.3 ggraph_2.1.0                      
 [85] tibble_3.2.1                       qvalue_2.32.0                     
 [87] ggplotify_0.1.2                    cli_3.6.2                         
 [89] munsell_0.5.0                      png_0.1-8                         
 [91] XML_3.99-0.16.1                    blob_1.2.4                        
 [93] DOSE_3.26.2                        bitops_1.0-7                      
 [95] viridisLite_0.4.2                  tidytree_0.4.5                    
 [97] scales_1.3.0                       purrr_1.0.2                       
 [99] crayon_1.5.2                       rlang_1.1.3                       
[101] cowplot_1.1.2                      fastmatch_1.1-4                   
[103] KEGGREST_1.42.0   
