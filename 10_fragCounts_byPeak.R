#Setup an interactive session
salloc --account=eon -t 0-04:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
library(Seurat)
library(Signac)
library(BiocManager)
library(BiocGenerics)
library(tidyr)
library(cowplot)
library(RColorBrewer)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)


# Load the ArchR project
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

# Check the project details and available matrices
projMCS7
getAvailableMatrices(projMCS7)

# Set up treatment groups by renaming sample IDs
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

# Assign the modified treatment names back to the project
projMCS7$treatment <- treatment

# Verify if the treatment column has been updated correctly
head(projMCS7$treatment)

# Identify marker peaks using the Wilcoxon test
markersPeaks <- getMarkerFeatures(
  ArchRProj = projMCS7,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

# Save marker peak results based on FDR and Log2FC cutoffs
markersPeaks_List <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")

# Export marker peaks for each cluster to a CSV file
for (i in seq_along(markersPeaks_List)) {
  clusterName <- names(markersPeaks_List)[i]
  write.csv(as.data.frame(markersPeaks_List[[i]]), 
            file = paste0("markersPeaks_List_", clusterName, "_2024-09-24.csv"), 
            row.names = FALSE)
}

# Extract group-wise fragment counts for peaks
fragmentCounts <- getGroupSE(
  ArchRProj = projMCS7,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters"
)

# Save the fragment counts (grouped by clusters) to a file
saveRDS(fragmentCounts, file = "fragmentCounts_PeakMatrix_2024-09-24.rds")


# Convert the SummarizedExperiment object to a data frame
# Extract the assay data (fragment counts) and rowData (peak information)
fragmentCounts_df <- as.data.frame(assay(fragmentCounts))  # Extracts the fragment count matrix
peakInfo <- rowData(fragmentCounts)  # Extracts the peak information

# Combine the peak information with the fragment counts
combined_df <- cbind(peakInfo, fragmentCounts_df)

# Save the combined data as a CSV file
write.csv(combined_df, file = "fragmentCounts_PeakMatrix_2024-09-24.csv", row.names = FALSE)


################################################
################################################
################################################
################################################
################################################

# The above code prints out the number of frags by cluster;
# the below code prints out the number of frags by treatment and cluster:


# Load the ArchR project
projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

# Check the project details and available matrices
projMCS7
getAvailableMatrices(projMCS7)

# Set up treatment groups by renaming sample IDs
treatment <- projMCS7$Sample

# Rename treatment groups based on sample IDs
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

# Assign the modified treatment names back to the project
projMCS7$treatment <- treatment

# Verify if the treatment column has been updated correctly
head(projMCS7$treatment)

# Combine clusters and treatments for grouping
projMCS7$Clusters_Treatment <- paste0(projMCS7$Clusters, "_", projMCS7$treatment)

# Extract group-wise fragment counts for peaks, grouped by clusters and treatment
fragmentCounts <- getGroupSE(
  ArchRProj = projMCS7,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters_Treatment"  # Group by both clusters and treatment
)

# Convert the SummarizedExperiment object to a data frame
fragmentCounts_df <- as.data.frame(assay(fragmentCounts))  # Extract fragment counts
peakInfo <- as.data.frame(rowData(fragmentCounts))  # Extract peak information

# Combine the peak information with fragment counts
combined_df <- cbind(peakInfo, fragmentCounts_df)

# Save the combined data as a CSV file
write.csv(combined_df, file = "fragmentCounts_Clusters_Treatment_2024-09-24.csv", row.names = FALSE)


################################################
################################################


## Plot the peak frag counts


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Check the column names first (for reference)
colnames(combined_df)

# Reshape the data to long format, ensuring to keep peak location information
long_df <- combined_df %>%
  pivot_longer(cols = starts_with("C"),  # Selecting columns that start with "C"
               names_to = "Treatment",
               values_to = "FragmentCount")

# Create a new data frame that includes peak information along with the long format data
long_df <- long_df %>%
  mutate(PeakLocation = paste(seqnames, start, end, sep = ":"))  # Combine seqnames, start, and end

# Create a distribution plot for each cluster and save each plot
clusters <- unique(long_df$Treatment)  # Treatment now serves as the cluster identifier

for (cluster in clusters) {
  # Create the plot for the current cluster
  plot <- ggplot(long_df %>% filter(Treatment == cluster), aes(x = PeakLocation, y = FragmentCount, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge") +  # Use bar plot for fragment counts
    labs(title = paste("Peak Fragment Counts for Cluster:", cluster),
         x = "Peak Location",
         y = "Fragment Count") +
    theme_minimal() +
    theme(legend.position = "top") +
    coord_flip()  # Flip coordinates for better readability
  
  # Save the plot
  ggsave(filename = paste0("Density_RankPlot_Cluster_", cluster, "_2024-09-26.png"),
         plot = plot,
         width = 10, height = 6, dpi = 300)  # Adjust dimensions as needed
}



################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################


## The following code rank orders the above plot by fragment counts in each cluster and by treatment group
# ERROR: This code doesn't actually plot the frags in rank order

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Check the column names first (for reference)
colnames(combined_df)

# Reshape the data to long format, ensuring to keep peak location information
long_df <- combined_df %>%
  pivot_longer(cols = starts_with("C"),  # Selecting columns that start with "C" (treatment columns)
               names_to = "Treatment",
               values_to = "FragmentCount")

# Create a new data frame that includes peak information along with the long format data
long_df <- long_df %>%
  mutate(PeakLocation = paste(seqnames, start, end, sep = ":"))  # Combine seqnames, start, and end

# Sort peaks by FragmentCount for each Treatment group
long_df <- long_df %>%
  group_by(Treatment) %>%
  arrange(desc(FragmentCount), .by_group = TRUE)  # Rank by FragmentCount in descending order

# Create a distribution plot for each cluster and save each plot
clusters <- unique(long_df$Treatment)  # Treatment now serves as the cluster identifier

for (cluster in clusters) {
  # Filter and reorder by fragment counts within the current cluster
  cluster_data <- long_df %>%
    filter(Treatment == cluster) %>%
    arrange(desc(FragmentCount))  # Ensure it's sorted for each plot
  
  # Create the plot for the current cluster
  plot <- ggplot(cluster_data, aes(x = reorder(PeakLocation, -FragmentCount), y = FragmentCount, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge") +  # Use bar plot for fragment counts
    labs(title = paste("Peak Fragment Counts for Cluster:", cluster),
         x = "Peak Location (Ordered by Fragment Count)",
         y = "Fragment Count") +
    theme_minimal() +
    theme(legend.position = "top") +
    coord_flip()  # Flip coordinates for better readability
  
  # Save the plot
  ggsave(filename = paste0("Density_Plot_Cluster_", cluster, "_Ranked_2024-09-25.png"),
         plot = plot,
         width = 10, height = 6, dpi = 300)  # Adjust dimensions as needed
}








################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################


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
 [1] nabor_0.5.0                 RColorBrewer_1.1-3         
 [3] cowplot_1.1.2               tidyr_1.3.1                
 [5] BiocManager_1.30.22         Signac_1.12.0              
 [7] Seurat_5.0.1                SeuratObject_5.0.1         
 [9] sp_2.1-2                    pheatmap_1.0.12            
[11] enrichplot_1.20.3           clusterProfiler_4.8.3      
[13] org.Mm.eg.db_3.17.0         AnnotationDbi_1.64.1       
[15] rhdf5_2.44.0                SummarizedExperiment_1.32.0
[17] Biobase_2.62.0              MatrixGenerics_1.14.0      
[19] Rcpp_1.0.12                 Matrix_1.6-4               
[21] GenomicRanges_1.54.1        GenomeInfoDb_1.38.5        
[23] IRanges_2.36.0              S4Vectors_0.40.2           
[25] BiocGenerics_0.48.1         matrixStats_1.2.0          
[27] data.table_1.14.10          stringr_1.5.1              
[29] plyr_1.8.9                  magrittr_2.0.3             
[31] ggplot2_3.4.4               gtable_0.3.4               
[33] gtools_3.9.5                gridExtra_2.3              
[35] ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.21                   splines_4.3.2                     
  [3] later_1.3.2                        BiocIO_1.12.0                     
  [5] bitops_1.0-7                       ggplotify_0.1.2                   
  [7] tibble_3.2.1                       polyclip_1.10-6                   
  [9] XML_3.99-0.16.1                    fastDummies_1.7.3                 
 [11] lifecycle_1.0.4                    globals_0.16.2                    
 [13] lattice_0.22-5                     MASS_7.3-60                       
 [15] plotly_4.10.3                      yaml_2.3.8                        
 [17] httpuv_1.6.13                      sctransform_0.4.1                 
 [19] spam_2.10-0                        spatstat.sparse_3.0-3             
 [21] reticulate_1.34.0                  pbapply_1.7-2                     
 [23] DBI_1.2.1                          abind_1.4-5                       
 [25] zlibbioc_1.48.0                    Rtsne_0.17                        
 [27] purrr_1.0.2                        ggraph_2.1.0                      
 [29] RCurl_1.98-1.14                    yulab.utils_0.1.0                 
 [31] tweenr_2.0.2                       GenomeInfoDbData_1.2.11           
 [33] ggrepel_0.9.4                      irlba_2.3.5.1                     
 [35] spatstat.utils_3.0-4               listenv_0.9.0                     
 [37] tidytree_0.4.5                     goftest_1.2-3                     
 [39] RSpectra_0.16-1                    spatstat.random_3.2-1             
 [41] fitdistrplus_1.1-11                parallelly_1.36.0                 
 [43] RcppRoll_0.3.0                     leiden_0.4.3.1                    
 [45] codetools_0.2-19                   DelayedArray_0.28.0               
 [47] DOSE_3.26.2                        ggforce_0.4.1                     
 [49] tidyselect_1.2.0                   aplot_0.2.2                       
 [51] farver_2.1.1                       viridis_0.6.4                     
 [53] spatstat.explore_3.2-5             GenomicAlignments_1.38.2          
 [55] jsonlite_1.8.8                     ellipsis_0.3.2                    
 [57] tidygraph_1.2.3                    progressr_0.14.0                  
 [59] ggridges_0.5.4                     survival_3.5-7                    
 [61] tools_4.3.2                        treeio_1.24.3                     
 [63] ica_1.0-3                          glue_1.7.0                        
 [65] SparseArray_1.2.3                  qvalue_2.32.0                     
 [67] dplyr_1.1.4                        withr_3.0.0                       
 [69] fastmap_1.1.1                      rhdf5filters_1.12.1               
 [71] fansi_1.0.6                        digest_0.6.33                     
 [73] R6_2.5.1                           mime_0.12                         
 [75] gridGraphics_0.5-1                 colorspace_2.1-0                  
 [77] scattermore_1.2                    GO.db_3.18.0                      
 [79] Cairo_1.6-2                        tensor_1.5                        
 [81] spatstat.data_3.0-3                RSQLite_2.3.5                     
 [83] utf8_1.2.4                         generics_0.1.3                    
 [85] rtracklayer_1.62.0                 htmlwidgets_1.6.4                 
 [87] graphlayouts_1.0.2                 httr_1.4.7                        
 [89] S4Arrays_1.2.0                     scatterpie_0.2.1                  
 [91] uwot_0.1.16                        pkgconfig_2.0.3                   
 [93] blob_1.2.4                         lmtest_0.9-40                     
 [95] XVector_0.42.0                     shadowtext_0.1.2                  
 [97] htmltools_0.5.7                    dotCall64_1.1-1                   
 [99] fgsea_1.26.0                       scales_1.3.0                      
[101] png_0.1-8                          ggfun_0.1.3                       
[103] rjson_0.2.21                       reshape2_1.4.4                    
[105] nlme_3.1-164                       cachem_1.0.8                      
[107] zoo_1.8-12                         KernSmooth_2.23-22                
[109] miniUI_0.1.1.1                     HDO.db_0.99.1                     
[111] restfulr_0.0.15                    pillar_1.9.0                      
[113] vctrs_0.6.5                        RANN_2.6.1                        
[115] promises_1.2.1                     xtable_1.8-4                      
[117] cluster_2.1.4                      Rsamtools_2.18.0                  
[119] cli_3.6.2                          compiler_4.3.2                    
[121] rlang_1.1.3                        crayon_1.5.2                      
[123] future.apply_1.11.0                fs_1.6.3                          
[125] stringi_1.8.3                      deldir_1.0-9                      
[127] viridisLite_0.4.2                  BiocParallel_1.36.0               
[129] munsell_0.5.0                      Biostrings_2.70.1                 
[131] lazyeval_0.2.2                     spatstat.geom_3.2-7               
[133] GOSemSim_2.26.1                    BSgenome_1.70.1                   
[135] RcppHNSW_0.5.0                     patchwork_1.1.3                   
[137] bit64_4.0.5                        future_1.33.0                     
[139] Rhdf5lib_1.22.1                    KEGGREST_1.42.0                   
[141] shiny_1.8.0                        ROCR_1.0-11                       
[143] igraph_1.5.1                       memoise_2.0.1                     
[145] ggtree_3.8.2                       fastmatch_1.1-4                   
[147] bit_4.0.5                          downloader_0.4                    
[149] BSgenome.Mmusculus.UCSC.mm10_1.4.3 ape_5.7-1                         
[151] gson_0.1.0  