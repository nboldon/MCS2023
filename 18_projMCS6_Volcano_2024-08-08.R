#Setup an interactive session
salloc --account=eon -t 0-06:00 --mem=64G --nodes=1 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

############################################
############################################

############################################
############################################

projMCS6 <- addHarmony(
  ArchRProj = projMCS6,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
projMCS6 <- addImputeWeights(projMCS6)
getImputeWeights(projMCS6)

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

## Subset by cluster

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
############################################

markerGenesC8 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC8

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC8, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C8-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC8 <- getMarkerFeatures(
  ArchRProj = C8ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC8

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC8, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C8-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################
############################################
############################################

markerGenesC21 <- getMarkerFeatures(
  ArchRProj = C21ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC21, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C21-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC21 <- getMarkerFeatures(
  ArchRProj = C21ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC21, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C21-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################
############################################
############################################

markerGenesC22 <- getMarkerFeatures(
  ArchRProj = C22ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t1",
  bgdGroups = "t3"
)

markerGenesC21

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC22, name = "t1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C22-T1vT3-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

############################################

markerGenesC22 <- getMarkerFeatures(
  ArchRProj = C22ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  useGroups = "t2",
  bgdGroups = "t4"
)

markerGenesC22

###############

# Volcano Plots

pv <- plotMarkers(seMarker = markerGenesC22, name = "t2", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
pv

plotPDF(pv, name = "C22-T2vT4-Volcano_2024-08-08", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)


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
 [3] harmony_1.2.0               pheatmap_1.0.12            
 [5] enrichplot_1.20.3           clusterProfiler_4.8.3      
 [7] org.Mm.eg.db_3.17.0         AnnotationDbi_1.64.1       
 [9] rhdf5_2.44.0                SummarizedExperiment_1.32.0
[11] Biobase_2.62.0              MatrixGenerics_1.14.0      
[13] Rcpp_1.0.12                 Matrix_1.6-4               
[15] GenomicRanges_1.54.1        GenomeInfoDb_1.38.5        
[17] IRanges_2.36.0              S4Vectors_0.40.2           
[19] BiocGenerics_0.48.1         matrixStats_1.2.0          
[21] data.table_1.14.10          stringr_1.5.1              
[23] plyr_1.8.9                  magrittr_2.0.3             
[25] ggplot2_3.4.4               gtable_0.3.4               
[27] gtools_3.9.5                gridExtra_2.3              
[29] ArchR_1.0.2                

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
 [31] labeling_0.4.3                     fansi_1.0.6                       
 [33] httr_1.4.7                         polyclip_1.10-6                   
 [35] abind_1.4-5                        compiler_4.3.2                    
 [37] bit64_4.0.5                        withr_3.0.0                       
 [39] downloader_0.4                     BiocParallel_1.36.0               
 [41] viridis_0.6.4                      DBI_1.2.1                         
 [43] ggforce_0.4.1                      MASS_7.3-60                       
 [45] DelayedArray_0.28.0                rjson_0.2.21                      
 [47] HDO.db_0.99.1                      tools_4.3.2                       
 [49] ape_5.7-1                          scatterpie_0.2.1                  
 [51] glue_1.7.0                         restfulr_0.0.15                   
 [53] nlme_3.1-164                       GOSemSim_2.26.1                   
 [55] rhdf5filters_1.12.1                shadowtext_0.1.2                  
 [57] reshape2_1.4.4                     fgsea_1.26.0                      
 [59] generics_0.1.3                     BSgenome_1.70.1                   
 [61] tidyr_1.3.1                        tidygraph_1.2.3                   
 [63] utf8_1.2.4                         XVector_0.42.0                    
 [65] ggrepel_0.9.4                      pillar_1.9.0                      
 [67] yulab.utils_0.1.0                  splines_4.3.2                     
 [69] dplyr_1.1.4                        tweenr_2.0.2                      
 [71] treeio_1.24.3                      lattice_0.22-5                    
 [73] rtracklayer_1.62.0                 bit_4.0.5                         
 [75] tidyselect_1.2.0                   GO.db_3.18.0                      
 [77] Biostrings_2.70.1                  RhpcBLASctl_0.23-42               
 [79] graphlayouts_1.0.2                 stringi_1.8.3                     
 [81] yaml_2.3.8                         lazyeval_0.2.2                    
 [83] ggfun_0.1.3                        codetools_0.2-19                  
 [85] BSgenome.Mmusculus.UCSC.mm10_1.4.3 ggraph_2.1.0                      
 [87] tibble_3.2.1                       qvalue_2.32.0                     
 [89] ggplotify_0.1.2                    cli_3.6.2                         
 [91] munsell_0.5.0                      png_0.1-8                         
 [93] XML_3.99-0.16.1                    blob_1.2.4                        
 [95] DOSE_3.26.2                        bitops_1.0-7                      
 [97] viridisLite_0.4.2                  tidytree_0.4.5                    
 [99] scales_1.3.0                       purrr_1.0.2                       
[101] crayon_1.5.2                       rlang_1.1.3                       
[103] cowplot_1.1.2                      fastmatch_1.1-4                   
[105] KEGGREST_1.42.0                   

