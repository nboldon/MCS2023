#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=64G --nodes=1 --ntasks-per-node=16

#Load required dependencies
module load miniconda3/4.12.0
conda activate /pfs/tc1/project/eon/archr_env

#Load libraries
R

library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtree)
library(aplot)
library(patchwork)
library(pheatmap)
library(Seurat)
library(Signac)
library(BiocManager)
library(BiocGenerics)
library(ggplot2)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/Analysis")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS2 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Analysis/Save-ProjMCS2", force = FALSE, showLogo = FALSE)

#Transferring files to local computer
scp nboldon@beartooth.uwyo.edu:/project/eon/nboldon/MCS2023/Analysis/AwesomeFile.pdf ~/Desktop/MCS2023/Analysis

######################################################

#Session Info 10/26/23 and 10/27/23

R version 4.2.1 (2022-06-23)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /pfs/tc1/project/eon/archr_env/lib/libopenblasp-r0.3.21.so

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

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] Rtsne_0.16                         igraph_1.5.0
 [3] scran_1.24.1                       scuttle_1.6.3
 [5] SingleCellExperiment_1.18.1        harmony_0.1.1
 [7] uwot_0.1.16                        nabor_0.5.1
 [9] SeuratObject_4.1.3                 Seurat_4.3.0.1
[11] BSgenome.Mmusculus.UCSC.mm10_1.4.3 BSgenome_1.64.0
[13] rtracklayer_1.56.1                 Biostrings_2.64.1
[15] XVector_0.36.0                     pheatmap_1.0.12
[17] enrichplot_1.16.2                  clusterProfiler_4.4.4
[19] org.Mm.eg.db_3.15.0                AnnotationDbi_1.58.0
[21] rhdf5_2.40.0                       SummarizedExperiment_1.26.1
[23] Biobase_2.56.0                     MatrixGenerics_1.8.1
[25] Rcpp_1.0.10                        Matrix_1.5-4.1
[27] GenomicRanges_1.48.0               GenomeInfoDb_1.32.4
[29] IRanges_2.30.1                     S4Vectors_0.34.0
[31] BiocGenerics_0.42.0                matrixStats_1.0.0
[33] data.table_1.14.8                  stringr_1.5.0
[35] plyr_1.8.8                         magrittr_2.0.3
[37] ggplot2_3.4.2                      gtable_0.3.3
[39] gtools_3.9.4                       gridExtra_2.3
[41] ArchR_1.0.2

loaded via a namespace (and not attached):
  [1] utf8_1.2.3                spatstat.explore_3.2-1
  [3] reticulate_1.30           tidyselect_1.2.0
  [5] RSQLite_2.3.1             htmlwidgets_1.6.2
  [7] BiocParallel_1.30.4       scatterpie_0.2.1
  [9] ScaledMatrix_1.4.1        munsell_0.5.0
 [11] codetools_0.2-19          ica_1.0-3
 [13] statmod_1.5.0             future_1.32.0
 [15] miniUI_0.1.1.1            withr_2.5.0
 [17] spatstat.random_3.1-5     colorspace_2.1-0
 [19] GOSemSim_2.22.0           progressr_0.13.0
 [21] ROCR_1.0-11               tensor_1.5
 [23] DOSE_3.22.1               listenv_0.9.0
 [25] labeling_0.4.2            GenomeInfoDbData_1.2.8
 [27] polyclip_1.10-0           bit64_4.0.5
 [29] farver_2.1.1              downloader_0.4
 [31] parallelly_1.36.0         vctrs_0.6.3
 [33] treeio_1.20.2             generics_0.1.3
 [35] R6_2.5.1                  graphlayouts_1.0.0
 [37] rsvd_1.0.5                locfit_1.5-9.8
 [39] spatstat.utils_3.0-3      bitops_1.0-7
 [41] rhdf5filters_1.8.0        cachem_1.0.8
 [43] fgsea_1.22.0              gridGraphics_0.5-1
 [45] DelayedArray_0.22.0       promises_1.2.0.1
 [47] BiocIO_1.6.0              scales_1.2.1
 [49] ggraph_2.1.0              beachmat_2.12.0
 [51] Cairo_1.6-0               globals_0.16.2
 [53] goftest_1.2-3             tidygraph_1.2.3
 [55] rlang_1.1.1               splines_4.2.1
 [57] lazyeval_0.2.2            spatstat.geom_3.2-1
 [59] abind_1.4-5               yaml_2.3.7
 [61] reshape2_1.4.4            httpuv_1.6.11
 [63] qvalue_2.28.0             tools_4.2.1
 [65] ggplotify_0.1.1           ellipsis_0.3.2
 [67] RColorBrewer_1.1-3        ggridges_0.5.4
 [69] sparseMatrixStats_1.8.0   zlibbioc_1.42.0
 [71] purrr_1.0.1               RCurl_1.98-1.12
 [73] deldir_1.0-9              pbapply_1.7-2
 [75] viridis_0.6.3             cowplot_1.1.1
 [77] zoo_1.8-12                ggrepel_0.9.3
 [79] cluster_2.1.4             scattermore_1.2
 [81] DO.db_2.9                 lmtest_0.9-40
 [83] RANN_2.6.1                fitdistrplus_1.1-11
[85] patchwork_1.1.2           mime_0.12
 [87] xtable_1.8-4              XML_3.99-0.10
 [89] compiler_4.2.1            tibble_3.2.1
 [91] KernSmooth_2.23-21        crayon_1.5.2
 [93] shadowtext_0.1.2          htmltools_0.5.5
 [95] ggfun_0.1.1               later_1.3.1
 [97] tidyr_1.3.0               aplot_0.1.10
 [99] DBI_1.1.3                 tweenr_2.0.2
[101] MASS_7.3-60               cli_3.6.0
[103] metapod_1.4.0             pkgconfig_2.0.3
[105] GenomicAlignments_1.32.1  sp_2.0-0
[107] spatstat.sparse_3.0-2     plotly_4.10.2
[109] ggtree_3.4.2              dqrng_0.3.0
[111] yulab.utils_0.0.6         digest_0.6.32
[113] sctransform_0.3.5         RcppAnnoy_0.0.20
[115] spatstat.data_3.0-1       leiden_0.4.3
[117] fastmatch_1.1-3           tidytree_0.4.2
[119] edgeR_3.38.4              DelayedMatrixStats_1.18.2
[121] restfulr_0.0.15           shiny_1.7.4
[123] Rsamtools_2.12.0          rjson_0.2.21
[125] lifecycle_1.0.3           nlme_3.1-162
[127] jsonlite_1.8.7            Rhdf5lib_1.18.2
[129] BiocNeighbors_1.14.0      limma_3.52.4
[131] viridisLite_0.4.2         fansi_1.0.4
[133] pillar_1.9.0              lattice_0.21-8
[135] KEGGREST_1.36.3           fastmap_1.1.1
[137] httr_1.4.6                survival_3.5-5
[139] GO.db_3.15.0              glue_1.6.2
[141] png_0.1-8                 bluster_1.6.0
[143] bit_4.0.5                 ggforce_0.4.1
[145] stringi_1.7.3             blob_1.2.4
[147] BiocSingular_1.12.0       memoise_2.0.1
[149] dplyr_1.1.2               irlba_2.3.5.1
[151] future.apply_1.11.0       ape_5.7-1
