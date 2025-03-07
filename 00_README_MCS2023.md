# MCS2023

# Libraries.R
# Primary libraries used
library(ArchR)
library(org.Mm.eg.db)
library(BiocManager)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

# Conda_Update_Dec2023
Updated Beartooth conda environment for interactive session

# environment_2024-02-16.yml
Beartooth conda computing environment


####################################################
####################################################
####################################################


## The following code scripts are used to complete the project:

##################################################


## 02. Demultiplex code

# Custom code from Qi Sun at Cornell University
# Used to convert custom scATAC-seq library prep fragments to workable Cellranger files

README_demultiplex
sciatac.py
get_uniqList1.py
get_uniqList2.py
cratac_curated.txt.gz
cellranger)_commands.py
samplelist
bedConvert.py

cellrangerRun.sh
- Slurm script to create Cellranger files

Rename_FragFiles.R
- Renames Cellranger files; can then be moved to new folder for downstream analysis


##############################################
##############################################
##############################################

## 03. Project MCS1
- Creates projMCS1 for downstream analysis


# Code used for projMCS1:

# 03_projMCS1_ArrowFiles.R
- Creates arrow files in ArchR to create project

# 03_projMCS1.R
- Adds doublet scores
- Creates dataframes for nFrags and TSS Enrichment
- Plots QC scores 
  - Density plot for nFrags vs TSS Enrichment
  - Ridge and violin plots for nFrags and TSS Enrichment individually plotted)

# QC files:
- 03_TSS-vs-Frags.pdf
- 03_TSS7_QC-MCS1.pdf
- 03_QC_FragSize-Distro_2024-02-29.pdf
- 03_Peak-Call-Summary.pdf


#############################################
#############################################
#############################################


## 04. Project MCS2


# projMCS2.R
- Filters doublets
- Adds Iterative LSI
- Adds Harmony
- Adds Clusters
- Creates cluster confusion matrix


# projMCS2_scEmbeddings.R
- UMAP using LSI by sample & cluster
- TSNE using LSI by sample & cluster
- UMAP using Harmony by sample & cluster
- TSNE using Harmony by sample & cluster

# Files created from projMCS2_scEmbedding.R
04_UMAP-Sample-Clusters.pdf
04_TSNE-Sample-Clusters.pdf
04_UMAP2Harmony-Sample-Clusters.pdf
04_TSNE2Harmony-Sample-Clusters.pdf


# projMCS2_MarkerGenes.R
# Group by: Cluster (no subset), FDR <= 0.01, Log2FC >= 1.25
- getMarkerFeatures
- Marker list by cluster (FDR & Log2FC)
- Heatmaps for marker features; cowplots for all genes
# NOTE: projMCS2_MarkerGenes.R did not calculate abs(Log2FC)

# projMCS5_MarkerGenes.R 
# Group by: Cluster (no subset), FDR <= 0.01, abs(Log2FC) >= 1.25
- Files saved for each cluster; ex: C1_MarkerGenes_2025-03-06.csv

# Summary files created from analysis results of projMCS2_MarkerGenes.R
- 04_Cluster_Analysis_2024-02-12.xlsm
    - Each tab summarizes cell count information by cluster, including:
      - cell counts by sample and cluster (incl. median TSS & nFrags per cell)
      - normalized cell abundance percentages by sample (cells in cluster by sample divided by total cells in  sample; incl. total cell counts by sample)
      - normalized cell abundance by treatment group (incl. total cell counts by cluster and treatment group)
      - cell counts by cluster and genotype
      - cell counts by cluster and sequencing lane
      - cell counts by PCR plate (library preparation day)
- 04_GeneMarkers_byCluster_1-31-2024.xlsx
  - Includes cell type-specific marker gene notes from literature reviews in tabs 1 & 2
  - There is a tab for each cluster that highlights:
    - the top 50 marker genes, based on FDR 
    - cell type-specific marker genes found in literature reviews
  - There are two additional tabs for C1, C3, C8, C18, and C21 that highlight:
    - One tab specifies the Panther GO family, molecular function, biological process, cellular component, protein class, and pathways affected by specified genes. 
    - The other tab displays seperate bargraphs for Panther GO categories associated with molecular function, biological process, cellular component, protein class, and pathways affected by specified genes.
- 04_MCS2023_1-31-2024.pptx
    - Summarizes key points of interest from the above analysis

# Heatmaps of projMCS2 by cluster from projMCS2_MarkerGenes.R
- 04_UMAP-MarkerGenes-WO-Imputation.pdf
- 04_GeneScores-Marker-Heatmap.pdf
    - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- 04_GeneMarkers_Top50_Heatmap_1-25-2024.pdf


# projMCS2_ArchRBrowser.R
- Track plotting with ArchRBrowser
- 04_Plot-Tracks-MarkerGenes_2024-04-04.pdf 
    - Creates browser tracks for the following genes:
    - "Aqp4", "Aldh1l1", "Mlc1", "Cbs", "Ppp1r3c", "Plcd4", "Dio2", #Astrocyte
    - "Cldn11", "Cx3cr1", "Csf1r", "Sparc", "Trem2", "Ccl4", "Cd14", "Tyrobp", "C1qa", #Microglia
    - "Olig1", "Mbp", "Opalin", "Mag", "Mog", "Cldn11", "Ugt8a", "Olig2", #Oligodendrocyte
    - "Spock3", "Gad1", "Grin3a", "Adarb2", "Grik1", "Lhx6", "Pvalb", "Gad2",  #GABAergic
    - "Sulf1", "Slc17a8", "Tshz2", "Slc17a6", "Neurod6", #Glutamatergic
    - "Cldn5" #"CD31"DoesNotExist #Endothelial
      


##############################################
##############################################
##############################################


## 05. Project MCS3

## projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


###############################################
###############################################
###############################################


## 06. Project MCS4 

## projMCS4.R

Add pseudobulk replicates & reproduceable peak set using MACS2


###############################################
###############################################
###############################################


## 07. Project MCS5

# projMCS5.R
- Adds peak matrix
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- Plots marker peaks
  - cutOff = "FDR<=0.01 & Log2FC>=1.25"
- Adds motif annotations
  - motifSet = "cisbp"
## NOTE: The 07_projMCS5.R code did not use abs(Log2FC)
## The 07_projMCS5-w-ABS-Log2FC.R code accommodates for that and reruns getMarkerFeatures before applying cutoffs

# 07_projMCS5-w-ABS-Log2FC.R
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
  - 07_Peak_markerList_2025-03-06.csv
- Plots marker peaks
  - cutOff = "FDR<=0.01 & abs(Log2FC)>=1.25"
  - 07_Peak-Marker-Heatmap_2025-03-06.pdf


######################################


# 07_projMCS5_ArchRBrowser.R
- Additional browser track gene regions of interest
- Subset by treatment group
- C18 Subcluster subset by treatment group
- Files generated using this code:
  - 07_Browser-Tracks_C18-grpByTx_2024-06-05.pdf
  - 07_Browser-Tracks_T3-grpByClusters_2024-06-05.pdf


###############################################
###############################################
###############################################


## 08. Cell Abundance

# 08_projMCS5_cellAbundance.R
- Cell abundance analysis uses 07_Cell-Abundance_2024-02-08.csv
- Performs ANOVA to determine significant differences between treatment groups
# 08_Abundance_Boxplots_2024-04-11.R
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved stats from 07_projMCS5_cellAbundance.R
- Uses ANOVA, t-tests, and Bonforonni correction
- C8 & C20 show significant differences
# 08_Cell_Abundance_Boxplots_2024-12-09.R
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved color palette than those run on 2024-04-11
- C8 shows significant differences by treatment
# 08_CellAbundance_Boxplots_wStats_2025-02-21
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved color palette and stats than those run on 2024-04-11 & 2024-12-09
- Uses ANOVA, then for those that are sig, uses t-tests, followed by Bonforonni correction
- See below for sig differences by comparison

# Files created by Cell_Abundance_Boxplots_2024-12-09.R (do not display stats):
- 08_CellAbund_Boxplots_2024-12-09.jpg
- 08_CellAbund_Boxplots_wMean_2024-12-09.jpg
- 08_CellAbund_Boxplots_ByCellType_2024-12-09.jpg
- 08_CellAbund_Boxplots_by_Genotype-CellType_2024-12-09.jpg
- 08_CellAbund_Boxplot_TsCombined_2024-12-09.pdf
- 08_CellAbund_Boxplot_2NCombined_2024-12-09.pdf

# Files created by Cell_Abundance_Boxplots_wStats_2025-02-21.R (incl. sig stats):
- 08_CellAbund_Boxplots_Stats_wMean_byCluster_2025-02-21.jpg
- 08_CellAbund_Boxplots_Stats_byCluster_2025-02-21.jpg
- 08_CellAbund_Genotype_Stats_2025-02-21.csv
- 08_CellAbund_Boxplots_by_Genotype-CellType_Stats_2025-02-21.pdf
- 08_CellAbund_Boxplots_by_Genotype-CellType_Stats_2025-02-21.jpg
  * P-adj pairwise comparisons show sig diff in Astrocytes comparing 2N & 2N+ vs. Ts & Ts+
- 08_CellAbund_Pairwise_Stats_byCluster_2025-02-21.csv
  * P-adj pairwise comparisons show sig diff in cluster specific comparisons:
    - C8 2N+ vs Ts, C8 2N+ vs Ts+, C14 2N+ vs Ts+, C20 2N vs Ts
- 08_CellAbund_Boxplots_2N_vs_2N+_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_2N_vs_2N+_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_2N_vs_2N+_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
  * No sig differences found in 2N vs 2N+ by cell type
- 08_CellAbund_Boxplots_2N_vs_Ts_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_2N_vs_Ts_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_2N_vs_Ts_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
  * P-adj pairwise comparisons show sig diff in Astrocytes comparing 2N vs Ts
- 08_CellAbund_Boxplots_2N_vs_Ts+_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_2N_vs_Ts+_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_2N_vs_Ts+_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
  * P-adj pairwise comparisons show sig diff in Astrocytes comparing 2N vs Ts+
- 08_CellAbund_Boxplots_2N+_vs_Ts_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_2N+_vs_Ts_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_2N+_vs_Ts_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
  * P-adj pairwise comparisons show sig diff in Astrocytes comparing 2N+ vs Ts
- 08_CellAbund_Boxplots_2N+_vs_Ts+_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_2N+_vs_Ts+_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_2N+_vs_Ts+_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
  * P-adj pairwise comparisons show sig diff in Astrocytes & GlutOligo comparing 2N+ vs Ts
- 08_CellAbund_Boxplots_Ts_vs_Ts+_cellType_2025-02-21.pdf
- 08_CellAbund_Boxplots_Ts_vs_Ts+_cellType_2025-02-21_pairwise_cellType_stats_2025-02-21.csv
- 08_CellAbund_Boxplots_Ts_vs_Ts+_cellType_2025-02-21_ANOVA_cellType_stats_2025-02-21.csv
  * No sig differences found in 2N vs 2N+ by cell type


###############################################
###############################################
###############################################


## 09. Marker genes subset by Treatment, grouped by cluster


09_projMCS5_geneMarkers_byClusterTxGrp_2024-04-10.R
- Subsets by treatment, group by cluster
- cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"

09_projMCS5_geneMarkers_byClusterTxGrp_2024-04-26.R
- code updated from the above
- cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
- Subsets by treatment, group by cluster
- File names ex: T1_C1_GeneMarkers_2024-04-10.csv
- Subsets by cluster, group by treatment
- File names ex: C1_GeneMarker_List_2024-04-08.csv & GeneMarker_List_Master_2024-04-08.csv

# Marker genes that were subset by treatment group, then groupBy cluster
# were combined into one spreadsheet using: 
09_Combine_Spreadsheets.R
- File names ex: C1_TxMarkers_Combined.csv


###############################################
###############################################
###############################################


## 10. Marker genes subset by cluster, grouped by treatment


10_projMCS5_markerList_byCluster-Tx.R
- Creates heatmap & z-scores for 4 treatment groups using all sig genes found in C18, subset by C18
- Generates file: 07_C18_All127_byTx_zscores_2024-03-13.csv


10_projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R
- Subsets by cluster
- Uses marker genes known to be significant for the cluster of interest
- Creates z-scores and heatmaps by treatment
- File names ex: C4_byTx_zscores_2024-03-21.csv & C4_TxComp_Heatmap_2024-03-21.pdf
10_projMCS5_C1_Sorted_Heatmap.R
- First sorted heatmaps created
10_projMCS7_GeneMarkers_byCluster-TxGrp_2024-10-17.R
- Upgrades from previous 2 code records include sorted heatmaps in 1 file
10_projMCS5_GeneMarkers_byClusterTxComp_Loop_2024-10-17.R
- Subsets by cluster, group by treatment; loop for all clusters
10_projMCS7_GeneMarkers_byCellType-TxGrp_2025-01-16.R
- Similar code to projMCS7_GeneMarkers_byCluster-TxGrp_2025-01-16.R but focuses on cell type groups
* - Files saved ex: astro_GeneMarker_List_2025-01-16.csv


###############################################
###############################################
###############################################







###############################################
###############################################
###############################################

## 08. Project MCS6

## ProjMCS6 - Harmony Added; Single cell embeddings & tx comparisons

projMCS6_2024-06-20.R

projMCS6_GeneMarkers_TxComp-byCluster_2024-06-26.R

projMCS6_Volcano_2024.08-08.R

T3vT124_byCluster_2024-08-09.R

T1-2-4vT3_TxComp_FDR-0-1_Log2FC-0-5_2024-07-31.R

Q4_GenotypeTxComps_2024-06-28.R

Q5_GenotypeClusterComps_2024_06-28.R

Q6_TvTbyClusterComps_2024-06-28.R

Volcano_byCluster-Tx.R


######################################################################################
######################################################################################
######################################################################################


## 09. Project MCS7 

## ProjMCS7 - Peak Matrix and motif analysis

projMCS7.R

normExp_projMCS7.R

peakVSmatrix.R


######################################################################################
######################################################################################
######################################################################################


## 10. Project MCS9

## ProjMCS9 - Peak & fragment count analysis

projMCS9.R

BetterBrowserTracks.R

Peak_TxComps-byCluster_2024-09-09.R

Peak_heatmaps_2024-09-11.R

fragCounts_GenomeWide_TxComps.R

fragCounts_byPeak.R

geneFrags.R

geneFragsLoop.R

getMatrix_MCS1.R

nFrags-TSS_byTx-Cluster.R

normPeakStats_2024-09-08.R

peakFrags.R

peakFragsLoop.R

peak_fragCounts_byTx_2024-09-18.R

peaksVSmotifs.R


######################################################################################
######################################################################################
######################################################################################


## 11. Cluster 18 Subset

C18_Subset.R

C18_Subset_TxComps_byCluster.R

C18_Subset_GO_Barplot.R

C18_TxComps_byCluster.R

C18_T1-2-4v3_byCluster_2024-07-31.R


######################################################################################
######################################################################################
######################################################################################


## 12. GO Analysis

GO_Analysis_2024-08-01.R

GO_Barplot.R

GO_Barplot_2.R

GO_Plots_Loop.R


######################################################################################
######################################################################################
######################################################################################


## 13. DotPlots

Dotplot.R

revised_DotPlot.R





######################################################################################

Code is adapted from: Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)
