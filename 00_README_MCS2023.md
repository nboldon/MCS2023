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

# cellrangerRun.sh
Slurm script to create Cellranger files

# Rename_FragFiles.R
Renames Cellranger files; can then be moved to new folder for downstream analysis


##############################################
##############################################
##############################################

## 03. Project MCS1
- Creates projMCS1 for downstream analysis


## Code used for projMCS1:

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


######################################################################################
######################################################################################
######################################################################################


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
UMAP-Sample-Clusters.pdf
TSNE-Sample-Clusters.pdf
UMAP2Harmony-Sample-Clusters.pdf
TSNE2Harmony-Sample-Clusters.pdf

# projMCS2_MarkerGenes.R
- getMarkerFeatures
- Marker list by cluster (FDR & Log2FC)
- Heatmaps for marker features; cowplots for all genes

# projMCS2_ArchRBrowser.R
- Track plotting with ArchRBrowser


######################################################################################
######################################################################################
######################################################################################


## 05. Project MCS3

## projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


######################################################################################
######################################################################################
######################################################################################


## 06. Project MCS4 

## projMCS4.R

Add reproduceable peak set using MACS2


######################################################################################
######################################################################################
######################################################################################


## 07. Project MCS5

## projMCS5.R

Add peak matrix

Identify marker peaks by cluster, etc.

Plotting marker peaks

Add motif annotations


## projMCS5_cellAbundance.R

Cell abundance analysis


## projMCS5_ArchRBrowser.R

Additional browser track regions of interest


## Cluster identification using specific gene markers

projMCS5_ClusterID_UMAP.R

projMCS5_markerList_byCluster-Tx.R

projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R

projMCS5_geneMarkers_byClusterTxGrp_2024-04-10.R

projMCS5_Pairwise_byClusterTxGrp_2024-05-21.R

projMCS5_GeneMarkers_TxComp-byCluster_2024-04-29.R

projMCS5_Zu-2023_geneMarker_UMAPs.R

projMCS5_C1_Sorted_Heatmap.R


######################################################################################
######################################################################################
######################################################################################


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
