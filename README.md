# MCS2023


## Libraries.R
Primary libraries:
library(ArchR)
library(org.Mm.eg.db)
library(BiocManager)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

## Conda_Update_Dec2023
Updated Beartooth conda environment for interactive session

## environment_2024-02-16.yml
Beartooth conda computing environment



# The following code scripts are used to complete the project:


## Demultiplex code

Custom code from Qi Sun at Cornell University

Used to convert custom scATAC-seq library prep fragments to workable Cellranger files

README_demultiplex

sciatac.py

get_uniqList1.py

get_uniqList2.py

cratac_curated.txt.gz

cellranger)_commands.py

samplelist

bedConvert.py



## cellrangerRun.sh

Slurm script to create Cellranger files


## Rename_FragFiles.R

Renames Cellranger files; can then be moved to new folder for downstream analysis


## projMCS1_ArrowFiles.R

Creates arrow files in ArchR for downstream analysis


## projMCS1.R

Add doublet scores

Create df for nFrags & TSS Enrichment

Plot QC scores (Unique frags vs TSS Enrichment, Ridge & Violin plots for unique frags & TSS Enrichment)


## projMCS2.R

Filter doublets

Add Iterative LSI

Add Harmony

Add Clusters

Create cluster confusion matrix



## projMCS2_scEmbeddings.R

UMAP using LSI by sample & cluster

TSNE using LSI by sample & cluster

UMAP using Harmony by sample & cluster

TSNE using Harmony by sample & cluster


## projMCS2_MarkerGenes.R

getMarkerFeatures

Marker list by cluster (FDR & Log2FC)

Heatmaps for marker features; cowplots for all genes


## projMCS2_ArchRBrowser.R

Track plotting with ArchRBrowser


## projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


## projMCS4.R

Add reproduceable peak set using MACS2


## projMCS5.R

Add peak matrix

Identify marker peaks by cluster, etc.

Plotting marker peaks

Add motif annotations


## projMCS5_cellAbundance.R

Cell abundance analysis


## Cluster identification using specific gene markers

projMCS5_ClusterID_UMAP.R

projMCS5_geneMarkers_byClusterTxGrp_2024-04-10.R

projMCS5_Pairwise_byClusterTxGrp_2024-05-21.R

projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R

projMCS5_GeneMarkers_TxComp-byCluster_2024-04-29.R

projMCS5_Zu-2023_geneMarker_UMAPs.R


## projMCS5_ArchRBrowser.R

Additional browser track regions of interest


## projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R


## ProjMCS6 - Harmony Added; Single cell embeddings & tx comparisons

projMCS6_2024-06-20.R

projMCS6_GeneMarkers_TxComp-byCluster_2024-06-26.R

projMCS6_Volcano_2024.08-08.R

T3vT124_byCluster_2024-08-09.R

T1-2-4vT3_TxComp_FDR-0-1_Log2FC-0-5_2024-07-31.R

Q4_GenotypeTxComps_2024-06-28.R

Q5_GenotypeClusterComps_2024_06-28.R

Q6_TvTbyClusterComps_2024-06-28.R


## ProjMCS7 - Peak Matrix and motif analysis

projMCS7.R

normExp_projMCS7.R

peakVSmatrix.R








Code is adapted from: Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)
