## MCS2023


# Libraries.R
Primary libraries:
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



## The following code scripts are used to complete the project:


# bedConvert.py

Custom code from Qi Sun at Cornell University

Used to convert custom scATAC-seq library prep fragments to workable Cellranger files


# cellrangerRun.sh

Slurm script to create Cellranger files


# Rename_FragFiles.R

Renames Cellranger files; can then be moved to new folder for downstream analysis


# projMCS1_ArrowFiles.R

Creates arrow files in ArchR for downstream analysis


# projMCS1.R

Add doublet scores

Create df for nFrags & TSS Enrichment

Plot QC scores (Unique frags vs TSS Enrichment, Ridge & Violin plots for unique frags & TSS Enrichment)


# projMCS2.R

Filter doublets

Add Iterative LSI

Add Harmony

Add Clusters

Create cluster confusion matrix



# projMCS2_scEmbeddings.R

UMAP using LSI by sample & cluster

TSNE using LSI by sample & cluster

UMAP using Harmony by sample & cluster

TSNE using Harmony by sample & cluster


# projMCS2_MarkerGenes.R

getMarkerFeatures

Marker list by cluster (FDR & Log2FC)

Heatmaps for marker features; cowplots for all genes


# projMCS2_ArchRBrowser.R

Track plotting with ArchRBrowser


# projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


# projMCS4.R

Add reproduceable peak set using MACS2


# projMCS5.R

Add peak matrix

Identify marker peaks by cluster, etc.

Plotting marker peaks

Add motif annotations


# projMCS5_cellAbundance.R

Cell abundance analysis


# projMCS5_ClusterID_UMAP.R

Cluster identification using specific gene markers


# projMCS5_ArchRBrowser.R

Additional browser track regions of interest


# projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R








Code is adapted from: Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)
