# MCS2023

# bedConvert.py
Custom code from Qi Sun at Cornell University
Used to convert custom scATAC-seq library prep fragments to workable Cellranger files

# Libraries.R
Primary libraries:
library(ArchR)
library(org.Mm.eg.db)
library(BiocManager)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)


# The following code scripts are used to setup the project:

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
Track plotting with ArchRBrowser

# projMCS3.R
Defining cluster ID with scRNA-seq 
Unable to complete; no RNA-seq data for cohort

projMCS4.R
# Add reproduceable peak set using MACS2

projMCS5.R
# Add peak matrix
# Identify marker peaks by cluster, etc.
# Plotting marker peaks
# Add motif annotations

projMCS5_cellAbundance.R





#Code is adapted from: Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)
