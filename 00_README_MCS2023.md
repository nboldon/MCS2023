# Project TIN CANN
## By: Naomi Boldon

# MCS2023


##################################################
##################################################
##################################################


# 01. The following software and package versions were used to complete the project:

##################################################

## 01_LICENSE
Specifies permissions granted and copyright notice. 

## 01_Libraries.R
Catelog of libraries used for the project. 

## 01_Conda_Update_Dec2023
Updated Beartooth conda environment for interactive session

## 01_beartooth_environment_2024-Feb.yml
Beartooth conda computing environment

## 01_beartooth_environment_2024-Oct
Beartooth conda computing environment

## 01_local_environment_2025-Mar
Naomi's personal computer computing environment


##################################################
##################################################
##################################################


# The following code scripts are used to complete the project:

##################################################


# 02. Demultiplex code

- Used to convert custom scATAC-seq library prep fragments to workable Cellranger files

- Custom code from Dr. Qi Sun at Cornell University


The script sciatac.py does two things:
1. demultiplexing samples
2. converting the lab custom barcodes into Cellranger 16bp barcodes (based on the Cellranger-atac whitelist file). 
The cellranger-atac requires three fastq files:
"sampleName_S1_L001_R1_001.fastq.gz"
"sampleName_S1_L001_R2_001.fastq.gz"
"sampleName_S1_L001_R3_001.fastq.gz"
(R1 and R3 are paired end reads that can be aligned to the reference genome. R2 contains the 16bp cell barcodes.)


Section 1. run sciatac.py to convert raw fastq data files to Cellranger-atac format 
a. prepare the samplelist file
See example file samplelist in the directory. It is a tab delimited text file with three columns: sampleName, bc1, bc2
(Make sure there is no space or funky characters in the smaple names. If a sample uses multiple barcodes, use multiple lines for the sample.)

b. run sciatac.py on a server with plenty of temporary storage (~2TB for ~100 billion read pairs), and >=20 cpu cores. If you are using a computer with less cpu, modify the scripts, set both "parallelSamples" and "parallelChunks" to number of CPUs. 

sciatac.py -1 raw_fq_R1.gz -2 raw_R2.fq.gz -s samplelist  -o outputDirName

(make sure that the file cratac_curated.txt.gz is located at the same directory as the script file sciatac.py)

Section 2. run cellranger 
basic cellranger command
cellranger-atac count --id=WT-8h_run \
                        --reference=/workdir/qisun/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                        --fastqs=/workdir/qisun/out/WT-8h \
                        --sample=WT-8h \
                        --localcores=32 \
                        --localmem=128
 
There is a script that can prepare the batch command file, where "-i" input directory set to the directory from previous command. It would create a batch script run.sh
cellranger_commands.py -i outputDirName -r /workdir/qisun/refdata-cellranger-arc-mm10-2020-A-2.0.0 -o run.sh

To run it on a large mem gen2 server (96 cores, and 500gb ram)
parallel -j3 < run.sh


Appendix. curate the cellranger barcode whitelist
The curated barcode file "cratac_curated.txt.gz" must be kept in the same directory as the script. This file is from the cellranger software directory, located in cellranger-atac-2.1.0/lib/python/atac/barcodes

The file 737K-cratac-v1.txt.gz must be scrambled before use
zcat 737K-cratac-v1.txt.gz | shuf |gzip -c > newfile.txt.gz 

After scrambling, the curated list was prepared with two scripts: get_uniqList1.py and get_uniqList2.py.  get_uniqList1.py prepare about ~30,000 barcodes with hamming distance <3,  get_uniqList2.py would scan the rest of post-shuf barcodes, and add extra barcodes in case more cell barcodes are needed.. 


02_README_demultiplex
02_sciatac.py
02_get_uniqList1.py
02_get_uniqList2.py
02_cratac_curated.txt.gz
02_cellranger_commands.py
02_samplelist
02_bedConvert.py

02_cellrangerRun.sh
Slurm script to create Cellranger files

02_Rename_FragFiles.R
Renames Cellranger files; can then be moved to new folder for downstream analysis


##############################################
##############################################
##############################################


# 03. Project MCS1

- Creates arrow files and project MCS1 for quality control and downstream analysis.

## 03_projMCS1_ArrowFiles.R
Creates arrow files in ArchR to create project

## 03_projMCS1.R
Adds doublet scores
Creates dataframes for nFrags and TSS Enrichment
Plots QC scores 
  Density plot for nFrags vs TSS Enrichment
  Ridge and violin plots for nFrags and TSS Enrichment individually plotted)

## QC files:
03_TSS-vs-Frags.pdf
03_TSS7_QC-MCS1.pdf
03_QC_FragSize-Distro_2024-02-29.pdf
03_Peak-Call-Summary.pdf


#############################################
#############################################
#############################################


# 04. Project MCS2

- Filters doublets, runs initial DA gene accessibility groupBy cluster (no subset), and creates initial visualizations for all single cells in the project, including UMAPs, t-SNE, heatmaps, browser tracks. 

## 04_projMCS2.R
- Filters doublets
- Adds Iterative LSI
- Adds Harmony
- Adds Clusters
- Creates cluster confusion matrix


## 04_projMCS2_scEmbeddings.R
- UMAP using LSI by sample & cluster
- TSNE using LSI by sample & cluster
- UMAP using Harmony by sample & cluster
- TSNE using Harmony by sample & cluster

## Files created from projMCS2_scEmbedding.R
04_UMAP-Sample-Clusters.pdf
04_TSNE-Sample-Clusters.pdf
04_UMAP2Harmony-Sample-Clusters.pdf
04_TSNE2Harmony-Sample-Clusters.pdf


## 04_projMCS2_MarkerGenes.R
- Group by: Cluster (no subset), FDR <= 0.01, Log2FC >= 1.25
- getMarkerFeatures
- Marker list by cluster (FDR & Log2FC)
- Heatmaps for marker features; cowplots for all genes
# NOTE: projMCS2_MarkerGenes.R did not calculate abs(Log2FC)
04_Cluster_GeneMarkers_2025-03-06.R
- Uses 04_projMCS2_MarkerGenes.R code, but calculates abs(Log2FC)

## 04_projMCS5_MarkerGenes.R 
- Group by: Cluster (no subset), FDR <= 0.01, abs(Log2FC) >= 1.25
- Files saved for each cluster; ex: C1_MarkerGenes_2025-03-06.csv

## Summary files created from analysis results of projMCS2_MarkerGenes.R
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

## Heatmaps of projMCS2 by cluster from projMCS2_MarkerGenes.R
- 04_UMAP-MarkerGenes-WO-Imputation.pdf
- 04_GeneScores-Marker-Heatmap.pdf
    - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- 04_GeneMarkers_Top50_Heatmap_1-25-2024.pdf

## 04_HighlightCells_2024-11-22.R
- Highlights treatment group cells of interest on UMAPs.
    - plotEmbedding uses colorBy = cellColData, name = treatment
    - Plots specific colors by treatment group. 

## 04_projMCS2_ArchRBrowser.R
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


# 05. Project MCS3

## 05_projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


###############################################
###############################################
###############################################


# 06. Project MCS4 

## 06_projMCS4.R

Add pseudobulk replicates & reproduceable peak set using MACS2


###############################################
###############################################
###############################################


# 07. Project MCS5

- Adds peak matrix and analyzes marker genes and peaks by cluster. 


## 07_projMCS5.R
- Adds peak matrix
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- Plots marker peaks
  - cutOff = "FDR<=0.01 & Log2FC>=1.25"
- Adds motif annotations
  - motifSet = "cisbp"
NOTE: The 07_projMCS5.R code did not use abs(Log2FC)
The 07_projMCS5-w-ABS-Log2FC.R code accommodates for that and reruns getMarkerFeatures before applying cutoffs

## 07_projMCS5_Top50.R
- No subset, groupBy cluster, cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
    - File name ex: GeneMarker_List_Master.csv & C1_GeneMarker_List.csv
- Filters top genes and plots heatmap.
    - File: GeneMarkers_Top50_1-25-2024.pdf

## 07_projMCS5-w-ABS-Log2FC.R
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
  - 07_Peak_markerList_2025-03-06.csv
- Plots marker peaks
  - cutOff = "FDR<=0.01 & abs(Log2FC)>=1.25"
  - '07_Peak-Marker-Heatmap_2025-03-06.pdf'

## 07_DESeq2.R
- Code does not run. 


## 07_projMCS5_ArchRBrowser.R
- Additional browser track gene regions of interest
- Subset by treatment group
- C18 Subcluster subset by treatment group
- Files generated using this code:
  - 07_Browser-Tracks_C18-grpByTx_2024-06-05.pdf
  - 07_Browser-Tracks_T3-grpByClusters_2024-06-05.pdf
  
## 07_Seurat_DotPlot.R
- Code does not run. 


###############################################
###############################################
###############################################


# 08. Cell Abundance

## 08_projMCS5_cellAbundance.R
- Cell abundance analysis uses 07_Cell-Abundance_2024-02-08.csv
- Performs ANOVA to determine significant differences between treatment groups
## 08_Abundance_Boxplots_2024-04-11.R
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved stats from 07_projMCS5_cellAbundance.R
- Uses ANOVA, t-tests, and Bonforonni correction
- C8 & C20 show significant differences
## 08_Cell_Abundance_Boxplots_2024-12-09.R
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved color palette than those run on 2024-04-11
- C8 shows significant differences by treatment
## 08_CellAbundance_Boxplots_wStats_2025-02-21
- Cell abundance analysis uses 07_Cell-Abund_Boxplot_2024-04-11.csv
- Improved color palette and stats than those run on 2024-04-11 & 2024-12-09
- Uses ANOVA, then for those that are sig, uses t-tests, followed by Bonforonni correction
- See below for sig differences by comparison

## Files created by Cell_Abundance_Boxplots_2024-12-09.R (do not display stats):
- 08_CellAbund_Boxplots_2024-12-09.jpg
- 08_CellAbund_Boxplots_wMean_2024-12-09.jpg
- 08_CellAbund_Boxplots_ByCellType_2024-12-09.jpg
- 08_CellAbund_Boxplots_by_Genotype-CellType_2024-12-09.jpg
- 08_CellAbund_Boxplot_TsCombined_2024-12-09.pdf
- 08_CellAbund_Boxplot_2NCombined_2024-12-09.pdf

## Files created by Cell_Abundance_Boxplots_wStats_2025-02-21.R (incl. sig stats):
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


# 09. Marker genes subset by Treatment, grouped by cluster


## 09_projMCS5_geneMarkers_byClusterTxGrp_2024-04-10.R
- Subsets by treatment, group by cluster
- cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
## 09_projMCS5_geneMarkers_byClusterTxGrp_2024-04-26.R
- code updated from the above
- cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
- Subsets by treatment, group by cluster
- File names ex: T1_C1_GeneMarkers_2024-04-10.csv
- Subsets by cluster, group by treatment
- File names ex: C1_GeneMarker_List_2024-04-08.csv & GeneMarker_List_Master_2024-04-08.csv

## 09_Combine_Spreadsheets.R
- Combines marker genes subset by treatment group groupBy cluster into one summary spreadsheet. 
- File names ex: C1_TxMarkers_Combined.csv

## 09_Combined_GeneMarkers_TxComps-byCluster-OR-CellType_Heatmap_2024-11-07.R
- File load ex: C1_Combined_GeneMarkers_Comps_2024-11-07.csv
- Plots a heatmap.
    - File name ex: C1_Combined_GeneMarkers_Comps_Heatmap_2024-11-07.pdf
- Loops through all clusters to plot heatmaps.
    - Fie name ex: C%d_Combined_GeneMarkers_Comps_Heatmap_2024-11-07.pdf
- Loops through all clusters, removes NA values, and plots heatmaps.
    - File name ex: C%d_Combined_GeneMarkers_Comps_Heatmap_NA-removed_2024-11-07.pdf
- Loops through cell types, removes NA values, and plots heatmaps.
    - File name ex: %s_Combined_GeneMarkers_Comps_Heatmap_NA-removed_2024-12-02.pdf
    *NOTE: Files saved under 11. Marker genes by cell type
- Plots heatmaps sorted by treatment group
    - Adapted from projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R
    
## 09_Combined_GeneMarkers_TxComps-byCluster_Heatmap_2024-12-19.R
- Runs the same code as 09_Combined_GeneMarkers_TxComps-byCluster-OR-CellType_Heatmap_2024-11-07.R but is updated for Viridis color palette ("#440154FF", "white", "#1F9E89FF"). 
- Code specifies clusters that did not have data to plot. 
    - File ex: C%d_Combined_GeneMarkers_Comps_Heatmap_2024-12-19.pdf & C%d_Combined_GeneMarkers_Comps_Heatmap_NA-removed_2024-12-19.pdf
- Includes heatmaps with labeled known markers. 
    - File ex: C18_All127_TxComp_Heatmap_2024-03-13.pdf
- Prints z-scores from the heatmap.
    - File ex: C18_All127_byTx_zscores_2024-03-13.csv
    
## 09_projMCS5_SingleTxGrp_Analysis.R
- Subsets by treatment group for downstream analysis. 
    - File name ex: T1_GeneMarker_Heatmap.pdf
    
## 09_GeneCounts_byTx-Comp_Barplot_2024-11-08.R
- Creates a barplot comparing the DA genes obtained by treatment group pairwise comparisons of subset clusters, and DA genes obtained by subsetting treatment groups and groupBy Cluster. 
    - File load ex: C3_TxComp_GeneUnions_2024-07-24.xlsx
    - Runs using a rainbow color palette, and Viridis color palette. 
    - File: gene_counts_by_comparison.jpg
- Loop runs through all clusters, code specifies when data was not found for a cluster.
    - File ex: C%d_gene_counts_by_comps_2024-11-08.jpg
    
## 09_Sorted_Heatmaps_2024-11-10.R
- Load file ex: C1_byTx_zscores_2024-03-21.csv
- Sorts data based on mean of the rows, updates labels, plots heatmap.
    - File ex: C1_byTx_sorted_heatmap_2024-10-16.pdf
- Ordered by treatment group. 
    - File ex: C1_byTx_sorted_heatmap_by_t1_no_dendrogram_2024-10-16.pdf
- Sorts data by 2 groups using means of groups of interest, plots heatmap.
    - File ex: C2_byTx_sorted_heatmap_by_t1_t2_2024-10-16.pdf
- Removes boxes and sorts based on Ts top 50%, and plots heatmap.
    - File ex: C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-05.pdf
- Hierarchical clustering with dendrogram using 1 group and top 50% of genes.
    - File ex: C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-06_with_dendrograms.pdf
- Hierarchical clustering with no dendrogram using multiple groups and top 50% of variance.
    - File ex: hierarchical_heatmap_C3-C6_top50byVariance_2024-11-10.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 50 of Ts genes.
    - File ex: hierarchical_heatmap_C3-C6_top50T3_2024-11-10.jpg
- Hierarchical clustering with no dendrogram using multiple groups and all genes.
    *NOTE: Due to differing numbers of genes in each comparison, this code limits the genes displayed to include only the genes that are common in both datasets.
    - File ex: combined_heatmap_common_genes_no_dendrogram.jpg
- Code to display multiple heatmaps on the same page. 
    - File: Comparison_Heatmaps_2024-11-06.pdf
- Code lays the heatmaps side by side with no spaces, gene names, titles, or legends using the top 50 genes from each comparison.
    - File: combined_heatmap_no_genes.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 90 of 2N genes.
    - File ex: hierarchical_heatmap_C3-C6_top90-T1_2024-12-19.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 90% of 2N genes.
    - Cannot run based on differing number of rows between clusters. 
## 09_Sorted_Heatmaps_2024-12-19.R
- Runs the same code as 09_Sorted_Heatmaps_2024-11-10.R but uses Viridis color palette. 
- Files saved using the same syntax, but with updated date (2024-12-19). 


###############################################
###############################################
###############################################


# 10. Marker genes subset by cluster, grouped by treatment


## 10_projMCS5_markerList_byCluster-Tx.R
- Creates heatmap & z-scores for 4 treatment groups using all sig genes found in C18, subset by C18
- Generates file: 07_C18_All127_byTx_zscores_2024-03-13.csv


## 10_projMCS5_GeneMarkers_TxComp-byCluster_Heatmap.R
- Subsets by cluster
- Uses marker genes known to be significant for the cluster of interest
- Creates z-scores and heatmaps by treatment
- File names ex: C4_byTx_zscores_2024-03-21.csv & C4_TxComp_Heatmap_2024-03-21.pdf
## 10_projMCS5_C1_Sorted_Heatmap.R
- First sorted heatmaps created
## 10_projMCS7_GeneMarkers_byCluster-TxGrp_2024-10-17.R
- Upgrades from previous 2 code records include sorted heatmaps in 1 file
## 10_projMCS5_GeneMarkers_byClusterTxComp_Loop_2024-10-17.R
- Subsets by cluster, group by treatment; loop for all clusters
 
## 10_projMCS6_GeneMarkers_TxComp-byCluster_2024-06-26.R
- Adds Harmony, subsets by cluster, groupBy treatment, uses MAGIC
- cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5
- Ex. file name: C25_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
- Subset by cluster, group by sample: all returns were NULL
- Uses top 200 marker genes to create heatmaps & z-scores 
- Heatmap cutOff = "FDR <= 1 & Log2FC >=0.01"
- Ex. file names saved: C4_TxComp_Heatmap_2024-03-21.pdf & C4_byTx_zscores_2024-03-21.csv


## 10_projMCS5_C1_Sorted_Heatmap.R
- File loaded: C1_byTx_zscores_2024-03-21.csv
- Creates a dataframe and plots a heatmap.
        - File name ex: C1_byTx_sorted-pheatmap_2024-03-25.pdf
## 10_projMCS7_C1_Sorted_Heatmap_2024-10-16.R
- Expands upon code from 10_projMCS5_C1_Sorted_Heatmap.R
- Files loaded ex: C1_byTx_zscores_2024-03-21.csv
- Sorts data and creates heatmap.
    - File name ex: C1_byTx_sorted_heatmap_2024-10-16.pdf
- Sorted by treatment group, plots heatmap.
    - File name ex: C1_byTx_sorted_heatmap_by_t1_no_dendrogram_2024-10-16.pdf
- Sorted by 2 treatment group, plots heatmap.
    - File name ex: C2_byTx_sorted_heatmap_by_t1_t2_2024-10-16.pdf
- Removes boxes, sorts by top 50% in Ts group, and plots heatmap.
    - File name ex: C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-05.pdf
- Hierarchical clustering with dendrogram using 1 group and top 50% of genes to plot heatmap.
    - File name ex: C3_byTx_sorted_heatmap_by_t3_top_50_percent_2024-11-06_with_dendrograms.pdf
- Hierarchical clustering with no dendrogram using 1 group and top 50% of genes to plot heatmap. 
    - File name ex: C3_byTx_hierarchical_heatmap_1-Comp_Top50Percent_noDendrogramGenes_2024-11-10.pdf.pdf
- Hierarchical clustering with no dendrogram using multiple groups and top 50% of variance to plot heatmap.
    - File name ex: hierarchical_heatmap_C3-C6_top50byVariance_2024-11-10.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 50 of T3 genes to plot heatmap.
    - File name ex: hierarchical_heatmap_C3-C6_top50T3_2024-11-10.jpg
- Hierarchical clustering with no dendrogram using multiple groups and all genes to plot heatmap.
    - File name ex: combined_heatmap_common_genes_no_dendrogram.jpg
- Plots multiple heatmaps on the same page.
    - File name ex: Comparison_Heatmaps_2024-11-06.pdf
- Plots heatmaps side by side with no spaces, gene names, titles, or legends using the top 50 genes from each comparison.
    - File name ex: combined_heatmap_no_genes.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 90 ot T1 genes to plot heatmap.
    - File name ex: hierarchical_heatmap_C3-C6_top200-T1_2024-11-10.jpg
- Hierarchical clustering with no dendrogram using multiple groups and top 90% of 2N genes.
    - Does not run based on differing number of rows between clusters. 
        
## 10_projMCS7_GeneMarkers_byCluster-TxGrp_2024-10-17.R
- getMarkerFeatures subset by cluster, groupBy Sample, maxCells = 45000, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")
    - File name ex: C3_markerList_2024-10-17.csv
- Loop to run above code for all clusters
    - File name ex: cluster, "_GeneMarker_groupBy-Sample_", Sys.Date(), ".csv"
- Loop through all clusters
    - Subset by cluster, groupBy Treatment, maxCells = 45000, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
    - File name ex: cluster, "_GeneMarker_groupBy-Treatment_", Sys.Date(), ".csv"
- Loop through all clusters
    - Subset by treatment, groupBy Cluster, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
    - File name ex: T1_GeneMarker_groupBy-Treatment_2024-10-17.csv
    *NOTE: These results are saved under 09. Marker genes subset by treatment, groupBy Cluster.
- Loop through all clusters
    - No subset, groupBy Clusters, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
    - File name ex: cluster, "_GeneMarker_List_", Sys.Date(), ".csv"
- Creates dataframes for top genes in each cluster by Log2FC
    - Subset by cluster, groupBy treatment, maxCells = 45000
    - Plots heatmap of top genes by cluster with labeled genes of interest
        - cutOff = "FDR <= 1 & Log2FC >=0.01"
        - File name ex: C25_TxComp_Heatmap_2024-10-18.pdf
- Prints z-scores for the above cluster heatmaps
    - cutOff = "FDR <= 1 & Log2FC >=.01"
    - File name ex: C25_byTx_zscores_2024-10-19.csv
- Plots sorted heatmap of descending z-scores by the 2N treatment group
    - File name ex: C25_byTx_sorted_heatmap_by_t1_2024-10-18.jpg
- Plots sorted heatmap of descending z-scores by the Ts treatment group
    - File name ex: C25_byTx_sorted_heatmap_by_t3_2024-10-18.jpg
- Plots sorted heatmap of descending z-scores by the 2N & 2N+ treatment groups
    - File name ex: C25_byTx_sorted_heatmap_by_t1-t2_2024-10-18.jpg
- Plots sorted heatmap of descending z-scores by the Ts & Ts+ treatment groups
    - File name ex: C25_byTx_sorted_heatmap_by_t3-t4_2024-10-18.jpg
- Plots sorted heatmap of descending z-scores by the 2N+ & Ts+ treatment groups
    - File name ex: C25_byTx_sorted_heatmap_by_t2-t4_2024-10-18.jpg
- Cluster z-scores were analyzed for significance between treatment groups
    - Files loaded name ex: C1_byTx_zscores_2024-10-18.csv
    - Means calculated for each group
    - Difference calculated between means; differences Log2FC abs>=1.25 were marked significant and saved to file.
        - File name ex: C1_byTx_zscores_2024-10-18_analyzed.csv
    - Loop to analyze 2N & 2N+ vs. Ts & Ts+ differences
        - All results 
            - File name ex: C%d_t1-t2_t3-t4_zscores_2024-10-18_analyzed.csv
        - If sig differences were found:
            - File name ex: C%d_significant_t1-t2_t3-t4_zscores_2024-10-18.csv
    - Loop to analyze 2N+ & Ts+ vs. 2N & Ts differences
        - All results 
            - File name ex: C%d_t2-t4_t1-t3_zscores_2024-10-18_analyzed.csv
        - If sig differences were found:
            - File name ex: C%d_significant_t2-t4_t1-t3_zscores_2024-10-18.csv
    - Loop to analyze Ts vs. Ts+ differences
        - All results 
            - File name ex: C%d_t3_t4_zscores_2024-10-18_analyzed.csv
        - If sig differences were found:
            - File name ex: C%d_significant_t3_t4_zscores_2024-10-18.csv
- Combines Ts vs. Ts+ spreadsheets
    - Files loaded name ex: C%d_byTx_zscores_2024-10-18.csv
    - Significant genes from all clusters were combined.
        - File: Combined_Cluster_significant_t3_t4_comparison_2024-10-18.csv
- Combines 2N+ & Ts+ vs. 2N & Ts spreadsheets
    - Files loaded name ex: C%d_byTx_zscores_2024-10-18.csv
    - Significant genes from all clusters were combined.
        - File: Combined_Cluster_significant_t2-t4_t1-t3_comparison_2024-10-18.csv
- Combines 2N & 2N+ vs. Ts & Ts+ spreadsheets
    - Files loaded name ex: C%d_byTx_zscores_2024-10-18.csv
    - Significant genes from all clusters were combined.
        - File: Combined_Cluster_significant_t1-t2_t3-t4_comparison_2024-10-18.csv
- Tallies total clusters that share DA gene markers between treatment comps
    - File loaded: Combined_Cluster_significant_t3_t4_comparison_2024-10-18.csv
    - Significant entries tallied and converted to a data frame. 
    - Original data combined with total counts.
    - File: Combined_Cluster_Total-sig_t3_t4_comparison_2024-11-10.csv
- Tallies the total number of genes considered significant; prints tally.
    - File loaded: Combined_Cluster_Total-sig_t1-t2_t3-t4_comparison_2024-10-18.csv
- Tallies total clusters that share DA gene markers between tx comps within cell types
    - File loaded: Combined_Cluster_significant_t3_t4_comparison_2024-10-18.csv
    - Total sig entries for each cell type were tallied & barplot created.
    *NOTE: File saved under 11. Marker genes by cell type
        - Files: Cluster_DA-Gene_Counts_byCellType_t3-t4_comparison_2024-11-10.csv & Cluster_TOTAL-DA-Gene_Counts_byCellType_t3-t4_comparison_2024-11-10.jpg
- Creates a barplot to display cell type tallies for the number of DA genes observed in 1 cluster, 2 clusters, or 3 clusters
    - File loaded: Combined_Cluster_significant_t3_t4_comparison_2024-10-18.csv
    - Loops through all cell types and creates a plot.
        - File: Gene_Counts_by_Cluster_t3-t4-comp_CellTypeID_2024-11-10.jpg
- Loop to create heatmaps and print out z-scores for all clusters
    - Does not run.
    
## 10_q-2a-2b-common-genes_2024-11-08.R
- Loops through each cluster using load file ex: ", treatment, "_", cluster, "_GeneMarkers_2024-04-10.csv" & ", cluster, "_byTx_zscores_2024-03-21.csv"
    - Extracts Log2FC and z-scores for common genes between the 2 files.
    - Tallies common genes by treatment and cluster. 
    - Code identifies for when no common genes are found by treatment in a cluster. 
    - Final results file ex: final_results, "q-2a-2b_common-genes_2024-11-08.csv"
    - Gene tally file ex: q-2a-2b_gene_tally_by_treatment_cluster_2024-11-08.csv
- Creates barplot for gene tally.
    - Loads vfile: q-2a-2b_gene_tally_by_treatment_cluster_2024-11-08.csv
    - Sorts clusters numerically and uses Viridis color palette. 
    - File: q-2a-2b_gene_tally_by_tx-Cluster_Barplot_2024-11-08.png
- Creates a new spreadsheet of gene names based on results from above comparisons.
    - Loops through each cluster for the current treatment.
    - Loads file ex: ", treatment, "_", cluster, "_GeneMarkers_2024-04-10.csv" & ", cluster, "_byTx_zscores_2024-03-21.csv"
    - Extracts gene name and common gene names in both files.
    - Compares genes between every pair of clusters.
        - File ex: common_genes_across_clusters, "common_genes_between_clusters_by_treatment_2024-11-08.csv"
    - Creates a heatmap of the above outputs by looping through clusters. 
        - Code uncludes gene name results by cluster and treatment group. 
        - File: q-2a-2b_final_gene_matrix_2024-11-08.csv"

## 10_Volcano_byCluster-Tx.R
- Add Harmony & impute weights
- Subset by cluster; group by tx
- Pairwise comparisons
- CutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"

## 10_Dotplot.R
- getMArkerFeatures groupBy Cluster (no subset)
- Uses cell type-specific markers defined in the literature. 
- Dotplot, file name ex: Ggtree_Zu-Glut2_2024-07-17.pdf
- Combines plots, file name ex: astro_dotPlot_2024-07-22.jpeg

## 10_revised_DotPlot.R
- getMArkerFeatures groupBy Cluster (no subset)
- Similar code to 10_Dotplot.R but displays dotplots using Log2FC.
    - File name ex: gene_expression_dotplot_with_dendrogram_log2fc.png
## 10_DotPlots_2025-01-11.R
- Similar code to 10_revised_DotPlot.R but updates code to minimize less significant findings and emphasize more impactful findings. 
    - File name ex: astrocyte_dotplot_with_dendrogram_2025-01-15.png
## 10_DotPlots_2025-01-16.R
- Similar code to 10_DotPlots_2025-01-11.R but includes cell type-specific genes of interest. 
    - File name ex: glut_allUMAP_dotplot-dendrogram_2025-01-16.png


###############################################
###############################################
###############################################


# 11. Marker genes subset by cell type, grouped by treatment

## 11_projMCS7_GeneMarkers_byCellType-TxGrp_2025-01-16.R
- Similar code to projMCS7_GeneMarkers_byCluster-TxGrp_2025-01-16.R but focuses on cell type groups
- Files saved ex: astro_GeneMarker_List_2025-01-16.csv

## 11_CellType_TxComps.R
- Code did not produce meaningful results





**** GitHub & publication folder verification stopped here


###############################################
###############################################
###############################################


# 12. Cluster specific marker genes


12_projMCS5_ClusterID_UMAP.R
- Additional UMAPs created using genes of interest
- File generated: 07_UMAP-Marker-Genes-W-Imputation_2024-04-04.pdf

12_projMCS5_Zu-2023_geneMarker_UMAPs_2024-07-10.R
- Additinal UMAPs created using genes from Zu et al. 2023 publication
- Files generated ex: microglia_geneMarkers-Zu-2023_2024-07-10.pdf & projMCS6_UMAP-Clusters_2024-07-10.pdf


###############################################
###############################################
###############################################


# 13. Project MCS6 - Pairwise Comparisons


## Pairwise Comparisons by Cluster

## 13_Clusters_TvsTComp_2024-06-28.R
- addHaromony and addImputeWeights applied
- Subsets by cluster, grpBy treament 
- Uses MAGIC to impute gene scores
- Makes pairwise comparisons between all treatment groups
- cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
- Ex. file name: C1_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
13_T3vT124_byCluster_2024-08-09.R
13_T1-2-4vT3_TxComp_FDR-0-1_Log2FC-0-5_2024-07-31.R
- Both of these files carry out similar functions to 13_Clusters_TvsTComp_2024-06-28.R
* Pairwise comparisons from the above code were combined into spreadsheets by cluster
- Ex. combined file name: C25_MCS6-TxComp_geneMarkersCombined_2024-07-06
* Pairwise comparisons from the above code were combined into 1 spreadsheet
- Combined file name: TxComp_Markers_2024-09-05.xlsx
* Pairwise comparison summary tables were combined (MCS2023 & C18 subset)
- Combined summary file name: SummaryTables_2024-07-31.xlsx
* Pairwise comparison gene unions by cluster and pairwise comparison were analyzed
- Gene union spreadsheet ex. file name: C20_TxComp_GeneUnions_2024-07-24
## 13_Q4_TxComps-byCluster_2024-07-02.R
- Uses same code as 13_Clusters_TvsTComp_2024-06-28.R but also includes:
    - Combines tx group cluster-specific comparisons
        - File name ex: C18_All-T3vT_TxMarkers_Combined_2024-07-06.csv
    - Combines all cluster specific pairwise comparisons into 1 document
        - File name ex: C1_MCS6-TxComp_geneMarkersCombined_2024-07-06.csv
    - Combines, by cluster, DA genes for all cells in the cluster, in addition to DA genes from each subset treatment group. 
        - File name ex: C1_MCS6-TxMarkers_Combined_2024-07-05.csv
    - Sorts rows in the results file
        - File name ex: C18-ordered_PairwiseCombined_2024-07-08.csv
    - Q6.2 combines all files for cluster-specific pairwise comparisons; spreadsheets are in descending order from the least number of "NA" column for DA pairwise co parisons in each row of gene markers.
        - File name ex: C3_Ordered_PairwiseCombined_2024-07-08.csv

13_Q4_GenotypeTxComps_2024-06-28.R
- Adds Harmony
- Groups by genotype
- cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
- Files: GenotypeTxComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv & GenotypeCompByCluster_FDR-0-1_Log2FC-0-5_2024-06-28.csv

13_Q5_GenotypeClusterComps_2024_06-28.R
- Adds Harmony
- Subsets by genotype, groups by cluster
- Ex file names: C2_GenotypeComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv

13_projMCS6_Volcano_2024.08-08.R
- Add Harmony & impute weights
- Subset by cluster; group by tx
- Pairwise comparisons
- CutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"

13_DA-Genes_Total-Cells_ByCluster.R
- Creates barplot summarizing DA genes by pairwise comparison, and total cells per cluster. 
13_DA-Genes_Total-Cells_ByCluster_2024-12-16.R
- Updates code from 13_DA-Genes_Total-Cells_ByCluster.R to plot using Viridis color palette. 


## Pairwise Comparisons by Cell Type

## 13_CellTypes_TvsTComp_2025-01-18.R
- Subset by cell type, groupBy treatment, pairwise comparisons between treatment groups. 
    - Does not loop code, but has all options specifically coded in the document. 
    - File name ex: astro_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv


## Pairwise Comparisons by Genotype

## 13_genotype-comps_byCluster_Boxplot_2024-11-08.R
- Load file ex: file_path, pattern = "GenotypeComp_FDR-0-1_Log2FC-0-5"
- Dataframe created to store results; loop runs through each file, counts genotype row data, and creates a barplot.
- File ex: genotyp_gene-counts_byCluster_Barplot_2024-11-08.jpg 
## 13_genotype-comps_byCluster_Boxplot_2024-12-16.R
- Runs the same code as 13_genotype-comps_byCluster_Boxplot_2024-11-08.R but updated for Viridis color palette. 


###############################################
###############################################
###############################################


# 14. Gene Ontology (GO) and KEGG differentially accessible marker gene analysis. 


## 14_GO_Analysis_2024-08-01.R
- getMarkerFeatures, groupBy Cluster, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
- Calls top marker genes in each cluster.
- cc.go run for Gene Ontology results by cluster
- cc.gos run for simplified Gene Ontology results by cluster
## 14_projMCS5_GO_Analysis.R
- Runs same code as 14_GO_Analysis_2024-08-01.R

## 14_GO_Barplot.R
- Reads in pairwise comparison results by cluster. 
- enrichGO used to create GO plots, file name ex: C3_T1vsT3_top20.png
## 14_GO_Barplot_2.R
- Uses same code as 14_GO_Barplot.R but includes results

## 14_C18-Subset_GO_Barplot.R
- Reads in pairwise comparison results by C18 subclusters. 
- enrichGO used to create GO plots, file name ex: C18.5_T124v3_top20.png

## 14_GO_Plots2.R
- topGOdata performs GO analysis
    - Runs enrichment analysis with Fisher's exact test
    - Summarizes results and extracts top GO terms
    - creates barplots and dotplots
- goseq performs GO analysis
    - Summarizes top GO terms; creates barplots and dotplots
- getBM performs GO analysis
    - Summarizes top GO terms; creates barplots. 
- Creates custom dotplots. 
## 14_GO_Plots_Loop.R
- Uses combined code from 14_GO_Plots2.R
- Loops through pairwise comparisons by cluster and treatment group.

## 14_C18_GO_Plots_Loop.R
- Loops through all C18 treatment group pairwise comparisons. 
- enrichGO used to run Gene Ontology enrichment using clusterProfiler, ont = "BP"
- Creates a dataframe and then a barplot of results. 
    - File name ex: _top20_GO.png
- topGOdata used to run GO enrichment analysis.
    - File name ex: _TopGO_Top10_Barplot.png
- goseq used to run GO enrichment analysis.
    - File name ex: _GOseq_Top10_Barplot.png
- getBM (BioMart) to run GO enrichment analysis.
    - File name ex: _BioMart_Top10_GO_Count.png

## 14_GO_Analysis_2025-01-18.R
- getMarkerFeatures groupBy cellTypes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
- Loops through each cell type and saves marker list as a seperate document.
    - File name ex: cellType, "_GeneMarker_List_", Sys.Date(), ".csv"
- clusterProfiler used to obtain GO enrichment "BP", creates dotplots. 
    - Files saved as .png & .pdf, name: GO_Enrichment_Plot_2025-01-18.png & GO_Simplified_Enrichment_Plot_2025-01-18.png
    - Results saved as .csv files, name: GO_Enrichment_Results_2025-01-18.csv & GO_Simplified_Enrichment_Results_2025-01-18.csv
- Filters top 200 GO results
    - Files saved as .png & .pdf, name: GO_top200_Plot_2025-01-18.png & GO_top200_Simplified_Plot_2025-01-18.png", dotplot_top200_simplified, width = 13, height =15)
    - Results saved as .csv files, name: GO_top200_Results_2025-01-18.csv & GO_top200_Simplified_Results_2025-01-18.csv
    
## 14_GO_CellTypes_PairwiseComps_2025-01-18.R
- Loads file for each cell type and treatment group pairwise comparison seperately. 
    - enrichGO used for GO BP enrichment
        - File name ex: glut_T1vsT4_GO_results_2025-01-19.csv
    - extracts the top 20 GO categories by cell type and pairwise comparison for plotting
        - File name ex: cell_type, "_", treatment_comparison, "_top20_2025-01-19.png"
    - Loop to run through all cell types and treatment group comparisons.
        - Saves GO results for each cell type and treatment group comparison
            - File name ex: _GO_results_2025-01-19.csv
        - Plots top 20 GO categories for each comparison 
            - File name ex: _top20_viridis_2025-01-19.png
    - Includes a list of all returned results; or when no GO results are found. 
- Runs GO enrichment and creates plots for sig DA genes in all cell types and grouped treatment comparisons (not pairwise comparisons).
    - Loads file ex: Combined_CellType_significant_t1-t2_t3-t4_comparison_2025-01-18.csv
        - Combined_CellType_significant_t2-t4_t1-t3_comparison_2025-01-18.csv
        - Combined_CellType_significant_t3_t4_comparison_2025-01-18.csv
        - Combined_CellType_significant_t3_vs_all_2025-01-19.csv
    - Loops through each cell type column, performs enrichGO "BP", and creates barplots for results. 
        - File name ex: cell_type, "_GO_results.csv" & cell_type, "_GO_enrichment_plot.pdf"
    - Records when no GO results were returned. 
    
## 14_GO_GenesOfInterest_byCluster_Barplot_2024-11-08.R
- Creates barplots for genes of interest count by cluster.
    - Load file ex: GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv
    - Identifies duplicate genes amongst clusters and plots results.
    - File ex: GO-BP-GenesOfInterest_Duplicate_Genes_t1-t2_t3-t4_", date_today, ".csv" & Cluster_Gene_Counts_t1-t2_t3-t4_", date_today, ".jpg
- Creates GO barplots for genes of interest.
    - Load file ex: GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv
    - enrichGO BP ran for genes of interest, pAdjustMethod = "BH", qvalueCutoff = 0.05.
    - GO enrichment barplot created when GO enrichment results are available.
    - Code specifies when gene mapping failed for a particular cluster and comparison.
    - File ex: col_name,"_GO_GenesOfInterest_byCluster_Barplots_t1-t2_t3-t4_2024-11-08.jpg
- Creates GO DotPlots for genes of interest.
    - Load file ex: GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv
    - enrichGO BP ran for genes of interest, pAdjustMethod = "BH", qvalueCutoff = 0.05.
    - GO enrichment dotplot created when GO enrichment results are available.
    - Code specifies when gene mapping failed for a particular cluster and comparison.
    - File ex: col_name, "_GO_GenesOfInterest_Dotplot_t1-t2_t3-t4_2024-11-08.png"

## 14_clusterProfiler_dotplots_2024-11-23.R
- Loaded ex: Combined_Cluster_significant_t2-t4_t1-t3_comparison_2024-10-18.csv
- enrichGO for BP, CC, and MF ontologies, pvalueCutoff = 0.05, qvalueCutoff = 0.05.
- Arranges file names for saving, saves output file and generates dotplot. 
    - File name ex: sanitized_name, "_TxComp_GO_results_2024-11-23.csv & sanitized_name, "_dotplot_2024-11-23.png
- Code includes results of pairwise comparisons and number of genes if applicable. 
## 14_clusterProfiler_Simplified_dotplots_2024-11-23.R
- enrichGO runs the same analysis as 14_clusterProfiler_dotplots_2024-11-23.R, but simplifies results for each comparison. 
    
## 14_GO_GenesOfInterest_Barplot_2024-11-24.R
- Creates barplots showing the number of GO genes of interest per cluster, by pairwise comparison. 
- Loads file ex: GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv
- Counts the number of counts in each column; removes rows with NA, empty, or NULL. 
- Creates barplot.
    - File name ex: GO_GenesOfInterest-byCluster_t3_t4_BarPlot_2024-11-24.png
    
## 14_GO_GenesOfInterest_byCluster_Barplot_2024-12-16.R
- Creats barplots for genes of interest count by cluster. 
- Loads file ex: GO_Genes-Of-Interest_t2-t4_t1-t3_2024-11-08.csv
    - Identifies duplicate genes between clusters and calculates total counts for unique genes in each cluster. 
    - File name ex: GO-BP-GenesOfInterest_Duplicate_Genes_t2-t4_t1-t3_", date_today, ".csv"
- Loops through all pairwise comparison GO genes of interest files, maps gene symbols to Entrez IDs, checks for GO enrichment results.
    - enrichGO BP ontology ran with pAdjustMethod = "BH" (Benjamini-Hochberg), qvalueCutoff = 0.05
    - Creates barplots for top 10 GO categories.
    - Code includes cluster pairwise comparisons that failed gene mapping. 
        - File name ex: col_name,"_GO_GenesOfInterest_byCluster_Barplots_t3_t4_2024-12-16.jpg"
- DotPlots created for GO genes of interest.
    - Loads file ex: GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv
    - Loops through each cluster and maps genes to Entrez IDs. 
    - enrichGO BP ontology ran with pAdjustMethod = "BH" (Benjamini-Hochberg), qvalueCutoff = 0.05
    - Code includes cluster pairwise comparisons that failed gene mapping. 
    - File name ex: col_name, "_GO_GenesOfInterest_Dotplot_t1-t2_t3-t4_2024-12-16.png"
    
## 14_GO_match-sig-genes_2024-10-29.R
- Combines GO results with treatment group accessibility levels, for one cluster
    - Loads 2 files, name ex: C1_significant_t1-t2_t3-t4_zscores_2024-10-18.csv & C1_GO_t1-t2_t3-t4_2024-10-18.csv
    - Joins data and saves file, name ex: C1_GO_t1-t2_t3-t4_matched_genes_2024-10-29.csv
- Loops through all clusters for specific pairwise comparisons: 2N/2N+ vs Ts/Ts+, 2N+/Ts+ vs 2N/Ts, Ts vs Ts+, 2N/2N+ vs Ts/Ts+, 2N+/Ts+ vs 2N/Ts, Ts vs Ts+
    - File ex: base_path, "C", i, "_GO_t1-t2_t3-t4_matched_genes_2024-10-29.csv"
    
## 14_GO_CellType-Analysis_2025-01-18
- getMarkerFeatures groupBy cellTypes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25". 
    - Loops through each cell type and saves marker list.
    - File name ex: cellType, "_GeneMarker_List_", Sys.Date(), ".csv"
- Loops through cell types 
    - clusterProfilier cc.go ran with pvalueCutoff = 0.05, pAdjustMethod = "BH"
    - clusterProfiler cc.gos simplified ran with cc.go, cutoff=0.05, by= "p.adjust"
    - Creates dotplots for top 10 categories, saved as both png and pdf. 
    - Plot file ex: GO_Enrichment_Plot_2025-01-18.png & GO_Simplified_Enrichment_Plot_2025-01-18.png
    - File ex: cc_go_df, "GO_Enrichment_Results_2025-01-18.csv" & "cc_gos_df, "GO_Simplified_Enrichment_Results_2025-01-18.csv"
- Obtains the top 200 markers for each cell type by Log2FC
    - runs the same code as loop but filters by top 200
    - File ex: GO Top 200 abs(Log2FC) Enrichment Across Cell Types & GO Top 200 abs(Log2FC) Simplified Enrichment Across Cell Types
    - Creates dotplots, file ex: GO_top200_Plot_2025-01-18.pdf & GO_top200_Simplified_Plot_2025-01-18.pdf
    - File ex: cc_go_top200_df, "GO_top200_Results_2025-01-18.csv" & GO_top200_Simplified_Results_2025-01-18.csv
    
## 14_GO_CellTypes_PairwiseComps_2025-01-18
- Completed for each cell type and treatment group pairwise comparison. 
    - Load file ex: glut_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2025-01-18.csv
    - Filters genes abs(Log2FC) >= 0.5
    - enrichGO ran for BP
    - File ex: GO_results_df, "glut_T1vsT4_GO_results_2025-01-19.csv"
    - Top 20 GO categories plotted
    - File ex: cell_type, "_", treatment_comparison, "_top20_2025-01-19.png"
- Loop runs through all cell types and treatment group comparisons
    - Code includes results of if genes met thresholds for cell type and tx comp. 
    - File ex: (treatment_comparison), "_GO_results_2025-01-19.csv"
    - Top 20 GO categories plotted, file ex: (treatment_comparison), "_top20_viridis_2025-01-19.png"
- Loop generates GO terms and plots for significant DA genes in all cell types and grouped treatment comparisons.
NOTE: Does not use pairwise comparisons. 
    - Load file ex: Combined_CellType_significant_t1-t2_t3-t4_comparison_2025-01-18.csv
        - Includes code for all comparisons of interest, including Combined_CellType_significant_t3_vs_all_2025-01-19.csv
    - Code includes results for cell type treatment group comps that did not return GO enrichment terms. 
    - File ex: sanitized_cell_type, "_GO_results.csv" & sanitized_cell_type, "_GO_enrichment_plot.pdf"
    
## 14_KEGG_Analysis_2025-01-25.R
- Analyzes cluster specific pairwise comarison results for KEGG enrichment pathways.
    - Load file ex: C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv
    - dplyr filters upregulated (Log2FC > 0) and downregulated (Log2FC < 0) genes, converts gene names to Entrez IDs. 
    - enrichKEGG run for both entrez_up and entrez_down; pvalueCutoff = 0.05, qvalueCutoff = 0.2.
    - File ex: C3_T1vsT4Comp_KEGG_Upregulated.csv & C3_T1vsT4Comp_KEGG_Downregulated.csv
    - KEGG Cytoscape-stype network plots, file ex: C3_T1vsT4Comp_Cytoscape_Network_", Sys.Date(), ".pdf"
- Loop runs through all pairwise comparisons in directory. 
    - Uses same code as above.
    - File ex: file_base, "_KEGG_Upregulated.csv" & file_base, "_KEGG_Downregulated.csv"
    - Generates Cytoscape network export using Rcy3
    * DOES NOT RUN
- Combines spreadsheets from loop into one document.
    - Code includes results or if no genes could be mapped for a cluster specific pairwise comparison. 
    - File: KEGG_Combined_Results.csv
    
## See summary files: 
- GO_Cluster_vs_CellType_Pairwise_TxComps_2025-01-23.xlsx
- GO_Cluster_vs_CellType_byTxComps_2025-01-22.xls
- KEGG_Results_byCluster_2025-01-25
- KEGG_Combined_Results.csv


###############################################
###############################################
###############################################


# 15. ProjMCS7 - Peak Matrix creation and motif annotations

## 15_projMCS7.R
Peak enrichment
- getPeakSet
- addPeakMatrix
- Creates projMCS7
- Marker peak cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
- Plots heatmap of all clusters, file: Peak-Marker-Heatmap_2024-08-07.pdf
- Plots MA and Volcano plots, file: C18-PeakMarkers-MA-Volcano_2024-08-07.pdf
- Plots Browser tracks, file: C18-Tracks-With-PeakFeatures_2024-08-07.pdf
- Pairwise testing between groups; Plots MA & Volcano plots
    - File: C18-vs-C22-Markers-MA-Volcano_2024-08-07.pdf
- addMotifAnnotations; motifSet = "cisbp"
- peakAnnoEnrichment; cutOff = "FDR <= 0.1 & Log2FC >= 0.5"(C18)/Log2FC <= -0.5"(C22)"
- Creates a dataframe object containing the motif names, corrected p-values, and significance rank. 
- Plots rank-sorted TF motifs, colored by significance of enrichment
    - File: C18-vs-C22-Markers-Motifs-Enriched_2024-08-07.pdf
- Plots an enrichment heatmap, file: Motifs-Enriched-Marker-Heatmap_2024-08-07.pdf
- addArchRAnnotations
- Plots top C18 (Log2FC >= 0.5) and C22 (Log2FC <= -0.5)
    - File: C18-vs-C22-Markers-Motifs-Enriched_2024-08-07.pdf
Motif Enrichment
- peakAnnoEnrichment; seMarker = markersPeaks, ArchRProj = projMCS7, peakAnnotation = "Motif"
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
    - Plots heatmap: Motifs-Enriched-Marker-Heatmap_2024-08-07.pdf
    - Code is not supported collections: EncodeTFBS, ATAC, Codex; or ChIP peakAnnotations
- addBgdPeaks; used in computing deviations. 
- addDeviationsMatrix; computes per-cell deviations across all motif annotations.
- getVarDeviations; used to access deviations & and return a ggplot object. 
    - File: Variable-Motif-Deviation-Scores_2024-08-07.pdf
- getFeatures; extracts motifs for downstreatm analysis. 
- grep; used to get just the features corresponding to z-scores. 
- Plots groups with impute weights, file: Plot-Groups-Deviations-w-Imputation_2024-08-07.pdf
- Plots z-score distributions on UMAPs.
- Cannot run code for markerRNA. 
## 15_projMCS7_2024-12-20.R
- Runs the same code as projMCS7.R but with a few notes added. 

## 15_DA_Motif_Plot_byCluster.R
- Loads motifs.csv to create a heatmap of DA Motif Count Comparison Across Clusters.
- File: DA_Motif_Count_Comparison_2024-09-27.png
## 15_DA_Motif_Plot_byCluster_2024-12-16.R
- Runs same code as 15_DA_Motif_Plot_byCluster.R
- File: DA_Motif_Count_Comparison_2024-12-16.png

## 15_Motifs_byTx.R
- Subsets projMCS7 by treatment, getMarkerFeatures groupBy Clusters, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
- Plots an enriched motifs heatmap for each treatment group subset.
    - File name ex: t4-Motifs-Enriched-Marker-Heatmap_2025-01-11
    
## 15_Motifs_byTx-CellType.R
- Assigns treatment and cell type
- getMarkerFeatures using PeakMatrix and groupBy cell_treatment
- Barplot displays single cell counts by cell type and treatment.
    - File: Cell-Cts_byCellType-Tx_2025-01-13.pdf & cellCount_byTx-CellType_2025-01-14.csv
- Creates boxplots for each .csv file with individual points, mean line, and Viridis color palette. 
    - Loaded file ex: microglia.csv
- Plots boxplots side by side & in a line
    - File: CellCount_byTx-CellType_2025-01-15.pdf
- Performs ANOVA tests for significant differences between treatment groups. 
    - File: CellCount_byTx-CellType_2025-01-15_with_ANOVA_p_values_only.pdf
- Makes treatment group pairwise comparisons by treatment group, and plots motifs.
    - FDR <= 0.1 & abs(Log2FC) >= 0.5
    - File name ex: t1_t3_glut_-Motifs-Enriched-Marker-Heatmap_2025-01-13
- Loop to process all cell types and treatment group comparisons
    *NOTE: No enrichments found in any comparisons made. 


15_normExp_projMCS7.R
- getMarkerFeatures uses PeakMatrix, groupBy clusters
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
- Creates Browser tracks of interesting genes 
    - groupBy Sample file: normTracksbySample-Interesting-With-PeakFeatures_2024-08-28.pdf
    - groupBy cluster file: normTracksbyClusters-Interesting-With-PeakFeatures_2024-08-28.pdf
    - groupBy treatment file: C3-normTracksbyTx-With-PeakFeatures_2024-08-28.pdf
    - Examines C18 subclusters


15_peakVSmatrix.R
- Follows code from projMCS7.R, but loops through subset clusters to make tx group pairwise comparisons
    - File ex: _Motif_FDR-0-1_2024-09-16.csv
- Creates Browser tracks for genes of interest
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
    - File: C18-Tracks-With-PeakFeatures_2024-08-07.pdf & Peak-C3-T1vsT3-Tracks_2024-09-16.pdf
    

15_Peak_TxComps-byCluster_2024-09-09.R
- Loops through all clusters to make treatment group pairwise comparisons
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
    - File name ex: cluster, "_", useGroup, "vs", bgdGroup, "_Peak_FDR-0-1_Log2FC-0-5_2024-09-09.csv"
- Combines treatment group cluster-specific comparisons
    - File: C18_All-T3vT_TxMarkers_Combined_2024-07-06.csv
    - File ex: C1_MCS6-TxComp_geneMarkersCombined_2024-07-06.csv
    - Sorts rows in results file, file: C18-ordered_PairwiseCombined_2024-07-08.csv
- Research Question 6.2 combines all files for C3 pairwise comparisons.
    - File: C3_Ordered_PairwiseCombined_2024-07-08.csv
    
15_Peak_heatmaps_2024-09-11.R
- subset by Cluster, groupBy treatment; peakAnnoErichment for heatmap did not run.

15_projMCS5_Subset_Clusters.R
- Tries to subset clusters and run DA analysis on marker peaks.
- Does not run. 


###############################################
###############################################
###############################################


## 16. Fragment count analysis of peak locations

16_peaksVSmotifs.R
- Similar to 15_peakVSmatrix.R code, but shows examples of removing elements from Browser track plotSummary, such as:
    - bulkTrack, scTrack, featureTrack, geneTrack
- Code to loop through all clusters, identifying specific genomic regions (ex: chr 16 & 17 mm10 triplicated regions)
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
    - File ex: Peak-C3_T1vsT3-C3_Tracks_2024-09-17.pdf
- Obtaining nFrag counts for regions of interest
    - Specifies regions; creates a dataframe. 
    - Uses indices of cells corresponding to each sample; extracts peak matrix from sample; gets row ranges; subsets peak ranges by regions of interest; gets indices of overlapping peaks and subsets peak matrix using these indices; calculates total fragments for each cell; sums total fragments across all cells for the current sample; adds results to the df.
    - File: region_name, "_peakFrag_counts_2024-09-17.csv
- Combines .csv files and merges them into a new dataframe. 
    - Creates a treatment column based on sample column; groups and reorders data by treatment.
    - File: combined_results, "combined_peakFrag_counts_byTx_2024-09-17.csv"
- Combines files and normalizes by fragment counts. 
    - Normalizes by dividing the number of gene frags by total number of frags in sample. 
    - File ex: sample_stats, "combined_normPeaks_statistics_with_treatment_2024-09-17.csv"
- Loop code does not run.


16_BetterBrowserTracks.R
- addCoAccessibility & color palette updates to previous browser track code
- Loops through all clusters and treatment comparison regions of interest. 
- plotBrowserTrack used to assess fragment coverage depth across different genomic regions.
- File name ex: "Coverage_", region_name, "_", Sys.Date(), ".pdf"
- getCounts() or getGroupBW() allows for direct extraction of fragment counts for regions of interest. 
- addGroupCoverages can precompute the group coverages for each treatment group.
    - This can then be visualized on Browser tracks. 
- getGroupBW can also be used to export BigWig files normalized by fragment count (normMethod = "nFrags"). 
    - BigWig files can then be loaded into genome browser (IGV) to view coverage depth. 
- getFraments is used to extract raw fragment data for genomic regions of interest. 


16_fragCounts_GenomeWide_TxComps.R
- Loads fragments by treatment using sample, "_fragments.tsv.gz" files.
    - Adds treatment label; loads fragments for each treatment group. 
    - Creates bin and fragment counts by chromosome.
    - Plots fragment count differences between treatment groups. 
    - Can adjust for multiple chromosomes. 
- Loop to make pairwise tratment group comparisons using chr 16 & 17.
    - File name ex: Fragment_Count_Difference_t1_vs_t2_", chr, "_2024-09-25.png
- Plots fragment counts of rank distributions by chromosome and treatment.
    - Saved by chr, file ex: "rank_distribution_", treatment, "_", chr, "_2024-09-25.png"
    - Rank distro by fragment counts and treatment, file ex: "Rank Distribution of Fragment Counts -", treatment"


16_fragCounts_byPeak.R
- To save the number of frags by cluster:
    - Uses getMarkerFeatures, groupBy Clusters, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
    - Creates dataframe of marker peaks by cluster, file name ex: markersPeaks_List_", clusterName, "_2024-09-24.csv
    - getGroupSE extracts group-wise fragment counts for peaks
    - Fragment counts, grouped by clusters, with peak information. 
        - File: fragmentCounts_PeakMatrix_2024-09-24.csv
- To save the number of frags by treatment and cluster
    - groupBy Clusters_Treatment
    - File: fragmentCounts_Clusters_Treatment_2024-09-24.csv
- Plots peak frag counts by cluster.
    - File: Density_RankPlot_Cluster_", cluster, "_2024-09-26.png
- Code to plot rank orders of fragment counts by cluster and treatment group
    - Doesn't actually plot the frags in rank order
    

16_peakFrags.R
- Creates an index of individual sample cells and defines regions of interest
- getMatrixFromProject extracts the PeakMatrix using rowRanges and subsetByOverlaps by regions of interest. 
- Subsets by sample and sums fragments
- Creates a dataframe by looping over each sample. 
    - Includes number of single cells, median fragments, median TSS Enrichment, total fragment counts.
    - File: sample_statistics.csv
- Loop first runs over all samples and 1 region of interest (file: app_fragment_counts.csv), then all samples and all chr 16 & 17 regions of interest.
    - File name ex: region_name, "_fragment_counts.csv"
16_geneFrags.R
- Duplicate code to 16_peakFrags.R
16_geneFragsLoop.R
- Alternate loop code to: 16_geneFrags.R
- No files generated, for code reference only.
16_peakFragsLoop.R
- Additional alternate loop code to: 16_peakFrags.R
- No files generated for code reference only.


16_getMatrix_MCS1.R
- Extracts and saves RDS files for tile and gene matrices. 
- Code to extract data from these matrices. 
- Combines extracted data and creates data frames. 
    - File: combined_TileData_projMCS1.csv & combined_GeneData_projMCS1.csv
    
    
## 16_nFrags-TSS_byTx-Cluster.R
- Subset by cluster, groupBy treatment.
- Violin plot colorBy cell data using log10(nFrags).
    - File: C25_nFrags-byTx_2024-08-13
- Violin plot colorBy cell data using TSS Enrichment.
    - File: C25_TSS-byTx_2024-08-13
## 16_nFrags-TSS_byTx-Cluster_2024-12-16.R
- Runs the same code as 16_nFrags-TSS_byTx-Cluster.R but also includes loop for all clusters, printed as 1 pdf. 
    - File: QC-TSS-UniqFrags_byTx-Cluster_ViolinPlots_2024-12-16.pdf
    

16_normPeakStats_2024-09-08.R
- Combines individual files into one to normalize fragment counts.
- Normalized by dividing the number of gene fragments by total number of fragments in sample. 
- Uses file pattern = "fragment_counts.*\\.csv") & sample_statistics.csv
- Not used to generate results; for code reference only. 

16_peak_fragCounts_byTx_2024-09-18.R
- Loop for treatment group pairwise comparison using multiple genomic regions of interest. 
- Loop creates a dataframe to store treatement group indexes, extracts PeakMatrix, rowRanges, subsets by peak overlaps, and sums total fragments. 
    - File name ex: region_name, "_peakFrag_counts_byTx_2024-09-18.csv"
- Loop calls previous files created, combines results by treatment group. 
    - File name ex: combined_results, "combined_peakFrag_counts_byTx_2024-09-17.csv"
- Loop calls previous files created and sample_statistics.csv; groups samples by treatment, and merges data from the 2 files. 
    - Ratio for each fragment count added to a new column in the dataframe and reordered by treatment. 
    - File: combined_normPeaks_statistics_with_treatment_2024-09-17.csv


16_CSVtoBED.R
- Converts .csv files to .bed files to load into IGV viewer for analysis. 
    

###############################################
###############################################
###############################################


# 17. ProjMCS9 Peak matrix pairwise comparisons

## 17_projMCS9.R
- Loops through all clusters making pairwise comparisons using peakMatrix
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
    - File ex: cluster, "_", useGroup, "vs", bgdGroup, "_Peak9_FDR-0-1_Log2FC-0-5_2024-09-12.csv"
- addGroupCoverages; to create pseudobulk replicates.
- addReproduciblePeakSet; getPeakSet
- Creates projMCS9_5
- Loops through all clusters; adds impute weights; loops through pairwise comparisons
    - cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
    - File name ex: cluster, "_", useGroup, "vs", bgdGroup, "_Peak9-5_FDR-0-1_Log2FC-0-5_2024-09-13.csv
- Code repeated to create projMCS9_7; reran pairwise comparison loops. 


###############################################
###############################################
###############################################


# 18. Cluster 18 Subset & Reclustering

## 18_C18_Subset.R
- subsetArchRProject subsets cluster 18 cells from all other cells. 
- addIterativeLSI, addHarmony, addClusters completed for the subset.
- Confusion matrix created, file: C18ArchRSubset_Clusters_2024-06-19.csv
- UMAP and tSNE single cell embeddings plotted. 
    - Files: UMAP-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf & TSNE-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf
- Harmony added to single cell embeddings.
    - Files: UMAP2Harmony-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf & TSNE2Harmony-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf
- getMarkerFeatures groupBy clusters; cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
    - File name ex: i, "_C18ArchRSubset_2024-06-19.csv"
- UMAPs for genes of interest created, file: GeneScores-Marker-Heatmap_C18ArchRSubset_2024-06-19.pdf
- UMAP for genes of interest before imputation, file: UMAP-MarkerGenes-WO-Imputation_C18ArchRSubset_2024-06-19.pdf.
- UMAP for genes of interest after applying MAGIC to impute gene scores, file: UMAP-MarkerGenes-W-Imputation_C18ArchRSubset_2024-06-19.pdf
- Additional UMAPs created for genes of interest: UMAP-RadialGliaMarkerGenes-W-Imputation_C18ArchRSubset_2024-06-21.pdf
- UMAP highlights Ts cells, file: UMAP_T3_C18ArchRSubset_2024-06-21.pdf.

## 18_C18Subcluster_Combine_Spreadsheets.R
- Loads file ex: C1_C18ArchRSubset_2024-06-19.csv
- Combines all subcluster results into 1 document.
    - File: C18Subcluster_TxMarkers_Combined_2024-06-19.csv

## 18_C18_Subset_TxComps_byCluster.R
- Harmony added, subset by subcluster, assigns treatment groups. 

## 18_C18_TxComps_byCluster_2024-07-24.R
- Addresses Research Question 3
    - Harmony added, subset by subcluster, treatment group assigned, getMarkerFeatures groupBytreatment, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5").
    - File name ex: C18.18_TxComp_FDR-0-1_Log2FC-0-5_2024-07-25.csv
- Addresses Research Question 5
    - Harmony added, subset by subcluster, treatment group assigned, getMarkerFeatures by treatment group pairwise comparison, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5").
    - File name ex: C18.18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv
- Reads pairwise comparison spreadsheets into a dataframe. 
    - File name ex: C18.18-TxComp_geneMarkersCombined_2024-07-25.csv
## 18_C18_T1-2-4v3_byCluster_2024-07-31.R
- Follows code workflow from C18_TxComps_byCluster pairwise comparison, but analyzes specific differences between combined groups T1, T2, and T4 compared to T3.

## 18_C18_DA-Genes_Total-Cells_ByCluster.R
- Creates barplot summarizing subset cluster 18 DA genes by pairwise comparison, and total cells per subcluster. 
## 18_C18_DA-Genes_Total-Cells_ByCluster_2024-12-16.R
- Updates code from 18_C18_DA-Genes_Total-Cells_ByCluster.R to plot using Viridis color palette. 


###############################################
###############################################
###############################################


# 19. Cytoscape

## 19_Cytoscape_Analysis_2025-01-25.R
- Code does not run. 

## 19_Cytoscape_Analysis_2025-01-31.R
- Code does not run. 


###############################################
###############################################
###############################################


# 20. Genomic-Behavioral Correlation Statistics

## 20_PearsonCorrelations_2024-11-14.R
- Plots pearson correlations comparisons septerately, includes all 4 treatment groups. 
- Defines behavior variable and gene variable; loops through each behavioral variable and gene variable; loops through each treatment group to calculate the correlation.
    - cor correlation calculation uses method = "pearson"
- Creates a scatterplot of correlations for the current behavior and gene pair across treatment groups.
    - File ex: C19_pearson-correlation_plots/", behavior_var, "_vs_", gene_var, "_2024-11-14.png
- Loop plots individual person correlation comparisons on the same page, includes all 4 treatment groups. 
    - File ex: C19_all_genes_correlation_plots
- Step by step debugging strategies.
- T1 vs T3 pearson correlations; plot with task on x-axis & gene accessibility on y-axis
    - File ex: C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_t3_gene_correlation.png"
    - Saves correlational results, file ex: C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_t3_correlations.csv
- Scatterplot code for T1 on x-axis and T3 on y-axis, comparing gene z-scores. 
    - File ex: C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_vs_t3_scatter_2024-11-14.png
    - Saves correlational results, file ex: C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_vs_t3_gene_accessibility.csv
- Scatterplot of T1 on x-axis and T3 on y-axis, comparing behavioral test scores, one plot for each comparison. 
    - File ex: C19_t1_t3_behavior_correlation_plots/", behavior_var, "_t1_vs_t3_behavior_scores.csv
- Scatterplot of T2 & T4 on x-axis and T1 & T3 on y-axis, comparing gene z-scores, one plot for all comparisons. 
- Loop breaks up and makes scatterplots of treatment group comparisons (instead of using all treatment groups, as above). 
    - Load file ex: C19_significant_t2-t4_t1-t3_zscores_2024-10-18.csv
    - Filters data to include only common genes.
    - File ex: scatterplot_", x_group, "_vs_", y_group, ".png
- Streamlines the above code to read in only the columns of interest and auto generates the gene-vars list, then prints all plots in one pdf document. 
    - File to save all plots, ex: C19_sig_t2-t4_t1-t3_zscores_scatterplots_2024-11-19.pdf
- Pearson correlation values between "Task_B_Correct" and gene expression values across treatments (t1 and t3), with t1 on the x-axis and t3 on the y-axis.
    - Does not run. 
    
## 20_Cluster_PairwiseComps_Pearson-withFDR_2025-02-10.R
- Combines spreadsheets for cluster treatment group pairwise comparison z-scores (C3, C8, C11, C14, C18, C22) and all behavior stats (A-45, B-59, C-48, D-149)
- Load file ex: TaskA-45_AllStats.csv & C3_byTx_zscores_2024-10-18.csv
- Creates scatterplots, calculates correlations and p-values using cor.test, calculates FDR-adjusted p-values using method = "BH".
    - File ex: C3_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv
- Plots modified to include FDR-adjusted p-values and groups results by cluster. 
    - File ex: C3-zscore_TaskD-149_scatterplots_withFDR_2025-02-12.pdf
- Loop code does not run. 
## 20_Cluster_Pearson-withFDR_2025-02-10.R
- Uses same code as 20_Cluster_PairwiseComps_Pearson-withFDR_2025-02-10.R but analyzes cluster treatment group comparisons (not pairwise comparisons). 
    - Load file ex: C3_TaskD-149_StatsCombined_2025-02-03.csv
    - File ex: C3_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv
    - Plot, file ex: C3-zscore_TaskD-149_scatterplots_withFDR_2025-02-12.pdf
    
## 20_Genescore-Behavior_Scatterplots_PearsonCorr_2025-02-03.R
- Combines spreadsheets for downstream anaysis.
    - Load file ex: TaskD-149_AllStats.csv & C22_byTx_zscores_2024-10-18.csv
    - File ex: C22_TaskD-149_StatsCombined_2025-02-03.csv
- Creates scatterplots, calculates correlation using cor.test.
    - Plot, file ex: C22-zscore_TaskD-149_scatterplots_2025-02-03.pdf
    - Summary file ex: C22_TaskD-149_CorrelationStats_byTreatment.csv
- General Interpretation of Pearson R results. 
- Correlations run for clusters (C3, C8, C11, C14, C18, C22) and behavior stats (A-45, B-59, C-48, D-149).
    - Files to combine ex: TaskA-45_WrongStats.csv & C3_byTx_zscores_2024-10-18.csv
    - Files ex: C22_TaskD-149_WrongStatsCombined_2025-02-03.csv
    - Creates scatterplots, calculates correlation using cor.test.
    - Plot, file ex: C22-zscore_TaskD-149-Wrong_scatterplots_2025-02-03.pdf
    - Summary file ex: C22_TaskD-149-Wrong_CorrelationStats_byTreatment.csv
- Searches through Pearson R correlation results for values abs(0.4). 
    - Load file ex: C3_TaskC-48-Wrong_CorrelationStats_byTreatment.csv
    - File ex: C3_high_value_genes.csv
- Loop runs through all files in the directory to search for "high value genes."
    - File ex: file_identifier, "_high_value_genes.csv"
    
## 20_Treatment-PearsonCorrelations.R
- Load file ex: C3_TaskD-149_StatsCombined_2025-02-03.csv
- Performs treatment-specific correlational analysis using cor.test. 
- File ex: C3_TaskD-149_TreatmentSpecific_CorrelationStats_2025-02-03.csv
- Plot file ex: C3-zscore_TaskD-149_TreatmentSpecific_scatterplots_2025-02-03.pdf


## 20_Genescore-Behavior_CellType_Scatterplots_PearsonCorr_2025-02-03.R
- Run for each cluster and "AllStats" behavioral study. 
- Combines spreadsheets for downstream analysis. 
    - Load files ex: TaskA-45_AllStats.csv & gaba_byTx_zscores_2025-01-16.csv
    - File ex: astro_TaskD-149_StatsCombined_2025-02-03.csv
- Creates scatterplots, uses cor.test to calculate correlation, then plots all genes and behavioral matrics. 
    - File ex: astro-zscore_TaskD-149_scatterplots_2025-02-03.pdf
    - Saves correlation summary, file ex: astro_TaskD-149_CorrelationStats_byTreatment.csv
    - Can print correlation summary with p-values, file ex: glut_TaskA-45_CorrelationStats-pValues_byTreatment.csv
- General interpretation guidelines for Pearson R correlation results. 
- Loop run for each cluster and "WrongStats" behavioral study
    - Load file ex: TaskA-45_WrongStats.csv & gaba_byTx_zscores_2025-01-16.csv
    - Saves combined file, ex: astro_TaskC-48_WrongStatsCombined_2025-02-03.csv
- cor.test calculates correlations and creates plots for all genes and behavioral metrics. 
    - File ex: astro-zscore_TaskC-48-Wrong_scatterplots_2025-02-03.pdf
    - Correlation summary, file ex: astro_TaskC-48-Wrong_CorrelationStats_byTreatment.csv
- Searches through Pearson R correlation results for values abs(0.4)
    - Load file ex: endo_TaskC-48-Wrong_CorrelationStats_byTreatment.csv
    - Results saved, file ex: high_value_rows, "high_value_genes.csv"
- Loop runs through all files to search for "high value genes"
    - Results saved, file ex: file_identifier, "_high_value_genes.csv"
- Code tries to save all "high value gene" results into 1 file, but it doesn't run correctly

## 20_CellType_Pearson-withFDR_2025-02-10.R
- Completed for each major cell type (astro, gaba, glut, oligo, endo, microglia) and behavioral test (TaskA-45, TaskB-59, TaskC-48, TaskD-149). 
- Load file ex: astro_TaskA-45_StatsCombined_2025-02-03.csv
- Creates scatterplot, calculates correlations using cor.test and adjusts p-values using method = "BH". 
- File ex: astro_TaskA-45_CorrelationStats_withFDR_2025-02-12.csv
- Plots modified to include FDR-adjusted p-values, file ex: astro-zscore_TaskA-45_scatterplots_withFDR_2025-02-12.pdf
- Loop code does not run.

## 20_Lasso_Correlation_Stas_2025-01-29.R
- Subsets projMCS7 by major cell types; many results are included in code. 
- getMarkerFeatures groupBy sample. 
- Loops through each sample's markers.
    - File ex: microglia_markers_bySample_2025-01-28.csv
- Prints total number of returns by sample and treatment group. 
    - Code includes results for each major cell type. 
    - File ex: endo_bySampleTx_withGene-Cell-Counts_2025-01-28.csv & endo_byTx_withGene-Cell-Counts_2025-01-28.csv
- Identifies samples with no marker returns by major cell type. 
- LASSO ran for all samples/treatment groups
    - Load file ex: microglia_markers_bySample_2025-01-28.csv & TaskA-45_AllStats.csv
    - cv.glmnet used for LASSO regression cross-validation to select the best lambda using alpha = 1 (alpha = 1 for Lasso (Ridge would be 0)) and family = "mgaussian" (changed to "mgaussian" for multivariate Gaussian). 
    - Plots LASSO path, gets predictions for behavioral responses, and prints correlations. 
- Loop runs through all behavioral tests and cell types. 
    - File ex: Lasso_Correlations/Correlations_", sub(".csv
- LASSO compares differences between treatment groups; recommended to include if you want a global model that adjusts for treatment. 
    - Does not run by treatment group. 
- Compares prediction accuracy across treatments; recommended to run if you want to see if the model works better for some treatments. 
    - Does not run by treatment group. 
- Computes correlations seperately for each treatment group. 
- Runs scatterplots, Pearson's R, and Spearman correlations by treatment group for each major cell type. 
    - File ex: microglia_Task45_correlationPlot_2025-01-30.pdf & microglia_Task45_correlationResults_2025-01-30.csv
- Loop runs code above using AllStats files. 
    - Plot file ex: cell_type, "_", test_name, "_AllScores_correlationPlot_2025-01-30.pdf
    - File ex: cell_type, "_", test_name, "_AllScores_correlationResults_2025-01-30.csv"
- Loop runs above code using WrongStats files. 
    - Plot file ex: cell_type, "_", test_name, "_WrongScores_correlationPlot_2025-01-30.pdf"
    - File ex: cell_type, "_", test_name, "_WrongScores_correlationResults_2025-01-30.csv"
- Includes additional code that was not used.
## 20_Lasso_2025-01-25.R was an earlier version of this code. 
## 20_Lasso_2025-01-27.R was an earlier version of this code. 

## 20_Lasso_Clusters_2025-02-04.R
- Completed for clusters 3, 8, 11, 14, 18, 22 & behavioral tests TaskA-45, TaskB-59, TaskC-48, TaskD-149
- Load file ex: C18_TaskD-149_StatsCombined_2025-02-03.csv
- Lasso regression performed for each behavioral metric.
- cv.glmnet fits LASSO with cross-validation, using alpha = 1, nfolds = 5.
- glmnet fits the model with optimal lambda, using alpha = 1, lambda = best_lambda.
- Get coefficients (excluding intercept) and calculates R-squared for the model. 
- Loops through all behavioral metrics and plots combined LASSO results.
    - File ex: C18-TaskD-149_lasso_analysis_2025-02-03.pdf
    - Saves a summary table, file ex: C18_TaskD-149_LASSO_Summary_2025-02-04.csv
- Loop runs through each "Wrong" behavioral result from TaskA-45, TaskB-59, TaskC-48, TaskD-149.
    - Load data ex: C18_TaskC-48_WrongStatsCombined_2025-02-03.csv
    - File ex: C18-TaskC-48-Wrong_lasso_analysis_2025-02-03.pdf
    - Summary file ex: C18_TaskC-48-Wrong_LASSO_Summary_2025-02-04.csv
    
## 20_Lasso_wFDR_2025-02-14.R
- Completed for each major cell type (astro, gaba, glut, oligo, endo, microglia) and behavioral test (TaskA-45, TaskB-59, TaskC-48, TaskD-149). 
- Load file ex: astro_TaskA-45_StatsCombined_2025-02-03.csv
    - cv.glmnet fits LASSO with cross-validation, using alpha = 1, nfolds = 5.
    - glmnet fits the model with optimal lambda, using alpha = 1, lambda = best_lambda.
    - Get coefficients (excluding intercept) and calculates R-squared for the model. 
    - Calculates p-values for LASSO coefficients using bootstrapping. 
    - Standardizes predictors; cv.glmnet calculates original fit; coef calculates original coefficients; LASSO fit to bootstrap sample; p-values calculated. 
- Includes an enhanced LASSO function with p-values that applies FDR correction across all results using p.adjust method = "BH". 
    - Enhanced plotting function includes FDR results. 
    - Results table, file: lasso_analysis_results_with_fdr.csv
- Loop to run through all files with the same file name ending. 
    - Requires that previous LASSO functions are loaded. 
    - Performs LASSO with cross validation, R-squared calculated as above, and bootstrapping, but includes modified p-value calculations for LASSO coefficients. 
    - P-value calculations modified to be less stringent for zero coefficients. 
- Loop provides modified FDR correction function with adjustable thresholds. 
- Loop runs modified main analysis using less astringent parameters. 
    - File ex: pattern = "*_StatsCombined_2025-02-03.csv"
    

###############################################
###############################################
###############################################


# 21. GO & KEGG Gene-Behavior Correlation Analysis

## 21_Clusters-CorrelationGenes_GO-KEGG_2025-02-13.R
- Loads file ex: C3_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv
    - Filters sig genes and converts gene symbols to Entrez IDs.
    - enrichGO run for BP, MF, and CC; pvalueCutoff = 0.05, qvalueCutoff = 0.05.
    - enrichKEGG ran; pvalueCutoff = 0.05, qvalueCutoff = 0.05; prints a summary of result GO terms and KEGG pathways. 
    - File name ex: GO_biological_process_results.csv, GO_molecular_function_results.csv, GO_cellular_component_results.csv, KEGG_pathway_results.csv
    - Plot file name ex: GO_biological_process_plot.pdf, GO_molecular_function_plot.pdf, GO_cellular_component_plot.pdf, KEGG_pathway_plot.pdf
    - Dotplot for KEGG sig pathways, plot file name: KEGG_pathway_dotplot.pdf
- Loop for all files with the same file ending in the directory
    - enrichGO BP, MF, and CC ran with pvalueCutoff = 0.05, qvalueCutoff = 0.05.
    - KEGG enrichment set high for testing at pvalueCutoff = 1, qvalueCutoff = 1; then ran with pvalueCutoff = 0.05, qvalueCutoff = 0.05.
    - GO results file name ex: dir_name, "KEGG_pathway_results.csv", dir_name, "KEGG_pathway_results_all.csv", dir_name, "GO_biological_process_results.csv", dir_name, "GO_molecular_function_results.csv", and dir_name, "GO_cellular_component_results.csv"
    - KEGG results file name ex: dir_name, "KEGG_pathway_results.csv"
    - GO result plots file name ex: GO_biological_process_plot.pdf, GO_molecular_function_plot.pdf, GO_cellular_component_plot.pdf
    - KEGG result plots file name ex: dir_name, "KEGG_pathway_barplot.pdf"
    - Saves the gene ID mapping for reference, file: dir_name, "gene_id_mapping.csv"
    - Creates and saves a summary for GO and KEGG results, file: pattern = "*CorrelationStats.*withFDR_2025-02-12\\.csv
    
## 21_CellType-CorrelationGenes_GO-KEGG_2025-02-13.R
- Runs same code as 21_Clusters-CorrelationGenes_GO-KEGG_2025-02-13.R but for cell types. 
- Loads file ex: gaba_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv
    

###############################################
###############################################
###############################################


Code is adapted from: Granja JM, Corces MR et al., ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nature Genetics (2021)
