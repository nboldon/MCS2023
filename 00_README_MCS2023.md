# Project TIN CANN
## By: Naomi Boldon

# MCS2023


##################################################
##################################################
##################################################


## 01. The following software and package versions were used to complete the project:

##################################################

01_LICENSE
Specifies permissions granted and copyright notice. 

01_Libraries.R
Catelog of libraries used for the project. 

01_Conda_Update_Dec2023
Updated Beartooth conda environment for interactive session

01_beartooth_environment_2024-Feb.yml
Beartooth conda computing environment

01_beartooth_environment_2024-Oct
Beartooth conda computing environment

01_local_environment_2025-Mar
Naomi's personal computer computing environment


##################################################
##################################################
##################################################


## The following code scripts are used to complete the project:

##################################################


## 02. Demultiplex code

- Used to convert custom scATAC-seq library prep fragments to workable Cellranger files

Custom code from Qi Sun at Cornell University


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

## 03. Project MCS1

- Creates arrow files and project MCS1 for quality control and downstream analysis.

03_projMCS1_ArrowFiles.R
Creates arrow files in ArchR to create project

03_projMCS1.R
Adds doublet scores
Creates dataframes for nFrags and TSS Enrichment
Plots QC scores 
  Density plot for nFrags vs TSS Enrichment
  Ridge and violin plots for nFrags and TSS Enrichment individually plotted)

QC files:
03_TSS-vs-Frags.pdf
03_TSS7_QC-MCS1.pdf
03_QC_FragSize-Distro_2024-02-29.pdf
03_Peak-Call-Summary.pdf


#############################################
#############################################
#############################################


# 04. Project MCS2

- Filters doublets, runs initial DA gene accessibility groupBy cluster (no subset), and creates initial visualizations for all single cells in the project, including UMAPs, t-SNE, heatmaps, browser tracks. 

04_projMCS2.R
Filters doublets
Adds Iterative LSI
Adds Harmony
Adds Clusters
Creates cluster confusion matrix


04_projMCS2_scEmbeddings.R
UMAP using LSI by sample & cluster
TSNE using LSI by sample & cluster
UMAP using Harmony by sample & cluster
TSNE using Harmony by sample & cluster

04_Files created from projMCS2_scEmbedding.R
    04_UMAP-Sample-Clusters.pdf
    04_TSNE-Sample-Clusters.pdf
    04_UMAP2Harmony-Sample-Clusters.pdf
    04_TSNE2Harmony-Sample-Clusters.pdf


04_projMCS2_MarkerGenes.R
GroupBy: Cluster (no subset), FDR <= 0.01, Log2FC >= 1.25
getMarkerFeatures
Marker list by cluster (FDR & Log2FC)
Heatmaps for marker features; cowplots for all genes
NOTE: projMCS2_MarkerGenes.R did not calculate abs(Log2FC)

## projMCS5_MarkerGenes.R 
Group by: Cluster (no subset), FDR <= 0.01, abs(Log2FC) >= 1.25
    Files saved for each cluster; ex: C1_MarkerGenes_2025-03-06.csv

Summary files created from analysis results of projMCS2_MarkerGenes.R
    04_Cluster_Analysis_2024-02-12.xlsm
        Each tab summarizes cell count information by cluster, including:
        cell counts by sample and cluster (incl. median TSS & nFrags per cell)
        Normalized cell abundance percentages by sample (cells in cluster by sample divided by total cells in  sample; incl. total cell counts by sample)
      Normalized cell abundance by treatment group (incl. total cell counts by cluster and treatment group)
      Cell counts by cluster and genotype
      Cell counts by cluster and sequencing lane
      Cell counts by PCR plate (library preparation day)
    04_GeneMarkers_byCluster_1-31-2024.xlsx
        Includes cell type-specific marker gene notes from literature reviews in tabs 1 & 2
        There is a tab for each cluster that highlights:
            The top 50 marker genes, based on FDR 
            Cell type-specific marker genes found in literature reviews
        There are two additional tabs for C1, C3, C8, C18, and C21 that highlight:
            One tab specifies the Panther GO family, molecular function, biological process, cellular component, protein class, and pathways affected by specified genes. 
            The other tab displays seperate bargraphs for Panther GO categories associated with molecular function, biological process, cellular component, protein class, and pathways affected by specified genes.
    04_MCS2023_1-31-2024.pptx
        Summarizes key points of interest from the above analysis

Heatmaps of projMCS2 by cluster from projMCS2_MarkerGenes.R
04_UMAP-MarkerGenes-WO-Imputation.pdf
04_GeneScores-Marker-Heatmap.pdf
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
04_GeneMarkers_Top50_Heatmap_1-25-2024.pdf


##############################################
##############################################
##############################################


# 05. Project MCS3

## projMCS3.R

Defining cluster ID with scRNA-seq 

Unable to complete; no RNA-seq data for cohort; study limitation


###############################################
###############################################
###############################################


# 06. Project MCS4 

## projMCS4.R

Add pseudobulk replicates & reproduceable peak set using MACS2


###############################################
###############################################
###############################################


# 07. Project MCS5

## projMCS5.R
- Adds peak matrix
- Identifies marker peaks by cluster
  - cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
- Plots marker peaks
  - cutOff = "FDR<=0.01 & Log2FC>=1.25"
- Adds motif annotations
  - motifSet = "cisbp"
NOTE: The 07_projMCS5.R code did not use abs(Log2FC)
The 07_projMCS5-w-ABS-Log2FC.R code accommodates for that and reruns getMarkerFeatures before applying cutoffs

## 07_projMCS5-w-ABS-Log2FC.R
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

## Marker genes that were subset by treatment group, then groupBy cluster
## were combined into one spreadsheet using: 
09_Combine_Spreadsheets.R
- File names ex: C1_TxMarkers_Combined.csv


###############################################
###############################################
###############################################


# 10. Marker genes subset by cluster, grouped by treatment


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
 
10_projMCS6_GeneMarkers_TxComp-byCluster_2024-06-26.R
- Adds Harmony, subsets by cluster, groupBy treatment, uses MAGIC
- cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5
- Ex. file name: C25_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv
- Subset by cluster, group by sample: all returns were NULL
- Uses top 200 marker genes to create heatmaps & z-scores 
- Heatmap cutOff = "FDR <= 1 & Log2FC >=0.01"
- Ex. file names saved: C4_TxComp_Heatmap_2024-03-21.pdf & C4_byTx_zscores_2024-03-21.csv



###############################################
###############################################
###############################################


# 11. Marker genes subset by cell type, grouped by treatment

11_projMCS7_GeneMarkers_byCellType-TxGrp_2025-01-16.R
- Similar code to projMCS7_GeneMarkers_byCluster-TxGrp_2025-01-16.R but focuses on cell type groups
* - Files saved ex: astro_GeneMarker_List_2025-01-16.csv
 
11_CellType_TxComps.R
- Code did not produce meaningful results


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


## ProjMCS6 - Harmony Added; Single cell embeddings & tx comparisons

13_Clusters_TvsTComp_2024-06-28.R
- Adds Harmony
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

13_Q4_GenotypeTxComps_2024-06-28.R
- Adds Harmony
- Groups by genotype
- cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5"
- Files: GenotypeTxComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv & GenotypeCompByCluster_FDR-0-1_Log2FC-0-5_2024-06-28.csv

13_Q5_GenotypeClusterComps_2024_06-28.R
- Adds Harmony
- Subsets by genotype, groups by cluster
- Ex file names: C2_GenotypeComp_FDR-0-1_Log2FC-0-5_2024-06-28.csv


###############################################
###############################################
###############################################


# 15. ProjMCS7 - Peak Matrix creation and motif annotations

15_projMCS7.R

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


15_CA_Motif_Plot_byCluster.R
- Uses motifs.csv to create a heatmap of DA Motif Count Comparison Across Clusters.
- File: DA_Motif_Count_Comparison_2024-09-27.png


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
    
    
16_nFrags-TSS_byTx-Cluster.R
- Subset by cluster, groupBy treatment.
- Violin plot colorBy cell data using log10(nFrags).
    - File: C25_nFrags-byTx_2024-08-13
- Violin plot colorBy cell data using TSS Enrichment.
    - File: C25_TSS-byTx_2024-08-13
    

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



###############################################
###############################################
###############################################


## 17. ProjMCS9 Peak matrix pairwise comparisons

projMCS9.R
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


## 18. Cluster 18 Subset

18_C18_Subset.R
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

C18_Subset_TxComps_byCluster.R
- Harmony added, subset by subcluster, assigns treatment groups. 

C18_TxComps_byCluster.R
- Addresses Research Question 3
    - Harmony added, subset by subcluster, treatment group assigned, getMarkerFeatures groupBytreatment, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5").
    - File name ex: C18.18_TxComp_FDR-0-1_Log2FC-0-5_2024-07-25.csv
- Addresses Research Question 5
    - Harmony added, subset by subcluster, treatment group assigned, getMarkerFeatures by treatment group pairwise comparison, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5").
    - File name ex: C18.18_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-24.csv
- Reads pairwise comparison spreadsheets into a dataframe. 
    - File name ex: C18.18-TxComp_geneMarkersCombined_2024-07-25.csv
    

C18_T1-2-4v3_byCluster_2024-07-31.R
- Follows code workflow from C18_TxComps_byCluster pairwise comparison, but analyzes specific differences between combined groups T1, T2, and T4 compared to T3.


###############################################
###############################################
###############################################

# 14. Volcano Plots

Volcano_byCluster-Tx.R
- Add Harmony & impute weights
- Subset byy cluster; group by tx
- Pairwise comparisons
- CutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"

14_projMCS6_Volcano_2024.08-08.R
- Add Harmony & impute weights
- Subset byy cluster; group by tx
- Pairwise comparisons
- CutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"


###############################################
###############################################
###############################################


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
