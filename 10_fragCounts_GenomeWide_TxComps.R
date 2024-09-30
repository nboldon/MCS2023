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


########################################
########################################
########################################


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


########################################
########################################


library(data.table)

# Define the path to the fragment files
fragment_dir <- "/project/eon/nboldon/MCS2023/TSS7_FDR-0.01_Log2FC-1.25"

# Create a function to load fragments for a specific treatment
load_fragments_by_treatment <- function(samples, treatment_label) {
  all_fragments <- NULL
  for (sample in samples) {
    file_path <- file.path(fragment_dir, paste0(sample, "_fragments.tsv.gz"))
    if (file.exists(file_path)) {
      fragments <- fread(file_path, sep="\t", header=FALSE)
      colnames(fragments) <- c("chromosome", "start", "end", "barcode", "metadata")
      fragments$treatment <- treatment_label
      all_fragments <- rbind(all_fragments, fragments)
    } else {
      message(paste0("File not found for sample: ", sample))
    }
  }
  return(all_fragments)
}

# Load fragments for each treatment group 
fragments_t1 <- load_fragments_by_treatment(c("C302", "C306", "C309","C318", "C323", "C328", "C332", 
                                              "C337", "C339", "C346", "C351", "C353", "C360"), "t1")
fragments_t2 <- load_fragments_by_treatment(c("C304", "C308", "C312", "C349", "C315", "C321", "C324", 
                                              "C355", "C327", "C330", "C333", "C358", "C336", "C342", "C348", "C362"), "t2")
fragments_t3 <- load_fragments_by_treatment(c("C305", "C307", "C313", "C350", "C316", "C320", "C322", 
                                              "C352", "C325", "C334", "C359", "C340", "C341", "C345", "C364"), "t3")
fragments_t4 <- load_fragments_by_treatment(c("C301", "C303", "C310", "C314", "C319", "C335", "C338", 
                                              "C344", "C354", "C356", "C361", "C363"), "t4")


########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################


## Bin and count fragments by chromosome

bin_size <- 5000  # 5kb bins

# Add a bin column to each fragment dataset
fragments_t1[, bin := floor(start / bin_size)]
fragments_t2[, bin := floor(start / bin_size)]
fragments_t3[, bin := floor(start / bin_size)]
fragments_t4[, bin := floor(start / bin_size)]

# Count the number of fragments per bin for each treatment
counts_t1 <- fragments_t1[, .N, by = .(chromosome, bin)]
counts_t2 <- fragments_t2[, .N, by = .(chromosome, bin)]
counts_t3 <- fragments_t3[, .N, by = .(chromosome, bin)]
counts_t4 <- fragments_t4[, .N, by = .(chromosome, bin)]

# Merge the counts for the treatments
merged_counts <- Reduce(function(x, y) merge(x, y, by=c("chromosome", "bin"), all=TRUE),
                        list(counts_t1, counts_t2, counts_t3, counts_t4))
colnames(merged_counts) <- c("chromosome", "bin", "count_t1", "count_t2", "count_t3", "count_t4")

# Replace NA values with 0
merged_counts[is.na(merged_counts)] <- 0



## Subtract fragment counts between treatment groups

# Subtract counts between treatment groups
merged_counts[, difference_t1_t2 := count_t1 - count_t2]
merged_counts[, difference_t3_t4 := count_t3 - count_t4]
merged_counts[, difference_t1_t3 := count_t1 - count_t3]
merged_counts[, difference_t1_t4 := count_t1 - count_t4]
merged_counts[, difference_t2_t3 := count_t2 - count_t3]
merged_counts[, difference_t2_t4 := count_t2 - count_t4]



########################################
########################################
########################################
########################################
########################################



## Plot fragment count differences 

library(ggplot2)

# Plot difference between t1 and t2 on a specific chromosome, e.g., chromosome 16
ggplot(merged_counts[chromosome == "chr16"], aes(x=bin, y=difference_t1_t2)) +
  geom_line() +
  labs(title="Fragment Count Difference (t1 vs t2) on Chromosome 16",
       x="Genomic Bin (5kb)", y="Fragment Count Difference") +
  theme_minimal()

# Plot difference between t3 and t4 on a specific chromosome
ggplot(merged_counts[chromosome == "chr16"], aes(x=bin, y=difference_t3_t4)) +
  geom_line() +
  labs(title="Fragment Count Difference (t3 vs t4) on Chromosome 16",
       x="Genomic Bin (10kb)", y="Fragment Count Difference") +
  theme_minimal()



########################################
########################################



## Adjust for multiple chromosomes


chromosomes <- unique(merged_counts$chromosome)
plot_list <- list()

for (chr in chromosomes) {
  plot_list[[chr]] <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t1_t2)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t1 vs t2) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
}

# Display one of the plots
print(plot_list[["chr16"]])




########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################


### Loop for all treatment group comparisons using chrom. 16 & 17

# Define the chromosomes to analyze
chromosomes <- c("chr16", "chr17")

# Create empty lists to store the plots for each comparison
plot_list_t1_t2 <- list()
plot_list_t3_t4 <- list()
plot_list_t1_t3 <- list()
plot_list_t1_t4 <- list()
plot_list_t2_t3 <- list()
plot_list_t2_t4 <- list()

# Specify the directory where you want to save the plots
output_dir <- "/project/eon/nboldon/MCS2023/ProjMCS7"  

# Loop through the chromosomes and generate and save plots for both comparisons
for (chr in chromosomes) {
  
  # Plot for t1 vs t2 on each chromosome
  plot_t1_t2 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t1_t2)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t1 vs t2) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t1_t2[[chr]] <- plot_t1_t2
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t1_vs_t2_", chr, "_2024-09-25.png"), 
         plot = plot_t1_t2, width = 8, height = 6)
  
  
  
  # Plot for t3 vs t4 on each chromosome
  plot_t3_t4 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t3_t4)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t3 vs t4) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t3_t4[[chr]] <- plot_t3_t4
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t3_vs_t4_", chr, "2024-09-25.png"), 
         plot = plot_t3_t4, width = 8, height = 6)
  
  
  
  
  # Plot for t1 vs t3 on each chromosome
  plot_t1_t3 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t1_t3)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t1 vs t3) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t1_t3[[chr]] <- plot_t1_t3
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t1_vs_t3_", chr, "_2024-09-25.png"), 
         plot = plot_t1_t3, width = 8, height = 6)
  
  
  
  # Plot for t1 vs t4 on each chromosome
  plot_t1_t4 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t1_t4)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t1 vs t4) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t1_t4[[chr]] <- plot_t1_t4
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t1_vs_t4_", chr, "_2024-09-25.png"), 
         plot = plot_t1_t4, width = 8, height = 6)
  
  
  
  # Plot for t2 vs t3 on each chromosome
  plot_t2_t3 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t2_t3)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t2 vs t3) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t2_t3[[chr]] <- plot_t2_t3
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t2_vs_t3_", chr, "_2024-09-25.png"), 
         plot = plot_t2_t3, width = 8, height = 6)
  
  
  
  # Plot for t2 vs t4 on each chromosome
  plot_t2_t4 <- ggplot(merged_counts[chromosome == chr], aes(x=bin, y=difference_t2_t4)) +
    geom_line() +
    labs(title=paste("Fragment Count Difference (t2 vs t4) on", chr),
         x="Genomic Bin (5kb)", y="Fragment Count Difference") +
    theme_minimal()
  
  # Store the plot in the list
  plot_list_t2_t4[[chr]] <- plot_t2_t4
  
  # Save the plot as a PNG file
  ggsave(filename = paste0(output_dir, "Fragment_Count_Difference_t2_vs_t4_", chr, "_2024-09-25.png"), 
         plot = plot_t2_t4, width = 8, height = 6)
  

}

# Display the plots for chromosome 16 and 17 for t1 vs t2
print(plot_list_t1_t2[["chr16"]])
print(plot_list_t1_t2[["chr17"]])

# Display the plots for chromosome 16 and 17 for t3 vs t4
print(plot_list_t3_t4[["chr16"]])
print(plot_list_t3_t4[["chr17"]])



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
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################



## To plot frag counts of chromosome/treatment by rank order in bins


library(data.table)
library(ggplot2)

# Load the fragments for each treatment (as in your original code)
# Assuming `fragments_t1`, `fragments_t2`, `fragments_t3`, `fragments_t4` are already loaded

# Set bin size
bin_size <- 5000  # 5kb bins

# Add a bin column to each fragment dataset
fragments_t1[, bin := floor(start / bin_size)]
fragments_t2[, bin := floor(start / bin_size)]
fragments_t3[, bin := floor(start / bin_size)]
fragments_t4[, bin := floor(start / bin_size)]

# Count the number of fragments per bin for each treatment
counts_t1 <- fragments_t1[, .N, by = .(chromosome, bin)]
counts_t2 <- fragments_t2[, .N, by = .(chromosome, bin)]
counts_t3 <- fragments_t3[, .N, by = .(chromosome, bin)]
counts_t4 <- fragments_t4[, .N, by = .(chromosome, bin)]

# Merge the counts for the treatments into one table
merged_counts <- Reduce(function(x, y) merge(x, y, by=c("chromosome", "bin"), all=TRUE),
                        list(counts_t1, counts_t2, counts_t3, counts_t4))
colnames(merged_counts) <- c("chromosome", "bin", "count_t1", "count_t2", "count_t3", "count_t4")

# Replace NA values with 0
merged_counts[is.na(merged_counts)] <- 0

# Calculate ranks within each chromosome for each treatment
merged_counts[, rank_t1 := rank(-count_t1), by = chromosome]  # Negative rank to rank in descending order
merged_counts[, rank_t2 := rank(-count_t2), by = chromosome]
merged_counts[, rank_t3 := rank(-count_t3), by = chromosome]
merged_counts[, rank_t4 := rank(-count_t4), by = chromosome]

# Reshape data to long format for easier plotting
rank_data <- melt(merged_counts, 
                  id.vars = c("chromosome", "bin"), 
                  measure.vars = c("rank_t1", "rank_t2", "rank_t3", "rank_t4"), 
                  variable.name = "treatment", 
                  value.name = "rank")

# Replace the treatment variable values for cleaner labels
rank_data[, treatment := factor(treatment, levels = c("rank_t1", "rank_t2", "rank_t3", "rank_t4"), 
                                labels = c("t1", "t2", "t3", "t4"))]

# Plot the rank distributions by chromosome and treatment
ggplot(rank_data, aes(x = chromosome, y = rank, color = treatment)) +
  geom_point(alpha = 0.5) +
  labs(title = "Rank Distribution of Fragment Counts by Chromosome and Treatment", 
       x = "Chromosome", 
       y = "Rank of Fragment Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("t1" = "blue", "t2" = "green", "t3" = "red", "t4" = "purple"))



################################################



## Save plots by chromosome

plot_rank_distribution_per_chromosome <- function(data, treatment) {
  chromosomes <- unique(data$chromosome)
  
  for (chr in chromosomes) {
    # Filter data for the current chromosome
    data_chr <- data[data$chromosome == chr]
    
    # Create the plot
    p <- ggplot(data_chr, aes(x = bin, y = N)) +
      geom_line(color = "blue") +
      labs(title = paste("Rank Distribution of Fragment Counts -", treatment, "Chromosome", chr),
           x = "Bins", y = "Fragment Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Save the plot using ggsave with a chromosome-specific file name
    ggsave(paste0("rank_distribution_", treatment, "_", chr, "_2024-09-25.png"), plot = p, width = 10, height = 6)
  }
}

# Generate and save plots for each chromosome in each treatment
plot_rank_distribution_per_chromosome(counts_t1, "t1")
plot_rank_distribution_per_chromosome(counts_t2, "t2")
plot_rank_distribution_per_chromosome(counts_t3, "t3")
plot_rank_distribution_per_chromosome(counts_t4, "t4")



plot_rank_distribution <- function(data, treatment, output_path) {
  ggplot(data, aes(x = bin, y = N, group = chromosome, color = chromosome)) +
    geom_line() +
    labs(title = paste("Rank Distribution of Fragment Counts -", treatment),
         x = "Bins", y = "Fragment Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~chromosome, scales = "free_x")  # Separate plot per chromosome with free x-axis scaling
  
  # Save the plot using ggsave
  ggsave(output_path, width = 12, height = 8)
}

# Plot and save rank distributions for each treatment
plot_rank_distribution(counts_t1, "t1", "rank_distribution_t1_2024-09-25.png")
plot_rank_distribution(counts_t2, "t2", "rank_distribution_t2_2024-09-25.png")
plot_rank_distribution(counts_t3, "t3", "rank_distribution_t3_2024-09-25.png")
plot_rank_distribution(counts_t4, "t4", "rank_distribution_t4_2024-09-25.png")

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
 [1] RColorBrewer_1.1-3          cowplot_1.1.2              
 [3] tidyr_1.3.1                 BiocManager_1.30.22        
 [5] Signac_1.12.0               Seurat_5.0.1               
 [7] SeuratObject_5.0.1          sp_2.1-2                   
 [9] pheatmap_1.0.12             enrichplot_1.20.3          
[11] clusterProfiler_4.8.3       org.Mm.eg.db_3.17.0        
[13] AnnotationDbi_1.64.1        rhdf5_2.44.0               
[15] SummarizedExperiment_1.32.0 Biobase_2.62.0             
[17] MatrixGenerics_1.14.0       Rcpp_1.0.12                
[19] Matrix_1.6-4                GenomicRanges_1.54.1       
[21] GenomeInfoDb_1.38.5         IRanges_2.36.0             
[23] S4Vectors_0.40.2            BiocGenerics_0.48.1        
[25] matrixStats_1.2.0           data.table_1.14.10         
[27] stringr_1.5.1               plyr_1.8.9                 
[29] magrittr_2.0.3              ggplot2_3.4.4              
[31] gtable_0.3.4                gtools_3.9.5               
[33] gridExtra_2.3               ArchR_1.0.2                

loaded via a namespace (and not attached):
  [1] fs_1.6.3                           spatstat.sparse_3.0-3             
  [3] bitops_1.0-7                       HDO.db_0.99.1                     
  [5] httr_1.4.7                         tools_4.3.2                       
  [7] sctransform_0.4.1                  utf8_1.2.4                        
  [9] R6_2.5.1                           lazyeval_0.2.2                    
 [11] uwot_0.1.16                        rhdf5filters_1.12.1               
 [13] withr_3.0.0                        progressr_0.14.0                  
 [15] textshaping_0.3.7                  cli_3.6.2                         
 [17] Cairo_1.6-2                        spatstat.explore_3.2-5            
 [19] fastDummies_1.7.3                  scatterpie_0.2.1                  
 [21] labeling_0.4.3                     spatstat.data_3.0-3               
 [23] ggridges_0.5.4                     pbapply_1.7-2                     
 [25] systemfonts_1.0.5                  Rsamtools_2.18.0                  
 [27] yulab.utils_0.1.0                  gson_0.1.0                        
 [29] DOSE_3.26.2                        R.utils_2.12.3                    
 [31] parallelly_1.36.0                  BSgenome_1.70.1                   
 [33] RSQLite_2.3.5                      generics_0.1.3                    
 [35] gridGraphics_0.5-1                 BiocIO_1.12.0                     
 [37] ica_1.0-3                          spatstat.random_3.2-1             
 [39] dplyr_1.1.4                        GO.db_3.18.0                      
 [41] fansi_1.0.6                        abind_1.4-5                       
 [43] R.methodsS3_1.8.2                  lifecycle_1.0.4                   
 [45] yaml_2.3.8                         qvalue_2.32.0                     
 [47] SparseArray_1.2.3                  Rtsne_0.17                        
 [49] blob_1.2.4                         promises_1.2.1                    
 [51] crayon_1.5.2                       miniUI_0.1.1.1                    
 [53] lattice_0.22-5                     KEGGREST_1.42.0                   
 [55] pillar_1.9.0                       fgsea_1.26.0                      
 [57] rjson_0.2.21                       future.apply_1.11.0               
 [59] codetools_0.2-19                   fastmatch_1.1-4                   
 [61] leiden_0.4.3.1                     glue_1.7.0                        
 [63] downloader_0.4                     ggfun_0.1.3                       
 [65] vctrs_0.6.5                        png_0.1-8                         
 [67] treeio_1.24.3                      spam_2.10-0                       
 [69] cachem_1.0.8                       S4Arrays_1.2.0                    
 [71] mime_0.12                          tidygraph_1.2.3                   
 [73] survival_3.5-7                     RcppRoll_0.3.0                    
 [75] ellipsis_0.3.2                     fitdistrplus_1.1-11               
 [77] ROCR_1.0-11                        nlme_3.1-164                      
 [79] ggtree_3.8.2                       bit64_4.0.5                       
 [81] RcppAnnoy_0.0.21                   irlba_2.3.5.1                     
 [83] KernSmooth_2.23-22                 colorspace_2.1-0                  
 [85] DBI_1.2.1                          tidyselect_1.2.0                  
 [87] bit_4.0.5                          compiler_4.3.2                    
 [89] BSgenome.Mmusculus.UCSC.mm10_1.4.3 DelayedArray_0.28.0               
 [91] plotly_4.10.3                      shadowtext_0.1.2                  
 [93] rtracklayer_1.62.0                 scales_1.3.0                      
 [95] lmtest_0.9-40                      digest_0.6.33                     
 [97] goftest_1.2-3                      spatstat.utils_3.0-4              
 [99] XVector_0.42.0                     htmltools_0.5.7                   
[101] pkgconfig_2.0.3                    fastmap_1.1.1                     
[103] rlang_1.1.3                        htmlwidgets_1.6.4                 
[105] shiny_1.8.0                        farver_2.1.1                      
[107] zoo_1.8-12                         jsonlite_1.8.8                    
[109] BiocParallel_1.36.0                R.oo_1.25.0                       
[111] GOSemSim_2.26.1                    RCurl_1.98-1.14                   
[113] GenomeInfoDbData_1.2.11            ggplotify_0.1.2                   
[115] dotCall64_1.1-1                    patchwork_1.1.3                   
[117] Rhdf5lib_1.22.1                    munsell_0.5.0                     
[119] ape_5.7-1                          viridis_0.6.4                     
[121] reticulate_1.34.0                  stringi_1.8.3                     
[123] ggraph_2.1.0                       zlibbioc_1.48.0                   
[125] MASS_7.3-60                        listenv_0.9.0                     
[127] ggrepel_0.9.4                      deldir_1.0-9                      
[129] Biostrings_2.70.1                  graphlayouts_1.0.2                
[131] splines_4.3.2                      tensor_1.5                        
[133] igraph_1.5.1                       spatstat.geom_3.2-7               
[135] RcppHNSW_0.5.0                     reshape2_1.4.4                    
[137] XML_3.99-0.16.1                    tweenr_2.0.2                      
[139] httpuv_1.6.13                      RANN_2.6.1                        
[141] purrr_1.0.2                        polyclip_1.10-6                   
[143] future_1.33.0                      scattermore_1.2                   
[145] ggforce_0.4.1                      xtable_1.8-4                      
[147] restfulr_0.0.15                    RSpectra_0.16-1                   
[149] tidytree_0.4.5                     later_1.3.2                       
[151] ragg_1.2.7                         viridisLite_0.4.2                 
[153] tibble_3.2.1                       aplot_0.2.2                       
[155] memoise_2.0.1                      GenomicAlignments_1.38.2          
[157] cluster_2.1.4                      globals_0.16.2  