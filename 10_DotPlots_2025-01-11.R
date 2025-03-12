
##################################################

#Load libraries
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

setwd("/Volumes/DataBox/ProjMCS6/Dotplots")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

##################################################

projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)
projMCS6
getAvailableMatrices(projMCS6)


##################################################
##################################################
#########################################
#########################################
#########################################
#########################################
#########################################


#To get a dataframe of the mean and relative accessibility of each gene for each cluster:

Cluster.markers <- getMarkerFeatures(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

assays(Cluster.markers)

#Take a look at how this is encoded
Cluster.markers@assays@data$Mean[1:10,1:3]
Cluster.markers@assays@data$MeanDiff[1:10,1:3]
Cluster.markers@assays@data$FDR>0

# Assuming Cluster.markers is already created

library(tidyverse)
library(reshape2)
library(ggdendro)
library(cowplot)
library(patchwork)
library(viridis)

# Data Preparation Steps
mean_acc <- Cluster.markers@assays@data$Mean
rel_acc  <- Cluster.markers@assays@data$MeanDiff
FDR_acc <- Cluster.markers@assays@data$FDR
Log2FC_acc <- Cluster.markers@assays@data$Log2FC

# Add rownames and column names
rownames(mean_acc) <- rownames(rel_acc) <- rownames(FDR_acc) <- rownames(Log2FC_acc) <- rowData(Cluster.markers)$name
colnames(mean_acc) <- colnames(rel_acc) <- colnames(FDR_acc) <- colnames(Log2FC_acc) <- colnames(Cluster.markers)

# Melt the data
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc), varnames = c("gene", "cluster"), value.name = "mean")
melted_rel_acc <- reshape2::melt(as.matrix(rel_acc), varnames = c("gene", "cluster"), value.name = "rel")
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc), varnames = c("gene", "cluster"), value.name = "FDR")
melted_Log2FC_acc <- reshape2::melt(as.matrix(Log2FC_acc), varnames = c("gene", "cluster"), value.name = "Log2FC")

# Combine all data
acc_comb <- melted_mean_acc %>%
  left_join(melted_rel_acc, by = c("gene", "cluster")) %>%
  left_join(melted_FDR_acc, by = c("gene", "cluster")) %>%
  left_join(melted_Log2FC_acc, by = c("gene", "cluster"))

# Define genes of interest
interesting_markerGenes <- c("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin")

# Filter data for genes of interest
acc_filtered <- acc_comb %>% 
  filter(gene %in% interesting_markerGenes)

unique(acc_filtered$gene)

# Prepare matrix for clustering
mat <- acc_filtered %>%
  select(gene, cluster, rel) %>%
  pivot_wider(names_from = cluster, values_from = rel, values_fill = 0) %>%
  column_to_rownames("gene")

# Perform clustering
clust <- hclust(dist(as.matrix(mat)))

# Get gene order from clustering
gene_order <- clust$labels[clust$order]

# Create dendrogram data
ddata <- dendro_data(clust, type = "rectangle")

# Create dendrogram plot
dendro <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse() +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

# Create gene labels plot
gene_labels <- ggplot(data.frame(y = seq_along(gene_order), label = gene_order), aes(x = 1, y = y, label = label)) +
  geom_text(hjust = 1, size = 3) +
  scale_y_continuous(breaks = seq_along(gene_order), expand = c(0.005, 0.005)) +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

# Modify dotplot
dotplot <- acc_filtered %>%
  filter(mean > 0, abs(rel) > 0.1) %>%
  ggplot(aes(x = cluster, y = factor(gene, levels = rev(gene_order)), color = rel, size = mean)) +
  geom_point() +
  scale_color_gradientn(colours = viridis(20), limits = c(-2, 2),
                        oob = scales::squish, name = 'MeanDiff') +
  scale_size_continuous(name = "Mean", range = c(1, 6)) +
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

# Combine plots
combined_plot <- (dendro + gene_labels + dotplot) +
  plot_layout(widths = c(1, 1, 4)) +
  plot_annotation(title = "Gene Expression Dotplot with Dendrogram",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# Display combined plot
print(combined_plot)

# Save plot
ggsave("gene_expression_dotplot_with_dendrogram.png", combined_plot, width = 12, height = 8, dpi = 300)


#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################
#########################################



## Log2FC Dotplot


library(tidyverse)
library(reshape2)
library(ggdendro)
library(cowplot)
library(patchwork)
library(viridis)

# Identify marker features

Cluster.markers <- getMarkerFeatures(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

assays(Cluster.markers)

# Data Preparation Steps
mean_acc <- Cluster.markers@assays@data$Mean
Log2FC_acc <- Cluster.markers@assays@data$Log2FC
FDR_acc <- Cluster.markers@assays@data$FDR

# Add rownames and column names
rownames(mean_acc) <- rownames(Log2FC_acc) <- rownames(FDR_acc) <- rowData(Cluster.markers)$name
colnames(mean_acc) <- colnames(Log2FC_acc) <- colnames(FDR_acc) <- colnames(Cluster.markers)

# Melt the data
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc), varnames = c("gene", "cluster"), value.name = "mean")
melted_Log2FC_acc <- reshape2::melt(as.matrix(Log2FC_acc), varnames = c("gene", "cluster"), value.name = "Log2FC")
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc), varnames = c("gene", "cluster"), value.name = "FDR")

# Combine all data
acc_comb <- melted_mean_acc %>%
  left_join(melted_Log2FC_acc, by = c("gene", "cluster")) %>%
  left_join(melted_FDR_acc, by = c("gene", "cluster"))

# Define genes of interest
glut_markerGenes <- c("Cux2","Neurod1", "Neurod2", "Neurod6", "Satb2", "Tbr1")
gaba_markerGenes <- c("Drd2", "Gad1", "Pvalb", "Sp8", "Sst", "Vip")
microglia_markerGenes <- c("C1qa", "Ccl4", "Cd14", "Cx3cr1", "Csf1r", "Trem2")
oligo_markerGenes <- c("Mag", "Mbp", "Olig1", "Opalin", "Pdgfra", "Ugt8a")
endo_markerGenes <- c("Anpep", "Cldn5", "Flt1", "Fn1", "Itm2a", "Slco1a4")
astro_markerGenes <- c("Aqp4", "Cbs", "Dio2", "Gfap", "Mlc1", "Ppp1r3c")

##############################################################


## Run the following for each set of genes of interest

# Filter data for genes of interest
acc_filtered <- acc_comb %>% 
  filter(gene %in% astro_markerGenes)

# Prepare matrix for clustering
mat <- acc_filtered %>%
  select(gene, cluster, Log2FC) %>%
  pivot_wider(names_from = cluster, values_from = Log2FC, values_fill = 0) %>%
  column_to_rownames("gene")

# Perform clustering
clust <- hclust(dist(as.matrix(mat)))

# Get gene order from clustering
gene_order <- clust$labels[clust$order]

# Create dendrogram data
ddata <- dendro_data(clust, type = "rectangle")

# Create dendrogram plot
dendro <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse() +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

# Create gene labels plot
gene_labels <- ggplot(data.frame(y = seq_along(gene_order), label = gene_order), aes(x = 1, y = y, label = label)) +
  geom_text(hjust = 1, size = 3) +
  scale_y_continuous(breaks = seq_along(gene_order), expand = c(0.005, 0.005)) +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



####################################################

# Create dotplot
#dotplot <- acc_filtered %>%
#  filter(mean > 0, abs(Log2FC) > 0.1) %>%  # Adjust this filter as needed
#  ggplot(aes(x = cluster, y = factor(gene, levels = rev(gene_order)), color = Log2FC, size = mean)) +
#  geom_point() +
#  scale_color_gradientn(colours = viridis(20), limits = c(-2, 2), # Adjust the limits to focus on the range of interest
#                        oob = scales::squish, name = 'Log2FC') +
#  scale_size_continuous(name = "Mean", range = c(1, 6)) +
#  cowplot::theme_cowplot() +
#  theme(axis.line = element_blank(),
#        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#        axis.text.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        axis.title.y = element_blank(),
#        plot.margin = unit(c(0,0,0,0), "cm"))

# "astro_dotplot_with_dendrogram_log2fc_2025-01-11.png"


#####################################################
#####################################################
#####################################################


## Updated plot code to minimize less significant findings and emphasize more impactful findings

# Create dotplot
dotplot <- acc_filtered %>%
  filter(mean > 0.1, abs(Log2FC) > 0.1) %>%  # Adjust this filter as needed
  ggplot(aes(x = cluster, y = factor(gene, levels = rev(gene_order)), color = Log2FC, size = mean)) +
  geom_point() +
  scale_color_gradientn(colours = viridis(20), limits = c(-2.5, 2.5), # Adjust the limits to focus on the range of interest
                        oob = scales::squish, name = 'Log2FC') +
  scale_size_continuous(name = "Mean", range = c(1, 6)) +
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))


##########

# Combine plots
combined_plot <- (dendro + gene_labels + dotplot) +
  plot_layout(widths = c(1, 1, 4)) +
  plot_annotation(title = "Astrocyte Gene Dotplot with Dendrogram",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# Display combined plot
print(combined_plot)

# Save plot
ggsave("astrocyte_dotplot_with_dendrogram_2025-01-15.png", combined_plot, width = 12, height = 8, dpi = 300)



























#####################################################
#####################################################
#####################################################


## Updated plot code to minimize less significant findings and emphasize more impactful findings, dendrogram removed

# Create dotplot
dotplot <- acc_filtered %>%
  filter(mean > 0.1, abs(Log2FC) > 0.1) %>%
  ggplot(aes(x = cluster, y = factor(gene, levels = rev(gene_order)), color = Log2FC, size = mean)) + 
  geom_point() +
  scale_color_gradientn(colours = viridis(20), limits = c(-2.5, 2.5), 
                        oob = scales::squish, name = 'Log2FC') +
  scale_size_continuous(name = "Mean", range = c(1, 6)) +
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),  # Adjust font size for clarity
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing = unit(0.5, "lines"),  # Adjust space between rows
        strip.text.y = element_blank())  # Remove facet strip text if used

##########

# Create gene_labels plot (make sure it doesn't include redundant gene names)
gene_labels <- ggplot(acc_filtered, aes(x = cluster, y = factor(gene, levels = rev(gene_order)))) +
  geom_blank() +  # This removes any redundant gene names in the gene_labels plot
  theme_void()

##########

# Combine plots
combined_plot <- (gene_labels + dotplot) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(title = "Astrocyte Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# Display combined plot
print(combined_plot)

# Save plot
ggsave("astro_dotplot_log2fc_2025-01-15.png", combined_plot, width = 12, height = 8, dpi = 300)


#########################################
#########################################
#########################################
#########################################
#########################################
