
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


## Mean Diff dotplot


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

#Add rownames
rownames(mean_acc) <- rowData(Cluster.markers)$name
rownames(rel_acc) <- rowData(Cluster.markers)$name
rownames(FDR_acc) <- rowData(Cluster.markers)$name
rownames(Log2FC_acc) <- rowData(Cluster.markers)$name

#Add column names
colnames(mean_acc) <- colnames(Cluster.markers)
colnames(rel_acc) <- colnames(Cluster.markers)
colnames(FDR_acc) <- colnames(Cluster.markers)
colnames(Log2FC_acc) <- colnames(Cluster.markers)

# Melt the data
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc))
melted_rel_acc <- reshape2::melt(as.matrix(rel_acc))
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc))
melted_Log2FC_acc <- reshape2::melt(as.matrix(Log2FC_acc))
colnames(melted_mean_acc) <- c("gene", "cluster", "mean")
colnames(melted_rel_acc) <- c("gene", "cluster", "rel")
colnames(melted_FDR_acc) <- c("gene", "cluster", "FDR")
colnames(melted_Log2FC_acc) <- c("gene", "cluster", "Log2FC")
melted_mean_acc$gene <- as.character(melted_mean_acc$gene)
melted_rel_acc$gene <- as.character(melted_rel_acc$gene)
melted_FDR_acc$gene <- as.character(melted_FDR_acc$gene)
melted_Log2FC_acc$gene <- as.character(melted_Log2FC_acc$gene)

# Combine all data
acc_comb <- melted_mean_acc %>%
  left_join(melted_rel_acc, by = c("gene", "cluster")) %>%
  left_join(melted_FDR_acc, by = c("gene", "cluster")) %>%
  left_join(melted_Log2FC_acc, by = c("gene", "cluster"))

## Mean Diff dotplot


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
glut_main_markerGenes <- c("Cux2","Neurod1", "Neurod2", "Neurod6", "Satb2", "Tbr1")
glut_allUMAP_markerGenes <- c("Adcyap1","Agmat","Arhgap25","Barx2","Calb1","Car12","Cdcp1",
                              "Cdh13","Chrna6","Cpa6","Cux2","Dmrtb1","Drd1","Emx1","Fezf2",
                              "Fgf6","Foxp2","Gm11549","Gpr88","Gpr139","Grp","Hs3st4","Kcng1",
                              "Krt80","Lamp5","Lhx5","Lrg1","Mup5","Nefl","Neurod1", "Neurod2",
                              "Neurod6","Ngb","Nov","Npsr1","Npw","Npy2r","Nr4a2","Ntng2","Nxph3",
                              "Olfm3","Ovol2","Pcdh19","Pcp4","Plch1","Pld5","Prss12","Ptgs2",
                              "Rorb","Rprm","Rspo1","Rxfp1","Satb2","Serpinb11","Slc17a6","Slc17a7",
                              "Slc17a8","Slco2a1","Sox5","Spon1","Sstr2","Sulf1","Synpr","Syt6","Syt17",
                              "Tbr1","Tcerg1l","Tgm6","Tle4","Tox","Trh","Tshz2","Tubb3","Ucma",
                              "Wfs1")

gaba_main_markerGenes <- c("Drd2", "Gad1", "Pvalb", "Sp8", "Sst", "Vip")
gaba_allUMAP_markerGenes <-c ("Adamts17", "Adamts19", "Adarb2", "Ankk1", "Aqp5", "Calb1", 
                              "Calb2", "Car4", "Car10", "Cbln4", "Chat", "Cnr1", "Cpne5",
                              "Crh", "Crhbp", "Cyp26a1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", 
                              "Drd1", "Drd2", "Flt3", "Gabrd", "Gad1", "Gad2", "Gpr6", 
                              "Gpx3", "Htr2c", "Htr3a", "Izumo1r", "Lamp5", "Lhx5", "Lhx6",
                              "Myl1", "Necab1", "Nkx2-1", "Nos", "Npy", "Npy2r", "Nr2e1",
                              "Nxph1", "Oxtr", "Parm1", "Penk", "Pnoc", "Pvalb", "Reln",
                              "Rspo2", "Satb1", "Sema3c", "Slc10a4", "Slc18a3", "Slc32a1",
                              "Slc6a1", "Sncg", "Sp8", "Sst", "Sv2c", "Synpr", "Tac1",
                              "Tacr3", "Th", "Tmem182", "Tpbg", "Vax1", "Vip", "Wt1")

microglia_main_markerGenes <- c("C1qa", "Ccl4", "Cd14", "Cx3cr1", "Csf1r", "Trem2")
microglia_allUMAP_markerGenes <- c("Aif1", "B2m", "C1qa", "C1qb", "C1qc", "Ccl4", "Cd14",
                                   "Cx3cr1", "Csf1r", "Ctss", "Dock8", "Hexb", "Lyz2", "Mrc1",
                                   "P2ry12", "Ptprc", "Siglech", "Smad3", "Tmem119", "Trem2", "Tyrobp")

oligo_main_markerGenes <- c("Mag", "Mbp", "Olig1", "Opalin", "Pdgfra", "Ugt8a")
oligo_allUMAP_markerGenes <- c("9630013A20Rik", "Aspa", "Cldn11", "C1ql1", "Cnksr3", "Dll1",
                               "Gjc3", "Gpr17", "Grm5", "Hapln2", "Lims2", "Mag", "Mal", "Mbp",
                               "Mobp", "Mog", "Neu4", "Nkx2-2", "Olig1", "Olig2", "Olig3", "Opalin",
                               "Pcdh15", "Pdgfra", "Plekhb1", "Plp1", "Ppp1r14a", "Qdpr", "Rnf122",
                               "Serpinb1a", "Sox8", "Sox10", "St18", "Stmn1", "Tcf7l2", "Ugt8a")

endo_main_markerGenes <- c("Anpep", "Cldn5", "Flt1", "Fn1", "Itm2a", "Slco1a4")
endo_allUMAP_markerGenes <- c("Anpep", "Bmx", "Bsg", "Cdh5", "Cldn5", "Ctla2a", "Emcn",
                              "Erg", "Flt1", "Fn1", "Isg15", "Itm2a", "Klf2", "Mcam",
                              "Pltp", "Ptprb", "Slco1a4", "Tek", "Vim", "Vtn", "Xdh")

astro_main_markerGenes <- c("Aqp4", "Cbs", "Dio2", "Gfap", "Mlc1", "Ppp1r3c")
astro_allUMAP_markerGenes <- c("Agt", "Aldh1l1", "Aldoc", "Apoe", "Aqp4", "Bmpr1b", "Cbs", 
                           "Clu", "Dio2", "Egfr", "Fam107a", "Gfap", "Gja1", "Gjb6", "Glis3",
                           "Glul", "Gpc5", "Hes5", "Id4", "Igfbp2", "Itih3", "Mfge8", "Mlc1",
                           "Plcd4", "Ppp1r3c", "Prdm16", "S1pr1", "Slc1a3", "Slco1c1", "Sox9",
                           "Sparcl1")

##############################################################


## Run the following for each set of genes of interest

# Filter data for genes of interest
acc_filtered <- acc_comb %>% 
  filter(gene %in% glut_allUMAP_markerGenes)

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
  plot_annotation(title = "Glutamatergic Gene Dotplot with Dendrogram",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# Display combined plot
print(combined_plot)

# Save plot
ggsave("glut_allUMAP_dotplot-dendrogram_2025-01-16.png", combined_plot, width = 12, height = 13, dpi = 300)



#########################################
#########################################
#########################################
#########################################
#########################################
