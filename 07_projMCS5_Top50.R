## 07_projMCS5_Top50.R
- No subset, groupBy cluster, cutOff = "FDR <= 0.01 & Log2FC >= 1.25"
    - File name ex: GeneMarker_List_Master.csv & C1_GeneMarker_List.csv
- Filters top genes and plots heatmap.
    - File: GeneMarkers_Top50_1-25-2024.pdf




#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

############################################
############################################

projMCS5

getAvailableMatrices(projMCS5)

############################################

##Get marker features
markerGS <- getMarkerFeatures(
  ArchRProj = projMCS5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 44500
)
## Increased maxCells from 4000

#Get te list of markers
markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

#Display marker genes & their info for cluster 6
markerList$C6

write.csv(markerList, file = "GeneMarker_List_Master.csv")
write.csv(markerList$C1, file = "C1_GeneMarker_List.csv")
write.csv(markerList$C2, file = "C2_GeneMarker_List.csv")
write.csv(markerList$C3, file = "C3_GeneMarker_List.csv")
write.csv(markerList$C4, file = "C4_GeneMarker_List.csv")
write.csv(markerList$C5, file = "C5_GeneMarker_List.csv")
write.csv(markerList$C6, file = "C6_GeneMarker_List.csv")
write.csv(markerList$C7, file = "C7_GeneMarker_List.csv")
write.csv(markerList$C8, file = "C8_GeneMarker_List.csv")
write.csv(markerList$C9, file = "C9_GeneMarker_List.csv")
write.csv(markerList$C10, file = "C10_GeneMarker_List.csv")
write.csv(markerList$C11, file = "C11_GeneMarker_List.csv")
write.csv(markerList$C12, file = "C12_GeneMarker_List.csv")
write.csv(markerList$C13, file = "C13_GeneMarker_List.csv")
write.csv(markerList$C14, file = "C14_GeneMarker_List.csv")
write.csv(markerList$C15, file = "C15_GeneMarker_List.csv")
write.csv(markerList$C16, file = "C16_GeneMarker_List.csv")
write.csv(markerList$C17, file = "C17_GeneMarker_List.csv")
write.csv(markerList$C18, file = "C18_GeneMarker_List.csv")
write.csv(markerList$C19, file = "C19_GeneMarker_List.csv")
write.csv(markerList$C20, file = "C20_GeneMarker_List.csv")
write.csv(markerList$C21, file = "C21_GeneMarker_List.csv")
write.csv(markerList$C22, file = "C22_GeneMarker_List.csv")
write.csv(markerList$C23, file = "C23_GeneMarker_List.csv")
write.csv(markerList$C24, file = "C24_GeneMarker_List.csv")
write.csv(markerList$C25, file = "C25_GeneMarker_List.csv")

#Get lengths of markerlist for each cluster
mlistlens <- lapply(markerList, nrow)

#Remove any clusters with zero markers
markerList <- markerList[mlistlens != 0]


filt_num <- 5000
top_genes <- lapply(markerList, function(x){ # for each element of markerList
  if(nrow(x) >=filt_num){ # if there are 5 entries
    x[1:filt_num,]$name   # get the top 5
  }else{ # if there are not 50 entries
    x[1:nrow(x),]$name # get all the entries
  }
})

#Get all the genes unlisted and only the unique ones
unique(unlist(top_genes))
length(unique(unlist(top_genes)))

str(top_genes)
write.table(top_genes, file = "Top_Genes_Master.csv")
#read.table(top_genes)
#Error in read.table(top_genes, header = TRUE, sep = "\t") : 
#  'file' must be a character string or connection

assays(markerGS)

#Take a look at how this is encoded
markerGS@assays@data$Mean[1:10,1:3]
markerGS@assays@data$MeanDiff[1:10,1:3]
markerGS@assays@data$FDR>0

#Get out mean and relative accessibility
mean_acc <- markerGS@assays@data$Mean
rel_acc  <- markerGS@assays@data$MeanDiff
FDR_acc <- markerGS@assays@data$FDR

#Add rownames
rownames(mean_acc) <- rowData(markerGS)$name
rownames(rel_acc) <- rowData(markerGS)$name
rownames(FDR_acc) <- rowData(markerGS)$name

#Add column names
colnames(mean_acc) <- colnames(markerGS)
colnames(rel_acc) <- colnames(markerGS)
colnames(FDR_acc) <- colnames(markerGS)

#Melt it
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc))
melted_rel_acc <- reshape2::melt(as.matrix(rel_acc))
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc))
colnames(melted_mean_acc) <- c("gene", "cluster", "mean")
colnames(melted_rel_acc) <- c("gene", "cluster", "rel")
colnames(melted_FDR_acc) <- c("gene", "cluster", "FDR")
melted_mean_acc$gene <- as.character(melted_mean_acc$gene)
melted_rel_acc$gene <- as.character(melted_rel_acc$gene)
melted_FDR_acc$gene <- as.character(melted_FDR_acc$gene)

#Merge em together
acc_comb1 <- merge(melted_mean_acc, melted_rel_acc)
acc_comb2 <- merge (acc_comb1, melted_FDR_acc)
tail(acc_comb2)

################################################
################################################

#https://davemcg.github.io/post/lets-plot-scrna-dotplots/

markers <- unique(unlist(top_genes))

markers

#Marker_ID <- c(
#  "Olig2", "Neurod6", "Gad2", "Dio2", "C1qa", "Anpep"
#)
#markers_1 <-c (
#  "Vip", "Sst", "Pvalb", "Pdgfra", "Olig2", "Olig1", "Neurod6",
#  "Neurod2", "Neurod1", "Mog", "Gfap", "Gad2", "Gad1", "Drd2",
#  "Drd1", "Dio2", "Cx3cr1", "C1qa", "Anpep"
#)


## This doesn't subset anything
# pdf(file="test_dotplot_smh.pdf", width=10, height=10)  
# ggplot(acc_comb2, aes(x=cluster, y = gene, color = rel, size = mean)) + 
#   geom_point() 
# dev.off()

#Remove dots where there is zero (or near zero) expression
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point()

#Better color, better theme, rotate x axis labels
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  scale_color_viridis_c(name = 'log2 (count + 1)') +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

#Tweak color scaling
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

acc_comb2

#########################################

# make data square to calculate euclidean distance
mat <- acc_comb2 %>%
  filter(gene %in% markers) %>%
  # filter(mean > 0, abs(rel) > 1.0) %>%
  select(-mean, -FDR) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = rel) %>%
  data.frame() # make df as tibbles -> matrix annoying

row.names(mat) <- mat$gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
mat <- as.matrix(mat)
clust <- hclust(dist(mat)) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

dotplot <- acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

cowplot::plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

dotplot <- acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  #gene = factor(gene, levels = clust$labels[clust$order]) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

cowplot::plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

#########################################
#########################################
#https://davemcg.github.io/post/lets-plot-scrna-dotplots/

#To group genes (cluster) by similar expression patterns and show the dendrogram

#hclust 
#mutate(gene = factor(markers, levels = YOURNEWGENEORDER))


#################################################
#################################################

pdf(file="GeneMarker_Top50_Markermap.pdf", width=15, height=6)
plotMarkerHeatmap(   ### Copied this over from the full function reference - leaving anything not commented at defaul$
  seMarker = markerGS, ## Set this to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  #Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)
dev.off()

##############################################

hmap_to_subset <- plotMarkerHeatmap(  #Anything not commented left as default
  seMarker = markerGS,   ##Set this to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  ##Label the markers of interest to plot
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

## Get the matrix of the heatmap to figure out row indices for the genes of interest

hmap_to_subset_MAT <- plotMarkerHeatmap(   #Anything not commented left as default
  seMarker = markerGS,   ##Set to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  ##Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  #Return the matrix rather than a heatmap
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)


hmap_row_idx <- which(rownames(hmap_to_subset_MAT) %in% unique(unlist(top_genes))) #Get the row numbers of the heatmap matrix that match the known markers

#Plot it out
pdf(file="GeneMarkers_Top50_1-25-2024.pdf", width=7, height=5)
hmap_to_subset[hmap_row_idx, ]
dev.off()

