#Setup an interactive session
salloc --account=eon -t 0-06:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtree)
library(aplot)
library(patchwork)
library(ggplot2)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ClusterID")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

############################################
############################################

# Assigning clusters with gene scores

proj <- addUMAP(ArchRProj = projMCS5, reducedDims = "IterativeLSI2", force = TRUE)

#proj <- addImputeWeights(proj)

# Overlay marker gene scores in 2D UMAP embedding

markerGenes  <- c(
    "Aqp4", "Aldh1l1", "Mlc1", "Cbs", "Ppp1r3c", "Plcd4", "Dio2", #Astrocyte
    "Cldn11", "Cx3cr1", "Csf1r", "Sparc", "Trem2", "Ccl4", "Cd14", "Tyrobp", "C1qa", #Microglia
    "Olig1", "Mbp", "Opalin", "Mag", "Mog", "Cldn11", "Ugt8a", "Olig2", #Oligodendrocyte
    "Spock3", "Gad1", "Grin3a", "Adarb2", "Grik1", "Lhx6", "Pvalb", "Gad2",  #GABAergic
    "Sulf1", "Slc17a8", "Tshz2", "Slc17a6", "Neurod6", #Glutamatergic
	"Cldn5" #"CD31"DoesNotExist #Endothelial
  )

# 04/09/2024 Additions
markerGenes <- c(
        "Sst","Frmd7", "Sp8", "Vipr2", "Etv1", "Ntsrl", "Lhx8", "Zic4", "Myo3a", "Lhx3", "Vip", "Drd2", "Drd1", #GABA$
        "Cux1", "Cux2", "Satb2", "Bcl11b", "Foxp2", "Pcp4", "Neurod1", "Neurod2", "Tbr1", "Synj1", #Glutamatergic
        "NF1", "Gfap", #Astrocyte
        "Pdgfra", #Oligodendrocyte
        "Anpep", "COUP-TFII" #Endothelial
)

p <- plotEmbedding(
    ArchRProj = projMCS5, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

# To plot a specific gene, subset the plot list using the gene name

p$Olig1

# To plot all genes, use cowplot to arrange the different plots together.
# Each of these marker genes lights up the corresponding cell clusters. 

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# Save an editable vectorized version of the plots

plotPDF(plotList = p, 
    name = "UMAP-Marker-Genes-W-Imputation_2024-04-04.pdf", 
    ArchRProj = projMCS5, 
    addDOC = FALSE, width = 15, height = 15)
