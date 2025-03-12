#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=64G --nodes=1 --ntasks-per-node=16

#Load required dependencies
module load miniconda3/4.12.0
conda activate /pfs/tc1/project/eon/archr_env

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/Analysis")
addArchRGenome("mm10")
addArchRThreads (threads = 16)

#Load project
projMCS2 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Analysis/Save-ProjMCS2", force = FALSE, showLogo = FALSE)

######################################
######################################

##Track plotting with ArchRBrowser

markerGenes  <- c(
    "Aqp4", "Aldh1l1", "Mlc1", "Cbs", "Ppp1r3c", "Plcd4", "Dio2", #Astrocyte
    "Cldn11", "Cx3cr1", "Csf1r", "Sparc", "Trem2", "Ccl4", "Cd14", "Tyrobp", "C1qa", #Microglia
    "Olig1", "Mbp", "Opalin", "Mag", "Mog", "Cldn11", "Ugt8a", "Olig2", #Oligodendrocyte
    "Spock3", "Gad1", "Grin3a", "Adarb2", "Grik1", "Lhx6", "Pvalb", "Gad2",  #GABAergic
    "Sulf1", "Slc17a8", "Tshz2", "Slc17a6", "Neurod6", #Glutamatergic
        "Cldn5" #"CD31"DoesNotExist #Endothelial
  )

p <- plotBrowserTrack(
	ArchRProj = projMCS5,
	groupBy = "Clusters",
	geneSymbol = markerGenes,
	upstream = 50000,
	downstream = 50000
)

#To track a specific gene
grid::grid.newpage()
grid::grid.draw(p$Olig2)

#To save a multi-page PDF with a single page for each gene locus in the plot
plotPDF(plotList = p,
	name = "Plot-Tracks-MarkerGenes_2024-04-04.pdf",
	ArchRProj = projMCS5,
	addDOC = FALSE, width = 5, height = 5)
