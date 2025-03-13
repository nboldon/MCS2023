## 07_DESeq2.R
- Code does not run. 




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
setwd("/Volumes/DataBox/ProjMCS6/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

projMCS6

######################################################
######################################################
######################################################

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#summarizedexperiment-input

# SummarizedExperiment input
# If one has already created or obtained a SummarizedExperiment, it can be easily input into DESeq2 as follows. 
# First we load the package containing the airway dataset.

library("airway")
data("airway")
se <- airway

# The constructor function below shows the generation of a DESeqDataSet from a RangedSummarizedExperiment se.

library("DESeq2")
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE


######################################################
######################################################
######################################################

# ArchR Functions

frags <- getFragmentsFromProject(
  ArchRProj = projMCS6,
  cellNames = getCellNames(ArchRProj = projMCS6)[which(projMCS6@cellColData$Clusters == "C1")]
)

unlist(frags)

getCellColData(projMCS6)

getCellNames(projMCS6)

getChromLengths(projMCS6)
getChromSizes(projMCS6)
getExons(projMCS6)

######################################################
######################################################
######################################################

## From Claude.com

# Use the getGroupSE function to obtain a SummarizedExperiment object for your clusters:
groupSE <- getGroupSE(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters"
)

assayNames(groupSE)
colnames(colData(groupSE))
head(colData(groupSE))

# create a cluster vector from the row names:
clusters <- rownames(colData(groupSE))
# add this information to the colData:
colData(groupSE)$Clusters <- clusters

# Check the dimensions of your data
dim(assay(groupSE, "GeneScoreMatrix"))
# Check the number of unique clusters
length(unique(clusters))
# Check the number of samples per cluster
table(clusters)

# Create a DESeqDataSet
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(assay(groupSE, "GeneScoreMatrix")),
  colData = colData(groupSE),
  design = ~ Clusters
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results for all pairwise comparisons
clusters <- levels(dds$Clusters)
results_list <- list()

for (i in 1:(length(clusters) - 1)) {
  for (j in (i + 1):length(clusters)) {
    comparison <- paste(clusters[i], "vs", clusters[j])
    results_list[[comparison]] <- results(dds, contrast = c("Clusters", clusters[i], clusters[j]))
  }
}

# Analyze and visualize the results:
# Example: Plot top differentially expressed genes for a comparison
library(EnhancedVolcano)

comparison <- "Cluster1 vs Cluster2"
EnhancedVolcano(results_list[[comparison]],
                lab = rownames(results_list[[comparison]]),
                x = "log2FoldChange",
                y = "padj",
                title = comparison
)

# You can also extract significant genes based on adjusted p-value and fold change
sigGenes <- subset(results_list[[comparison]], padj < 0.05 & abs(log2FoldChange) > 1)
