
#Setup an interactive session
salloc --account=eon -t 0-10:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
library(tibble)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/Tx_Comp")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

projMCS5

getAvailableMatrices(projMCS5)

############################################

mat<-read.csv("./C1_byTx_zscores_2024-03-21.csv",header=TRUE,row.names=1,
              check.names=FALSE)

head(mat)

################

# Makes a heatmap
#pdf(file="C1_byTx_pheatmap_2024-03-25.pdf", width=15, height=30)
#pheatmap(mat,scale="row",
#         color=colorRampPalette(c("navy", "white", "red"))(50))
#dev.off()

###############

# Create data frame for annotations

# Assuming mat is your matrix and you've correctly set its column names
# Create dfh with treatment and dex columns
dfh <- data.frame(treatment = as.character(colnames(mat)), dex = "Genotype")

# Convert 'treatment' column to row names
dfh <- tibble::column_to_rownames(dfh, var = "treatment")

# Update 'dex' column based on condition
dfh$dex <- ifelse(rownames(dfh) %in% c("t1", "t2"), "2N", "Ts")

# View the updated data frame
dfh


# Plot the pheatmap

pdf(file="C1_byTx_sorted-pheatmap_2024-03-25.pdf", width=7, height=25)
pheatmap(mat,scale="row", annotation_col = dfh,
         annotation_colors=list(dex=c("2N"="orange","Ts"="black")),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2, cutree_rows=2,
         main="Accessibility and clustering of C1 top DE genes",
         fontsize=11, cellwidth=35, cellheight=10.25)
dev.off()

