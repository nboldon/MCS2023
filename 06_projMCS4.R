Add pseudobulk replicates & reproduceable peak set using MACS2


#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=64G --nodes=1 --ntasks-per-node=16

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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS2 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS2", force = FALSE, showLogo = FALSE)

############################################
##############################################

## Pseudo-bulk replicates in ArchR

projMCS4 <- addGroupCoverages(ArchRProj = projMCS2, groupBy = "Clusters") 

# 10 Calling Peaks w/ Macs2

pathToMacs2 <- findMacs2()
projMCS4 <- addReproduciblePeakSet(
    ArchRProj = projMCS4, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

getPeakSet(projMCS4)

saveArchRProject(ArchRProj = projMCS4, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS4", load = FALSE)

