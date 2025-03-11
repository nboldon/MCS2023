#Setup an interactive session
salloc --account=eon -t 0-04:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

conda activate seurat_2024-07

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

setwd("/project/eon/nboldon/MCS2023/Subset/C18_Subset")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

##############

#Load project
C18ArchRSubset <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Subset/Save-C18ArchRSubset",
                                   force = FALSE, showLogo = FALSE)

getAvailableMatrices(C18ArchRSubset)

table(C18ArchRSubset$Clusters)

######################################
######################################

C18ArchRSubset <- addHarmony(
  ArchRProj = C18ArchRSubset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

############################################

# Subset by Cluster
C18.1ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C1"]
C18.2ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C2"]
C18.3ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C3"]
C18.4ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C4"]
C18.5ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C5"]
C18.6ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C6"]
C18.7ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C7"]
C18.8ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C8"]
C18.9ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C9"]
C18.10ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C10"]
C18.11ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C11"]
C18.12ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C12"]
C18.13ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C13"]
C18.14ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C14"]
C18.15ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C15"]
C18.16ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C16"]
C18.17ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C17"]
C18.18ArchRSubset <- C18ArchRSubset[C18ArchRSubset$Clusters %in% "C18"]

############################################
############################################
############################################

## Cluster 1 

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C18.4ArchRSubset$Sample

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

# Check to make sure it worked
unique(treatment)

# Assign the treatment to the actual project
C18.4ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(C18.4ArchRSubset$treatment)

############
