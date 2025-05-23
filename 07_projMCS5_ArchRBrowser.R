## 07_projMCS5_ArchRBrowser.R
- Additional browser track gene regions of interest
- Subset by treatment group
- C18 Subcluster subset by treatment group
- Files generated using this code:
  - 07_Browser-Tracks_C18-grpByTx_2024-06-05.pdf
  - 07_Browser-Tracks_T3-grpByClusters_2024-06-05.pdf





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
library(rhandsontable)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ClusterID")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

######################################
######################################

##Track plotting with ArchRBrowser

markerGenes <- c(
        "Gad1", "Gad1", "Girk1", "Spock3", "Cux2", "Lhx6", "Pvalb", "Grin3a", "Adarb2", 
	"Slc17A6", "Slc17a8", "Sulf1", "Rorb", "Tle4", "Tshz2", "S100b", "Gfap", "Mag",
	"Mog", "Mbp", "Oligo1", "Opalin", "C1qa", "CLDN5", "CX3CR1", "AQP4", "Sparc",
	"Ppp1r3c", "Aldh1l1", "Plcd4", "Mlc1", "Cbs", "Trem2", "Ccl4", "Csf1r", "Cd14",
	"Tyrobp", "Cldn11", "Ugt8a", "Reln", "Tek", "Cdh5"
)

# 04/09/2024 Additions
markerGenes <- c(
	"Sst","Frmd7", "Sp8", "Vipr2", "Etv1", "Ntsrl", "Lhx8", "Zic4", "Myo3a", "Lhx3", "Vip", "Drd2", "Drd1", #GABAergic
	"Cux1", "Cux2", "Satb2", "Bcl11b", "Foxp2", "Pcp4", "Neurod1", "Neurod2", "Tbr1", "Synj1", #Glutamatergic
	"NF1", "Gfap", #Astrocyte
	"Pdgfra", #Oligodendrocyte
	"Anpep", "COUP-TFII" #Endothelial
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
        name = "Plot-Tracks-GeneMarkers_2024-03-19.pdf",
        ArchRProj = projMCS5,
        addDOC = FALSE, width = 5, height = 5)

######################################
######################################

# Subsetting by Tx group for comparison - Option 1

t1ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C302_", "C306_", "C309_", "C318_", "C323_", "C328_", "C332_", "C337_", "C339_",
        "C346_", "C351_", "C353_", "C360_"),]

t1ArchRSubset
#NumberOfCells: 24969
#medianTSS: 14.968
#medianFrags: 7696

#####

t2ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C304_", "C308_", "C312_", "C315_", "C321_", "C324_", "C327_", "C330_", "C333_",
        "C336_", "C342_", "C348_", "C349_", "C355_", "C358_", "C362_"),]

t2ArchRSubset
#NumberOfCells: 30149
#medianTSS: 14.749
#medianFrags: 6860

#####

t3ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C305_", "C307_", "C313_", "C316_", "C320_", "C322_", "C325_", "C334_", "C340",
        "C341_", "C345_", "C350_", "C352_", "C359_", "C364_"),]

t3ArchRSubset
#NumberOfCells: 25923
#medianTSS: 14.29
#medianFrags: 6507

#####

t4ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C301_", "C303_", "C310_", "C314_", "C319_", "C335_", "C338_", "C344_", "C354_",
        "C356_", "C361_", "C363"),]

t4ArchRSubset
#NumberOfCells: 19728      
#medianTSS: 14.0195      
#medianFrags: 5697.5

###############

##Track plotting with ArchRBrowser

markerGenes <- c(
        "Tmem176a", "Mir6238", "Sp5", "Ces1f", "Smoc2", "Dio2", "Svs6")

p <- plotBrowserTrack(
        ArchRProj = t3ArchRSubset,
        groupBy = "Clusters",
        geneSymbol = markerGenes,
        upstream = 2000,
        downstream = 2000
)

# 6/4/24 up/downstream set to 50000
# 6/5/24 up/downstream set to 2000

#To track a specific gene
grid::grid.newpage()
grid::grid.draw(p$Olig2)

#To save a multi-page PDF with a single page for each gene locus in the plot
plotPDF(plotList = p,
        name = "Browser-Tracks_T3-grpByClusters_2024-06-05.pdf",
        ArchRProj = t3ArchRSubset,
        addDOC = FALSE, width = 5, height = 5)

######################################
######################################

# Subsetting by Tx group for comparison - Option 2
## This option is recommended by Sean when comparing treatment groups

C18ArchRSubset <- projMCS5[projMCS5$Clusters %in% "C18"]

# Specify which treatment group each sample is in:

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C18ArchRSubset$Sample

treatment <- gsub("C302_", "t1", treatment)
treatment <- gsub("C306_", "t1", treatment)
treatment <- gsub("C309_", "t1", treatment)
treatment <- gsub("C318_", "t1", treatment)
treatment <- gsub("C323_", "t1", treatment)
#C323 is T1 (not T3)
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
C18ArchRSubset$treatment <- treatment

# Check that this worked - if not, make sure the previous line was run successfully
head(C18ArchRSubset$treatment)

########################

##Track plotting with ArchRBrowser

markerGenes <- c(
        "Tmem176a", "Mir6238", "Sp5", "Ces1f", "Smoc2", "Dio2", "Svs6")

p <- plotBrowserTrack(
        ArchRProj = C18ArchRSubset,
        groupBy = "treatment",
        geneSymbol = markerGenes,
        upstream = 2000,
        downstream = 2000
)

# 6/4/24 up/downstream set to 50000
# 6/5/24 up/downstream set to 2000

#To track a specific gene
grid::grid.newpage()
grid::grid.draw(p$Olig2)

#To save a multi-page PDF with a single page for each gene locus in the plot
plotPDF(plotList = p,
        name = "Browser-Tracks_C18-grpByTx_2024-06-05.pdf",
        ArchRProj = C18ArchRSubset,
        addDOC = FALSE, width = 5, height = 5)

######################################
######################################

