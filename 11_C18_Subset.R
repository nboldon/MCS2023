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

setwd("/project/eon/nboldon/MCS2023/Subset/C18_Subset")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

##############

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)
C18ArchRSubset <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Subset/Save-C18ArchRSubset", force = FALSE, showLogo = FALSE)

getAvailableMatrices(C18ArchRSubset)

table(projMCS5$Clusters)

######################################
######################################

## Subset cluster/samples of interest

C18ArchRSubset <- projMCS5[
	projMCS5$Clusters %in% c("C18"),]
C18ArchRSubset

# numberOfCells(1): 44381
# medianTSS(1): 13.187
# medianFrags(1): 8504

subsetArchRProject(
  ArchRProj = projMCS5,
  cells = getCellNames(C18ArchRSubset),
  outputDirectory = "C18Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

#################################

getAvailableMatrices(C18ArchRSubset)
# "GeneScoreMatrix" "PeakMatrix"      "TileMatrix"  

#################################

## Recluster

#Dimensionality reduction with ArchR

C18ArchRSubset <- addIterativeLSI(
        ArchRProj = C18ArchRSubset,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        iterations = 4,
        clusterParams = list (
                resolution = c(0.1, 0.2, 0.4),
                sampleCells = 10000,
                n.start = 10
        ),
	varFeatures = 20000,
        dimsToUse = 1:30,
        seed = 123
)
# These are the default resolution parameters from ArchR

# When the iterative LSI approach isn't enough of a correction for strong batch effect differences, the batch effect correction tool called Harmony can be used. 
# This creates a new reducedDims object called "Harmony" in the MCS2 object

C18ArchRSubset <- addHarmony(
        ArchRProj = C18ArchRSubset,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
	force = TRUE
)
# Harmony converged after 3 iterations

########################################

#Clustering using Seurat's FindClusters() function

C18ArchRSubset <- addClusters(
        input = C18ArchRSubset,
        reducedDims = "Harmony",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        seed = 123,
	force = TRUE
)

# Number of nodes: 44381
# Number of edges: 1524932
# Maximum modularity in 10 random starts: 0.8597
# Number of communities: 18
# 2024-06-19 14:30:29.542865 : Testing Biased Clusters, 1.117 mins elapsed.
# 2024-06-19 14:30:29.71117 : Testing Outlier Clusters, 1.12 mins elapsed.
# 2024-06-19 14:30:29.716876 : Assigning Cluster Names to 18 Clusters, 1.12 mins elapsed.

#To access these clusters, use the $ accessor which shows the cluster ID for each single cell
head(C18ArchRSubset$Clusters)

#To tabulate the number of cells present in each cluster
table(C18ArchRSubset$Clusters)

#To better understand which samples reside in which clusters, create a cluster confusion matrix across each sample
cM <- confusionMatrix(paste0(C18ArchRSubset$Clusters), paste0(C18ArchRSubset$Sample))
cM

#To save the matrix as a .csv file
write.csv(cM, "C18ArchRSubset_Clusters_2024-06-19.csv")

#To plot the confusion matrix as a heatmap
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = "black"
)
p

#################################################
#################################################

##Single Cell Embeddings - UMAP

#To run a Uniform Manifold Approximation and Projection (UMAP)
C18ArchRSubset <- addUMAP(
        ArchRProj = C18ArchRSubset,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine",
	force = TRUE,
        seed = 123
)

#To plot the UMAP results by sample
p1 <- plotEmbedding(ArchRProj = C18ArchRSubset, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#To plot the UMAP results by cluster
p2 <- plotEmbedding(ArchRProj = C18ArchRSubset, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#To visualize these 2 plots side by side
ggAlignPlots(p1, p2, type = "h")

#To save an editable vectorized version of the plots:
plotPDF(p1,p2, name = "UMAP-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf", ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#################################################

##Single cell embeddings - tSNE

#To run a t-Stocastic Neighbor Embedding (tSNE)
C18ArchRSubset <- addTSNE(
        ArchRProj = C18ArchRSubset,
        reducedDims = "IterativeLSI",
        name = "TSNE",
        perplexity = 30,
	force = TRUE,
        seed = 123
)

#To plot the tSNE (the same parameters apply to colorBy and name regardless of the type of embedding used)
p1 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "TSNE"
)

p2 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "TSNE"
)

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "TSNE-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf", ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#######################################
#######################################

##Dimensionality reduction after Harmony for UMAP

#Used to assess the effects of Harmony by visualizing the embedding using UMAP or tSNE and comparing this to the embeddings visualized previously for iterative LSI.

#Repeat the UMAP embedding with the same parameters but for the "Harmony" reduced Dims object
C18ArchRSubset <- addUMAP(
        ArchRProj = C18ArchRSubset,
        reducedDims = "Harmony",
        name = "UMAPHarmony",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine",
	force = TRUE,
        seed = 123
)

p3 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "UMAPHarmony"
)

p4 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "UMAPHarmony"
)

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "UMAP2Harmony-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf", ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

##Dimensionality reduction after Harmony for tSNE

#Follow similar steps to those used for UMAP

C18ArchRSubset <- addTSNE(
        ArchRProj = C18ArchRSubset,
        reducedDims = "Harmony",
        name = "TSNEHarmony",
        perplexity = 30,
	force = TRUE,
        seed = 123
)

p7 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "TSNEHarmony"
)

p8 <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "TSNEHarmony"
)

ggAlignPlots(p7, p8, type = "h")

plotPDF(p1,p2,p7,p8, name = "TSNE2Harmony-Sample-Clusters_C18ArchRSubset_2024-06-19.pdf", ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#########################################
#########################################

##Gene scores and marker genes

#markerGS <- getMarkerFeatures(
#        ArchRProj = projMCSv2,
#        useMatrix = "GeneScoreMatrix",
#        groupBy = "Clusters",
#        bias = c("TSSEnrichment", "log10(nFrags)"),
#        testMethod = "wilcoxon"
#)

#markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

#Marker list by cluster

#markerList$C6

Cluster.markers <- getMarkerFeatures(
        ArchRProj = C18ArchRSubset,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clusters",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
#        maxCells = 4000
)

#markerList.Cluster.markers <- data.frame(getMarkers(Cluster.markers, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"))

#cluster.markers.test <- markerList.Cluster.markers

#To get a list of dataframe objects, one for each cluster, containing the relevant marker features
clusterMarkerList <- getMarkers(Cluster.markers, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

#Marker list by cluster
clusterMarkerList$C6

for(i in names(clusterMarkerList)) {
        write.csv(clusterMarkerList[[i]], file = paste(i, "_C18ArchRSubset_2024-06-19.csv", sep = ""))
}

#########################################

##Heatmaps for marker features

#To visualize all of the marker features simultaneously
markerGenes <- c(
        "Olig2", "Neurod6", "Gad2", "Dio2", "C1qa", "Cd31"
)

heatmapGS <- plotMarkerHeatmap(
        seMarker = Cluster.markers,
        cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
        labelMarkers = markerGenes,
        transpose = TRUE
)

#To plot the heatmap
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#To save the heatmap
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap_C18ArchRSubset_2024-06-19.pdf", width = 8, height = 6, ArchRProj = C18ArchRSubset, addDOC = FALSE)

#############################################

##Visualizing marker genes on an embedding
markerGenes <- c(
        "Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "GeneScoreMatrix",
        name = markerGenes,
        embedding = "UMAP",
        quantCut = c(0.01, 0.95),
        imputeWeights = NULL
)

#To plot a specific gene
p$Olig2

#To plot all genes, use cowplot to arrange the various marker genes into a single plot
p2 <- lapply(p, function(x){
        x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0,0,0,0), "cm")) +
        theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
        )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#To save an editable vectorized version of the plot
plotPDF(plotList = p,
        name = "UMAP-MarkerGenes-WO-Imputation_C18ArchRSubset_2024-06-19.pdf",
        ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#######################################################

##Use MAGIC to impute gene scores by smoothing signal across nearby cells

#Impute weights to the ArchRProject
C18ArchRSubset <- addImputeWeights(C18ArchRSubset)

#The impute weights can then be passed to plotEmbedding() when plotting gene scores overlayed on the UMAP embedding
markerGenes <- c(
        "Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "GeneScoreMatrix",
        name = markerGenes,
        embedding = "UMAP",
        imputeWeights = getImputeWeights(C18ArchRSubset)
)

#To subset the plot list to select a specific gene
p$Olig2

#To plot all the marker genes at once using cowplot
#Rearrange for grid plotting
p2<- lapply(p, function(x){
        x + guides(color = FALSE, fill = FALSE) +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0,0,0,0), "cm")) +
        theme(
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
        )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

#To save an editable vectorized version of the plot
plotPDF(plotList = p,
        name = "UMAP-MarkerGenes-W-Imputation_C18ArchRSubset_2024-06-19.pdf",
        ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#######################################################
#######################################################

#Save project
saveArchRProject(ArchRProj = C18ArchRSubset, outputDirectory = "/project/eon/nboldon/MCS2023/Subset/Save-C18ArchRSubset", load = FALSE)

#######################################################
#######################################################

# Markers from 5-21-2024 and  5-23-2024 gene lists

markerGenesAstro <- c(
        "Agt", "Aldh1l1", "Aldoc", "Apoe", "Aqp4", "Bmpr1b", "Cbs", "C4B", "CCK", "Clu", "Cryab",
        "Dio2", "Egfr", "Fam107a", "Fosb", "Gfap", "Gja1", "Gjb6", "Glis3", "Glul", "Gpc5", "Hes5", "Id4", "Igfbp2",
        "Irx1", "Itih3", "Mfge8", "Mlc1", "Nr2f2", "Nrp1", "Nrp2", "Nts", "Plcd4", "Ppp1r3c",
        "Prdm16", "S1pr1", "Serpinf1", "Slc1a3", "Slco1c1", "Sox9", "Sparcl1", "Stox1","Tbx5"
)
# Mett17a does not exist

markerGenesGlut <- c(
        "Adamts2", "ADCYAP1", "Agmat", "Arf5", "Arhgap25", "BARX2", "Batf3",
        "BCL11B", "Brinp3", "Calb1", "Car3", "Car12", "CDCP1", "Cdh13", "Chrna6", "Col18a1", "Col27a1",
        "Cpa6", "Cpne7", "Crh", "Ctgf", "Ctxn3", "Cux2", "Dmrtb1", "Drd1", "EMX1", "Endou", "Etv1",
        "Fam84b", "Fezf2", "FGF6", "Foxp2", "Fst", "Gm11549", "Gpr88", "Gpr139", "Grp", "GSX2",
        "Hpgd", "Hs3st4", "Hsd11b1", "Hsd17b2", "Kcng1", "Krt80", "Lamp5", "Lgr5", "Lhx5", "Lrg1",
        "Lrrk1", "Met", "Mgp", "Mup5", "Nefl", "Neurod1", "Neurod2", "Neurod6", "Neurog2", "Ngb",
        "Nov", "Npsr1", "Npw", "NPY2R", "Nr4a2", "Ntng2", "Nxph1", "Nxph3", "Olfm3", "Oprk1", "Ovol2",
        "P2ry12", "Pcdh19", "Pcp4", "Pde1c", "Plch1", "Pld5", "PROX1", "Prss12", "Ptgfr", "Ptgs2",
        "Rapgef3", "Rell1", "Rgs12", "Rorb", "Rprm", "Rrad", "Rspo1", "Rxfp1", "Satb2", "Scnn1a",
        "Serpinb11", "Sla", "SLC17A6", "Slc17a7", "Slc17a8", "Slco2a1", "SOX2", "SOX5", "Spon1", "Sstr2",
        "Sulf1", "Synpr", "Syt6", "Syt17", "Tbr1", "Tcerg1l", "Tgfb1", "TGM6", "Tle4", "Tnc", "Tox",
        "Tox2", "Tph2", "Trh", "Tshz2", "Tubb3", "Ucma", "VWA2", "Wfs1", "Wls"
)
# 3110035E14Rik,CTGL,Gaml3,ITL6GL,NKX2.1,NPGL,Pdzen4,Ptrf,Tbr2 do not exist

markerGenesGABA <- c(
        "4930553C11Rik", "ADAMTS17", "Adamts19", "ADARB2", "Ankk1", "Aqp5", "Arx", "B3gat2", "Btg2", "C1ql1",
        "Calb1", "Calb2", "Car4", "Car10", "Cbln4", "Cdk6", "Chat", "CHODL", "Chst9", "Cnr1", "Cpne5",
        "Crh", "CRHBP", "Cxcl14", "CYP26A1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Drd1", "Drd2", "Efemp1", "Eomes",
        "Esm1", "Etv1", "Eya1", "Fam114a1", "Fibin", "FLT3", "FOXP2", "Fut9", "Gabrd", "Gabrg1", "GAD1",
        "GAD2", "Glra3", "Gpc3", "Gpr6", "Gpr50", "Gpx3", "Hes5", "Hmcn1", "Htr1f", "Htr2c", "Htr3a",
        "Id2", "Islr", "Itga4", "Itih5", "IZUMO1R", "Kank4", "Kcne4", "Krt73", "Lamp5", "Lhx5", "Lhx6",
        "Lsp1", "Manea", "MC4R", "MEIS2", "Mybpc1", "Myh4", "Myh8", "Myl1", "Ndnf", "Necab1", "Neurog2", "Nkx2-1",
        "Nos1", "Nptx2", "Npy", "Npy2r", "NR2E1", "NR2F2", "Nts", "NXPH1", "Obox3", "Oxtr", "Parm1",
        "PAX6", "Pde1a", "Pdlim5", "Penk", "Pkp2", "Plekhh1", "Pltp", "Pnoc", "Prdm8", "Ptgdr", "Pvalb", "Reln",
        "Rspo2", "Satb1", "Sema3c", "SIX3", "SLC10A4", "Slc17a8", "Slc18a3", "SLC32A1", "Slc6a1", "Smad3",
        "Sncg", "Sox2", "Sp8", "SST", "Sv2c", "Synpr", "Tac1", "Tacr3", "Tacstd2", "Th",
        "Tmem182", "Tpbg", "Unc5B", "VAX1", "Vip", "Vipr2", "Wt1", "ZAR1"
)
# Clm1,Fam159b,Hltr1d,Lgtp,Nxph1/2,SFTA3,Tac3,Tactd2 do not exist

markerGenesOligo <- c(
        "9630013A20Rik", "Aspa", "Ccnb1", "Cldn11", "C1ql1", "Cnksr3", "DLL1", "EGFR", "Gjc3", "Gpr17", "Grm5",
        "Hapln2", "Idh1", "Itpr2", "Kif5a", "Lims2", "Mag", "Mal", "Mbp", "Mobp", "Mog", "MKI67", "Neu4",
        "Nkx2-2", "Olig1", "Olig2", "Olig3","OPALIN", "PCDH15", "Pdgfra", "Plekhb1",
        "Plp1", "Ppp1r14a", "Prom1", "Qdpr", "Rassf10", "RBFOX1", "Rnf122", "Serpinb1a", "SOX8", "SOX10",
        "ST18", "Stmn1","Synpr", "Tcf7l2", "Ugt8a"
)
# Olig4,Olig5,Olig6 do not exist

markerGenesEndo <- c(
        "Acta2", "Anpep","Bmx", "Bsg","CCNE2", "Cdh5", "Cldn5", "Ctla2a", "Cytl1",
        "Emcn", "ERG", "Flt1", "FN1", "Ifit3", "Isg15", "Itm2a", "Klf2", "Lgals1", "Ly6c1", "Mcam",
        "Pltp", "Ptprb", "RHOH", "SELL", "Slco1a4", "Tek", "Tyrobp", "VIM", "Vtn", "Xdh"
)
# Cd31,Co11a1 do not exist 

markerGenesVasc <- c(
        "ACTA2", "Myl9", "RUNX3", "Cd74", "Col15a1", "Hs3st6", "Mc5r"
)
# GFI2 does not exist 

markerGenesPericyte <- c(
        "Aif1", "Cald1", "Cspg4", "Cx3cr1", "Flt1", "Gng11", "GRM8", "Higd1b", "Igfbp7",
        "Kcnj8", "Lama2", "Lyl1", "Lyve1", "Mrc1", "Myl9", "Ndufa4l2", "Pdgfrb", "Rgs5", "Spic"
)

markerGenesMicroglia <- c(
        "Aif1", "B2m", "C1qa", "C1qb", "C1QC", "Ccl4", "Cd14", "Cx3cr1", "Csf1r", "Ctss", "DOCK8",
        "EGFR", "FOXO1", "GATA4", "Hexb", "Lyz2", "MKI67", "Mrc1", "P2ry12", "Ptprc", "Siglech",
        "SMAD3", "STAT1", "Tmem119", "Trem2", "Tyrobp"
)
# FYB1 does not exist 

markerGenesRadialGlia <- c(
        "Dbi", "Fabp7", "HES1", "Hes5", "Hopx", "Mlc1", "Slc1a3", "Thrsp"
)
# GLAST,RC2 do not exist

###########################################

##Use MAGIC to impute gene scores by smoothing signal across nearby cells

#Impute weights to the ArchRProject
C18ArchRSubset <- addImputeWeights(C18ArchRSubset)

#The impute weights can then be passed to plotEmbedding() when plotting gene scores overlayed on the UMAP embedding
p <- plotEmbedding(
        ArchRProj = C18ArchRSubset,
        colorBy = "GeneScoreMatrix",
        name = markerGenesRadialGlia,
        embedding = "UMAPHarmony",
        imputeWeights = getImputeWeights(C18ArchRSubset)
)

#To save an editable vectorized version of the plot
plotPDF(plotList = p,
        name = "UMAP-RadialGliaMarkerGenes-W-Imputation_C18ArchRSubset_2024-06-21.pdf",
        ArchRProj = C18ArchRSubset, addDOC = FALSE, width = 5, height = 5)

#######################################################
#######################################################

## Highlight cells of interest

################

# Subset samples of interest, add embedding and impute weights

t3ArchRSubset <-C18ArchRSubset[
        C18ArchRSubset$Sample %in% c("C305_", "C307_", "C313_", "C316_", "C320_", "C322_", "C325_", "C334_", "C340",
        "C341_", "C345_", "C350_", "C352_", "C359_", "C364_"),]

t3ArchRSubset

#################

p2<- plotEmbedding(
  ArchRProj = C18ArchRSubset,
  embedding = "UMAPHarmony",
  colorBy = "cellColData",
  name = "Sample",
  size = 1,
  sampleCells = NULL,
  highlightCells = getCellNames(ArchRProj = C18ArchRSubset)[which(C18ArchRSubset@cellColData$Sample %in% c("C305_", 
	"C307_", "C313_", "C316_", "C320_", "C322_", "C325_", "C334_", "C340",
        "C341_", "C345_", "C350_", "C352_", "C359_", "C364_"))],
  baseSize = 10,
  plotAs = "points")

#To save an editable vectorized version of the plot

plotPDF(p2, name = "UMAP_T3_C18ArchRSubset_2024-06-21.pdf", width = 5, height = 5, ArchRProj = C18ArchRSubset, addDOC = FALSE)
