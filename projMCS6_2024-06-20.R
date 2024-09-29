#Setup an interactive session
salloc --account=eon -t 09:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

###########################################
###########################################

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

############################################
############################################

##Dimensionality reduction with ArchR

projMCS6 <- addIterativeLSI(
        ArchRProj = projMCS5,
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
# This creates a new reducedDims object called "Harmony" in the MCS5 object

projMCS6 <- addHarmony(
        ArchRProj = projMCS6,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
	force = TRUE
)
# Harmony converged after 3 iterations

#############################################

#To access these clusters, use the $ accessor which shows the cluster ID for each single cell
head(projMCS6$Clusters)

#To tabulate the number of cells present in each cluster
table(projMCS6$Clusters)

#To better understand which samples reside in which clusters, create a cluster confusion matrix across each sample
cM <- confusionMatrix(paste0(projMCS6$Clusters), paste0(projMCS6$Sample))
cM

#To save the matrix as a .csv file
write.csv(cM, "ProjMCS6_IterativeLSI-wHarmony_2024-06-20.csv")

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
projMCS6 <- addUMAP(
        ArchRProj = projMCS6,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine",
	force = TRUE,
        seed = 123
)

#To plot the UMAP results by sample
p1 <- plotEmbedding(ArchRProj = projMCS6, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

#To plot the UMAP results by cluster
p2 <- plotEmbedding(ArchRProj = projMCS6, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#To visualize these 2 plots side by side
ggAlignPlots(p1, p2, type = "h")

#To save an editable vectorized version of the plots:
plotPDF(p1,p2, name = "UMAP-Sample-Clusters_projMCS6_2024-06-20.pdf", ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

#################################################

##Single cell embeddings - tSNE

#To run a t-Stocastic Neighbor Embedding (tSNE)
projMCS6 <- addTSNE(
        ArchRProj = projMCS6,
        reducedDims = "IterativeLSI",
        name = "TSNE",
        perplexity = 30,
	force = TRUE,
        seed = 123
)

#To plot the tSNE (the same parameters apply to colorBy and name regardless of the type of embedding used)
p1 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "TSNE"
)

p2 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "TSNE"
)

ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "TSNE-Sample-Clusters_projMCS6_2024-06-20.pdf", ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

#######################################
#######################################

##Dimensionality reduction after Harmony for UMAP

#Used to assess the effects of Harmony by visualizing the embedding using UMAP or tSNE and comparing this to the embeddings visualized previously for iterative LSI.

#Repeat the UMAP embedding with the same parameters but for the "Harmony" reduced Dims object
projMCS6 <- addUMAP(
        ArchRProj = projMCS6,
        reducedDims = "Harmony",
        name = "UMAPHarmony",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine",
	force = TRUE,
        seed = 123
)

p3 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "UMAPHarmony"
)

p4 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "UMAPHarmony"
)

ggAlignPlots(p3, p4, type = "h")

plotPDF(p1,p2,p3,p4, name = "UMAP2Harmony-Sample-Clusters_projMCS6_2024-06-20.pdf", ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

##Dimensionality reduction after Harmony for tSNE

#Follow similar steps to those used for UMAP

projMCS6 <- addTSNE(
        ArchRProj = projMCS6,
        reducedDims = "Harmony",
        name = "TSNEHarmony",
        perplexity = 30,
	force = TRUE,
        seed = 123
)

p7 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Sample",
        embedding = "TSNEHarmony"
)

p8 <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "cellColData",
        name = "Clusters",
        embedding = "TSNEHarmony"
)

ggAlignPlots(p7, p8, type = "h")

plotPDF(p1,p2,p7,p8, name = "TSNE2Harmony-Sample-Clusters_projMCS6_2024-06-20.pdf", ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

#########################################
#########################################

##Gene scores and marker genes

# Here we are asking what genes, that exceed thresholds, are different comparing one cluster to all other clusters

Cluster.markers <- getMarkerFeatures(
        ArchRProj = projMCS6,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clusters",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
#        maxCells = 4000
)

#To get a list of dataframe objects, one for each cluster, containing the relevant marker features
clusterMarkerList <- getMarkers(Cluster.markers, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

#Marker list by cluster
clusterMarkerList$C6

for(i in names(clusterMarkerList)) {
        write.csv(clusterMarkerList[[i]], file = paste(i, "_projMCS6_2024-06-20.csv", sep = ""))
}

#########################################

##Heatmaps for marker features

#To visualize all of the marker features simultaneously
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

markerGenesGABA	<- c(
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

heatmapGS <- plotMarkerHeatmap(
        seMarker = Cluster.markers,
        cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
        labelMarkers = markerGenes,
        transpose = TRUE
)

#To plot the heatmap
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#To save the heatmap
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap_projMCS6_2024-06-20.pdf", width = 8, height = 6, ArchRProj = projMCS6, addDOC = FALSE)

#############################################

##Visualizing marker genes on an embedding
markerGenes <- c(
        "Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "GeneScoreMatrix",
        name = markerGenesAstro,
        embedding = "UMAPHarmony",
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
        name = "UMAP-MarkerGenes-WO-Imputation_projMCS6_2024-06-20.pdf",
        ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

#######################################################

##Use MAGIC to impute gene scores by smoothing signal across nearby cells

#Impute weights to the ArchRProject
projMCS6 <- addImputeWeights(projMCS6)

#The impute weights can then be passed to plotEmbedding() when plotting gene scores overlayed on the UMAP embedding
markerGenes <- c(
        "Olig2", "Neurod6", "Gad2"
)

p <- plotEmbedding(
        ArchRProj = projMCS6,
        colorBy = "GeneScoreMatrix",
        name = markerGenesRadialGlia,
        embedding = "UMAPHarmony",
        imputeWeights = getImputeWeights(projMCS6)
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
        name = "UMAP-RadialGliaMarkerGenes-W-Imputation_projMCS6_2024-06-21.pdf",
        ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)

#######################################################
#######################################################

#######################################################
#######################################################

## Highlight cells of interest

################

# Subset samples of interest, add embedding and impute weights

# T1 subset
projMCS6@cellColData$Sample %in% c("C302_", "C306_", "C309_", "C318_", "C323_", "C328_", "C332_", "C337_", "C339_",
        "C346_", "C351_", "C353_", "C360_"),]

# T2 subset
projMCS6@cellColData$Sample %in% c("C304_", "C308_", "C312_", "C315_", "C321_", "C324_", "C327_", "C330_", "C333_",
        "C336_", "C342_", "C348_", "C349_", "C355_", "C358_", "C362_"),]

# T3 subset
projMCS6@cellColData$Sample %in% c("C305_", "C307_", "C313_", "C316_", "C320_", "C322_", "C325_", "C334_", "C340",
        "C341_", "C345_", "C350_", "C352_", "C359_", "C364_"),]

# T4 subset
projMCS6@cellColData$Sample %in% c("C301_", "C303_", "C310_", "C314_", "C319_", "C335_", "C338_", "C344_", "C354_",
        "C356_", "C361_", "C363"),]

#################

p2<- plotEmbedding(
  ArchRProj = projMCS6,
  embedding = "UMAPHarmony",
  colorBy = "cellColData",
  name = "Sample",
  size = 1,
  sampleCells = NULL,
  highlightCells = getCellNames(ArchRProj = projMCS6)[which(projMCS6@cellColData$Sample %in% c("C305_",
        "C307_", "C313_", "C316_", "C320_", "C322_", "C325_", "C334_", "C340",
        "C341_", "C345_", "C350_", "C352_", "C359_", "C364_"))],
  baseSize = 10,
  plotAs = "points")

#To save an editable vectorized version of the plot

plotPDF(p2, name = "UMAP_T3_projMCS6_2024-06--.pdf", width = 5, height = 5, ArchRProj = projMCS6, addDOC = FALSE)

#######################################################
#######################################################

#Save project
saveArchRProject(ArchRProj = projMCS6, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", load = FALSE)
