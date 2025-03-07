#!/usr/bin/env Rscript

#Setup an interactive session
salloc --account=eon -t 0-08:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/Tx_Comp")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

############################################
############################################

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

############################################
############################################

getAvailableMatrices(projMCS5)

table(projMCS5$Clusters)

############################################
############################################

C3ArchRSubset <- projMCS5[projMCS5$Clusters %in% "C3"]

C18ArchRSubset <- projMCS5[projMCS5$Clusters %in% "C18"] 

############################################
############################################

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
treatment <- gsub("C323_", "t3", treatment)
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
C3ArchRSubset$treatment <- treatment

# Check that this worked - if not, make sure the previous line was run successfully
head(C3ArchRSubset$treatment)

##########################################
##  Get marker genes for each cluster:
##########################################


## Get marker features
markerGS <- getMarkerFeatures(
  ArchRProj = C18ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 45000
)

#Get the list of markers
markerList <- getMarkers(markerGS, cutOff = "FDR <= 0.01 & abs(Log2FC >= 1.25)")

#############################

# Write markerList to a CSV file: use for groupBy = Cluster
write.csv(markerList, file = "C18_markerGS_bySample_2024-04-2.csv", row.names = FALSE)

##############################

#Use for groupBy = treatment

unique_treatments <- unique(markerList$treatment)
for (group in unique_treatments) {
  group_markers <- markerList[markerList$treatment == group, ]
  file_name <- paste("markerList_", group, ".csv", sep="")
  write.csv(group_markers, file = file_name, row.names = FALSE)
}

##############################

## To get names of all marker genes identified across all of the clusters:
rowData(markerGS)$name

###

#Make a dataframe of some marker genes 
known_markers_C3 <-c("Olig1", "Mbp", "Opalin", "Spock3", "S100b", "Mag", "Mog", "Cldn11", "Ugt8a")
known_markers_C8 <- c("Olig1", "Gfap", "Aqp4", "Sparc", "Ppp1r3c", "Aldh1l1", "Plcd4", "Mlc1", "Cbs")
known_markers_C10 <- c("C1qa", "Cx3cr1", "Sparc", "Trem2", "Ccl4", "Csf1r", "Cd14", "Tyrobp")
known_markers_C11 <- c("C1qa", "Cx3cr1", "Sparc", "Trem2", "Ccl4", "Csf1r", "Cd14", "Tyrobp", "Cdh5")
known_markers_C16 <- c("Slc17a8", "Tshz2")
known_markers_C18 <- c("Cux2")
known_markers_C18Ex <- c("Wfdc15a", "Ms4a6c", "Wfdc15b", "Dgat2l6", "Fam183b", "Spata31d1a", "Spag11b", "Serpina1f", "Olfr1441", "Olfr464") 
known_markers_C22 <- c("Gad1os", "Gad1", "Gad2", "Grik1", "Reln", "Lhx6", "Pvalb")

known_markers_Top61_C18 <- c("Wfdc15a", "Ms4a6c","Wfdc15b","Dgat2l6","Fam183b","Spata31d1a",
	"Spag11b","Serpina1f","Olfr1441","Olfr464","Mir5107","Panx3","Gm8369","Vmn1r233",
	"Cd163","S100a7a","Trpc5os","Olfr16","Retnlb","Olfr1391","Krt12","Olfr1427","1700009C05Rik",
	"Mir8106","Gm11938","Svs6","Mir215","Gm11937","Lrit2","Olfr317","Tmem213","Defb34","4921517D22Rik",
	"Olfr316","Gzmn","Krtap4-1","Lce1m","Olfr10","Olfr1340","Olfr808","Krtap4-2","Olfr789",
	"Olfr809","Krtap3-2","Olfr56","Krtap1-5","Gm6654","Olfr1428","Slfn10-ps","Krtap3-1",
	"Apol7b","Olfr120","Krtap4-7","Krtap1-4","Olfr1383","Mir216c","Olfr314","Gm11568","Gm11559",
	"Svs2","Rhox4a"
)

known_markers_All127_C18 <- c("Wfdc15a", "Ms4a6c","Wfdc15b","Dgat2l6","Fam183b","Spata31d1a",
        "Spag11b","Serpina1f","Olfr1441","Olfr464","Mir5107","Panx3","Gm8369","Vmn1r233",
        "Cd163","S100a7a","Trpc5os","Olfr16","Retnlb","Olfr1391","Krt12","Olfr1427","1700009C05Rik",
        "Mir8106","Gm11938","Svs6","Mir215","Gm11937","Lrit2","Olfr317","Tmem213","Defb34","4921517D22Rik",
        "Olfr316","Gzmn","Krtap4-1","Lce1m","Olfr10","Olfr1340","Olfr808","Krtap4-2","Olfr789",
        "Olfr809","Krtap3-2","Olfr56","Krtap1-5","Gm6654","Olfr1428","Slfn10-ps","Krtap3-1",
        "Apol7b","Olfr120","Krtap4-7","Krtap1-4","Olfr1383","Mir216c","Olfr314","Gm11568","Gm11559",
        "Svs2","Rhox4a","Ptgs2","4930452N14Rik","4930532M18Rik","Samd3","Themis","Glt8d2","Olfr1392",
	"Car4","Kcnh4","Crhr1","Sntg2","Slc24a4","Serpina1e","Eef1e1","A430090L17Rik","Gm8267","Nrl",
	"Mir6239","1700006F04Rik","Agxt2","Kcnv1","Sstr3","Nptxr","Rtn4r","Mir802","Dact2","Smoc2",
	"Psors1c2","Cdsn","4933400B14Rik","Gm4719","Sh3rf2","Cerkl","Gm13710","Rtn4rl2","Itpka",
	"Mir1952","Ovol2","Arsj","Clca3a2","Calb1","Tmem215","4933428C19Rik","Mpl","Tmco2","Rspo1",
	"Vwa5b1","Ubxn10","Pla2g2c","Slc30a3","Trim54","Cux2","Ocm","Neurod6","Unc5d","Prss54",
	"Snai3","Vsig2","Olfr877","Olfr145","Exph5","Sh2d7","Lingo1","Uchl4","Grm2","Ackr2"
)

known_markers_Top200_C3 <- c("S100b","Snx1","Picalm","S100a16","Trim13","Art4","Rpl10l","Cdk5rap2",
	"Mapk8ip1","Bche","Gm10400","Trpm1","Cep97","Card10","Pcbp4","Cdk18","Sox8","Smad7","Nkx6-2",
	"Matn4","Ankrd13a","4930581F22Rik","Mir365-2","1810010D01Rik","Syne3","Ccdc83","Bcar1",
	"1700056E22Rik","Elovl6","Tbc1d16","4931419H13Rik","Chdh","9130209A04Rik","Itgb4","Spdl1",
	"Tspan32","Slc15a4","Gm21671","Il1rap","Fzd7","Cdc42se2","AW046200","Dnajc7","Pls1","Col2a1",
	"Gm9733","Chmp7","Lifr","Pstpip2","Abtb2","Pcdh17","Arrdc2","4931406P16Rik","Mid1ip1",
	"4930413E15Rik","Dock2","R3hcc1","Gdpd1","Stambp","Crybg2","Hid1","Qdpr","Kcnrg","4930544G11Rik",
	"Tcf12","Mettl16","Ifnk","Tmem196","2310010J17Rik","Flnc","Gm31641","Rnf122","Gm9895","Serinc5",
	"Atp1b3","Ppp1r17","Ano1","Itih5","Tppp","Rgma","Foxi2","Fgd3","Gstp1","Pcolce2","Lhfpl2","Polm",
	"Lbh","Btla","Prss48","Enpp2","Trp53bp2","Nat8f4","4930426L09Rik","Tma16","Ccdc13","Man2a2","Pigz",
	"Tmem125","D7Ertd443e","Pou3f1","Lcp1","Lpo","Ick","Adamts2","BC055402","Cryaa","Pla2g3","Mir1961",
	"Gng7","Nipal4","Lbr","45359","B230344G16Rik","Zdhhc20","4930515G01Rik","Svip","2700070H01Rik",
	"Cep72","Sppl2c","Mtmr2","Gstm7","Ippk","Pcyt1a","Ndst4","Spag16","Gm4668","Mob3b","Tnfaip6",
	"9230105E05Rik","Aif1l","Slc20a2","Prdm5","Trmt13","Sema5a","Osbpl7","Ccp110","Zfp276","2210409D07Rik",
	"Iqck","Ttll7","Cdc42ep2","Rnf141","Chst8","Rps6ka1","Gas2","Pde1c","Commd9","Glt6d1","Slc25a48",
	"Polr2f","BC065397","Hdac11","Actbl2","4930554C24Rik","Utp11","Klf15","Tjap1","Clrn3","2610528J11Rik",
	"Olfml1","Mecom","Plekha2","2010010A06Rik","Gpr62","Halr1","Gsn","Bet1","Vldlr","Fgf7","Kat2b",
	"Mink1","Sh3gl3","Sec11c","Gm15417","Cnn3","Ripk4","Pdlim5","Slc5a11","Smim5","Pick1","Tmprss12",
	"Zbtb7b","Hcn2","Lrig3","4933417E11Rik","Arhgap31","Defb30","Gm15713","St18","Emilin2","Fam71e2",
	"Nrap","Mir6896","Gm12530","Capn2","Lin28a","Zfp42","Tfpi2","A930017K11Rik","Ankib1"
)


###

#Subset the markerGS object to just these known markers
se_idx <- which(rowData(markerGS)$name %in% known_markers_All127_C18)  

#Get an index of these from the summarized experiment
subset_markerGS <- markerGS[se_idx,]  

##Plot it out:

pdf(file="C18_All127_TxComp_Heatmap_2024-03-13.pdf", width=15, height=30)
plotMarkerHeatmap(   ### Copied this over from the full function reference - leaving anything not commented at defaults
  seMarker = subset_markerGS,   ## Set this to the correct object containing the markers
  cutOff = "FDR <= 1 & Log2FC >=0.01",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = known_markers_All127_C18,  #Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)
dev.off()

######################################

# To print the Z-scores for the above heatmap

z_scores <- plotMarkerHeatmap(
  seMarker = subset_markerGS,
  cutOff = "FDR <= 1 & Log2FC >=.01",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = known_markers_All127_C18,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  # Change this to TRUE
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

# Print z-scores from plotMarkerHeatmap
write.csv(z_scores, "C18_All127_byTx_zscores_2024-03-13.csv", row.names = TRUE)


###

#C18 & C18Ex: Error in plotMarkerHeatmap(seMarker = subset_markerGS, cutOff = "FDR <= 1 & Log2FC >=0.001",  : 
#  No Makers Found!


######################################
######################################
######################################
######################################

## Could not run the following: 

#Option 2:

        mark_test <- getMarkerFeatures( # get the marker features for the comparison
          ArchRProj = C8ArchRSubset,  ## Here, use the subsetted project for only the cluster of interest
          useMatrix = "GeneScoreMatrix",
          groupBy = "treatment",   
          testMethod = "wilcoxon",
          bias = c("TSSEnrichment", "log10(nFrags)"),
          useGroups = "t3",
          bgdGroups = "t1", 
          maxCells = 8000)
        
#        markers <- getMarkers(mark_test, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25") # get the markers
        

#####    Try making a whole heatmap and then subsetting it
##############################################################################

c8_known_markers <- c("Sparc")

hmap_to_subset <- plotMarkerHeatmap(  #Anything not commented left as default
  seMarker = mark_test,   ##Set this to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = TRUE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = c8_known_markers,  ##Label the markers of interest to plot
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap_2024-02-29")
)

## Get the matrix of the heatmap to figure out row indices for the genes of interest

hmap_to_subset_MAT <- plotMarkerHeatmap(   #Anything not commented left as default
  seMarker = markerGS,   ##Set to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = known_markers,  ##Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  #Return the matrix rather than a heatmap
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)


hmap_row_idx <- which(rownames(hmap_to_subset_MAT) %in% known_markers) #Get the row numbers of the heatmap matrix that match the known markers

#Plot it out
pdf(file="hmap_excitatory_subset.pdf", width=7, height=5)
hmap_to_subset[hmap_row_idx, ]
dev.off()




####################################
###################################

        # make a volcano plot:
        volc <- plotMarkers(seMarker = mark_test, name = colnames(mark_test), cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25", plotAs = "Volcano")
        
# if(nrow(markers[[1]]) > 0){
#          df_labeled <- cbind(comparison_name, markers[[1]])
#          samp_comps_labs <- c(samp_comps_labs, df_labeled)
#          names(samp_comps_labs)[[length(samp_comps_labs)]] <- comparison_name # put a name on each comparison - here we index by length(samp_comps_labs) to put the name on the last element of the list, because it no longer matches i 
 
         #Save a pdf of the Volcano plot
          pdf(file=paste0(j, "_", comparison_name, "_Volcano_2024-02-29.pdf"), width=5, height=6)
          print(volc)
          dev.off()
#        }else{
#          empty_comps <- c(empty_comps, comparison_name)
#        }
#      },
#      error=function(cond){
#        errors <- paste(comparison_name, "failed. Error is: ", cond)
#        return(errors)
#      },
#      warning=function(cond){
#        warns <- paste(comparison_name, "Had warning warning is: ", cond)
#        return(warns)
#      }
#    )
#  }

#if(length(empty_comps) > 0){
#    write.csv(empty_comps, file=paste0(j, "_EMPTY_sample_comps.csv"))
#  }
#  if(length(error_comps) > 0){
#    write.csv(error_comps, file=paste0(j, "_ERRORS_sample_comps.csv"))
#  }
  
  
#  if(length(samp_comps_labs) > 0){
#    # Combine the ones that aren't empty
#    combined_sample_comps <- do.call("rbind", samp_comps_labs)
#    write.csv(combined_sample_comps, file=paste0(j, "_non_empty_sample_comps_2024-02-29.csv"))

######################################
#####################################    
    
    #New heatmap using ggplot
    best_df1 <-as.data.frame(combined_sample_comps)
    
    # Split out things by which sample is first
    split_comps <- strsplit(best_df1$comparison_name, split = "_")
    first_samp <- sapply(split_comps, function(x) x[[1]])
    best_df <- cbind(best_df1, first_samp)
    
    best_df$first_samp <- factor(best_df$first_samp, levels = c("C337", "C353", "C346", "C340", "C345", "C336", "C355", "C348", "C335", "C356", "C344"))
    
    p <- ggplot(markers, aes(x = name, y = comparison_name, fill = Log2FC)) +
      geom_tile(color = "black", lwd = 0.2, linetype = 1) +
      scale_fill_gradient2(low = "blue",  mid = "white", high = "red") +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_grid(cols = vars(seqnames), rows = vars(first_samp), scales = "free", space = "free",) +
      ggtitle(paste0(j, " heatmap")) +
      xlab("Gene name") + ylab("Samples")

######################################
#####################################   

# Get some colors:

make_cols1 <- sort(unique(C8ArchRSubset$treatment))
make_cols1 <- gsub("C302_", "t5", make_cols1)
make_cols1 <- gsub("C306_", "t5", make_cols1)
make_cols1 <- gsub("C309_", "t5", make_cols1)
make_cols1 <- gsub("C318_", "t5", make_cols1)
make_cols1 <- gsub("C328_", "t5", make_cols1)
make_cols1 <- gsub("C332_", "t5", make_cols1)
make_cols1 <- gsub("C337_", "t5", make_cols1)
make_cols1 <- gsub("C339_", "t5", make_cols1)
make_cols1 <- gsub("C346_", "t5", make_cols1)
make_cols1 <- gsub("C351_", "t5", make_cols1)
make_cols1 <- gsub("C353_", "t5", make_cols1)
make_cols1 <- gsub("C360_", "t5", make_cols1)
make_cols1 <- gsub("C304_", "t5", make_cols1)
make_cols1 <- gsub("C308_", "t5", make_cols1)
make_cols1 <- gsub("C312_", "t5", make_cols1)
make_cols1 <- gsub("C349_", "t5", make_cols1)
make_cols1 <- gsub("C315_", "t5", make_cols1)
make_cols1 <- gsub("C321_", "t5", make_cols1)
make_cols1 <- gsub("C324_", "t5", make_cols1)
make_cols1 <- gsub("C355_", "t5", make_cols1)
make_cols1 <- gsub("C327_", "t5", make_cols1)
make_cols1 <- gsub("C330_", "t5", make_cols1)
make_cols1 <- gsub("C333_", "t5", make_cols1)
make_cols1 <- gsub("C358_", "t5", make_cols1)
make_cols1 <- gsub("C336_", "t5", make_cols1)
make_cols1 <- gsub("C342_", "t5", make_cols1)
make_cols1 <- gsub("C348_", "t5", make_cols1)
make_cols1 <- gsub("C362_", "t5", make_cols1)
make_cols1 <- gsub("C305_", "t6", make_cols1)
make_cols1 <- gsub("C307_", "t6", make_cols1)
make_cols1 <- gsub("C313_", "t6", make_cols1)
make_cols1 <- gsub("C350_", "t6", make_cols1)
make_cols1 <- gsub("C316_", "t6", make_cols1)
make_cols1 <- gsub("C320_", "t6", make_cols1)
make_cols1 <- gsub("C322_", "t6", make_cols1)
make_cols1 <- gsub("C352_", "t6", make_cols1)
make_cols1 <- gsub("C323_", "t6", make_cols1)
make_cols1 <- gsub("C325_", "t6", make_cols1)
make_cols1 <- gsub("C334_", "t6", make_cols1)
make_cols1 <- gsub("C359_", "t6", make_cols1)
make_cols1 <- gsub("C340_", "t6", make_cols1)
make_cols1 <- gsub("C341_", "t6", make_cols1)
make_cols1 <- gsub("C345_", "t6", make_cols1)
make_cols1 <- gsub("C364_", "t6", make_cols1)
make_cols1 <- gsub("C301_", "t6", make_cols1)
make_cols1 <- gsub("C303_", "t6", make_cols1)
make_cols1 <- gsub("C310_", "t6", make_cols1)
make_cols1 <- gsub("C314_", "t6", make_cols1)
make_cols1 <- gsub("C319_", "t6", make_cols1)
make_cols1 <- gsub("C335_", "t6", make_cols1)
make_cols1 <- gsub("C338_", "t6", make_cols1)
make_cols1 <- gsub("C344_", "t6", make_cols1)
make_cols1 <- gsub("C354_", "t6", make_cols1)
make_cols1 <- gsub("C356_", "t6", make_cols1)
make_cols1 <- gsub("C361_", "t6", make_cols1)
make_cols1 <- gsub("C363_", "t6", make_cols1)

  
    # convert to colors:
    make_cols2 <- gsub("t1", "blue", make_cols1)
    make_cols2 <- gsub("t2", "red", make_cols2)
    make_cols2 <- gsub("t3", "yellow", make_cols2)
    make_cols2 <- gsub("t4", "pink", make_cols2)
    
    
    g <- ggplot_gtable(ggplot_build(p))
    stripr <- which(grepl('strip-r', g$layout$name))
    fills <- make_cols2
    k <- 1
    for (i in stripr) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    pdf(file=paste0(j, "_ggplot_Test_2024-02-29.pdf"), width=length(unique(best_df$name))/2.7, height=length(unique(best_df$comparison_name))/4)
    print(grid::grid.draw(g))
    dev.off()
  }else{
    error_comps <- c(error_comps, paste0("cluster ", j, " identified no marker genes for any individual comparison"))
    write.csv(error_comps, file=paste0(j, "_ERRORS_sample_comps_2024-02-29.csv"))
  }
}









##########################################
##########################################
##########################################
##########################################

