Setup an interactive session
salloc --account=eon -t 0-05:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

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
setwd("/project/eon/nboldon/MCS2023/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS6)
table(projMCS6$Clusters)

############

projMCS6 <- addHarmony(
        ArchRProj = projMCS6,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
        force = TRUE
)

############################################
############################################

C1ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C1"]
C2ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C2"]
C3ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C3"]
C4ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C4"]
C5ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C5"]
C6ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C6"]
C7ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C7"]
C8ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C8"]
C9ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C9"]
C10ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C10"]
C11ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C11"]
C12ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C12"]
C13ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C13"]
C14ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C14"]
C15ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C15"]
C16ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C16"]
C17ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C17"]
C18ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C18"]
C19ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C19"]
C20ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C20"]
C21ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C21"]
C22ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C22"]
C23ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C23"]
C24ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C24"]
C25ArchRSubset <- projMCS6[projMCS6$Clusters %in% "C25"]

#Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!

############################################
############################################

# Specify which treatment group each sample is in:

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

treatment <- C25ArchRSubset$Sample

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
C25ArchRSubset$treatment <- treatment

# Check that this worked - if not, make sure the previous line was run successfully
head(C25ArchRSubset$treatment)

############

C25ArchRSubset <- addHarmony(
        ArchRProj = C25ArchRSubset,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "treatment",
        force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
C25ArchRSubset <- addImputeWeights(C25ArchRSubset)
getImputeWeights(C25ArchRSubset)

############################################

##########################################
##  Get marker genes for each treatment:
##########################################


## Get marker features
c25Markers <- getMarkerFeatures(
  ArchRProj = C25ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  #maxCells = 45000
)

#Get the list of markers
c25markerList <- getMarkers(c25Markers, cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")

#cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25"
#cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1.15

#############################

# Write markerList to a CSV file: 
write.csv(c25markerList, file = "C25_TxComp_FDR-0-1_Log2FC-0-5_2024-06-26.csv", row.names = FALSE)


##########################################################
##########################################################
##########################################################

# Using ArchRProj = projMCS5, groupBy = Sample
## markerList$C1 = all clusters returned NULL values

write.csv(markerList, file = "GeneMarker_List_Master_2024-04-26.csv")
write.csv(markerList$C1, file = "C1_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C2, file = "C2_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C3, file = "C3_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C4, file = "C4_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C5, file = "C5_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C6, file = "C6_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C7, file = "C7_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C8, file = "C8_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C9, file = "C9_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C10, file = "C10_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C11, file = "C11_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C12, file = "C12_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C13, file = "C13_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C14, file = "C14_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C15, file = "C15_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C16, file = "C16_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C17, file = "C17_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C18, file = "C18_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C19, file = "C19_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C20, file = "C20_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C21, file = "C21_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C22, file = "C22_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C23, file = "C23_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C24, file = "C24_GeneMarker_List_2024-04-26.csv")
write.csv(markerList$C25, file = "C25_GeneMarker_List_2024-04-26.csv")

##############################

## To get names of all marker genes identified across all of the clusters:
rowData(markerGS)$name

###

#Make a dataframe of some marker genes 

#C1 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C1 <-c("Nat8f6", "H2-Ea-ps", "Cts8-ps", "Olfr699", "Olfr663", "E330017A01Rik",
	"Olfr829", "Vmn2r13", "Cts3", "Mrgprb2", "Glt6d1", "Scarletltr", "Cts7", "Hoxb6", "Mkrn3",
	"Nlrp2", "Cmtm5", "Vmn2r111", "Rnase1", "Mrgprb3", "Vmn2r61", "Haglr", "Mug1", "Lyg1",
	"Il25", "Gm5168", "Cltrn", "Cyp2d10", "Mir5135", "Defa21", "4930544G11Rik", "Vmn2r70",
	"Mettl11b", "I730030J21Rik", "Gm595", "S100a1", "Cyp4f17", "Nat8", "Halr1", "Foxn4", "Vmn2r32",
	"Ccdc70", "Lpo", "Nat8f2", "Dcxr", "Gm6083", "Ermn", "Tdpoz3", "Scgb2b21", "Nat8f5", "Adam18",
	"Mag", "Rnase6", "Efcab3", "Usp17le", "Mir6236", "Gm6093", "Tmem81", "Vmn2r25", "Arhgef16", 
	"Gstm6", "Eln", "5930403L14Rik", "Vmn1r174", "Defb30", "Mir6924", "Pmel", "Mrgprb4", "Ttyh2",
	"Anln", "Hoxb5", "D630033O11Rik", "Fbxl21", "Vmn2r71", "Pgk2", "Fam107a", "Cd9", "Atp6v1e2",
	"Selenop", "Serpinb1a", "Cyp2d9", "Lcn3", "B230206H07Rik", "Fbxw19", "Cyp4f15", "Gm4632",
	"1190005I06Rik", "Taar9", "Apoc1", "Bpifc", "BC064078", "Vim", "Olfr24", "Micall2", 
	"2700070H01Rik", "Mrgpra3", "Vmn1r78", "Folh1", "Pax8", "Acod1", "Pp2d1", "Gstm3", "Ugt1a2",
	"Ugt1a6a", "Vmn2r74", "Ugt1a7c", "Fa2h", "Padi2", "Als2cl", "Mir8110", "Ugt1a10", "Pdk4",
	"Car3", "Fzd2", "Vmn1r77", "Vmn2r78", "Pmp2", "Klhdc7a", "Tmem82", "Vmn2r104", "A530053G22Rik",
	"Tmem95", "Itih3", "Hhipl2", "Vmn1r167", "Fga", "4930447J18Rik", "Asap3", "Ankrd66", "Celsr1",
	"1700003P14Rik", "4930432M17Rik", "Vmn2r117", "Tead3", "Cyp2j9", "Bpifa6", "Gm3219", "Endou",
	"2810047C21Rik1", "Tktl2", "Ugt1a9", "Plac1", "Mir511", "Plin1", "Gjb6", "Gm10863", "Gbp5",
	"Mir7662", "Gm20219", "Spc25", "Kcnj10", "Anxa8", "Hoxb5os", "Nutm2", "Mir683-1", "AY761185",
	"5031434C07Rik", "Kirrel2", "Vmn2r45", "Esp23", "Klkb1", "Prdm16", "Castor1", "4732490B19Rik",
	"Cryaa", "Fbxw15", "4930550C17Rik", "Ankub1", "Cyp2j5", "Cxcl14", "Mrgprb8", "Apoe", "Mir7001",
	"Mir7229", "Mir344f", "Rgr", "Fgfr2", "4930425K10Rik", "Gm1979", "Vmn2r43", "Daam2", "Rasl12",
	"Ppp1r14a", "Eci1", "Mug2", "Plpp3", "Sox10", "Gm13293", "Olig1", "9230117E06Rik", "Lgals2",
	"4933431E20Rik", "Rpl10l", "Dmrta2", "3930402G23Rik", "Hoxa4", "Eva1a", "4930515L03Rik",
	"Gm28453", "Plxnb3")

#C2 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C2 <-c("Cd300c2", "Cd33", "Clec1b", "Cd22", "Olfr1193", "Mkrn3", "Tmem125", "Peg12",
	"Olfr498", "Clec4a3", "Gpr160", "Vmn2r77", "Gm6093", "Vmn2r100", "Smim5", "Zc3hav1",
	"4933406J10Rik", "Tmem119", "Fga", "Dennd1c", "Fgl2", "Esp24", "Myoz2", "2210414B05Rik",
	"Mir5135", "Gjc3", "Tex45", "Gm12530", "Dyrk4", "Cpt2", "Bcas1", "Cx3cr1", "Bin2", "Gpr17",
	"Smtnl2", "Mfng", "Mir6961", "Rnase1", "Cyp3a13", "Speer3", "Has2os", "Gng11", "Vmn2r97",
	"Gm6607", "Nkx2-2", "Cdh1", "Hfe", "C1qc", "Mir19a", "Mir18", "Mir20a", "Mir6370", "Ccr1",
	"4930591E09Rik", "Tgif1", "Mir344i", "Bcas1os2", "Sox10", "Zfp773", "Lims2", "Bmp4", "Pp2d1",
	"Nox4", "Siglecg", "Mir351", "Dcpp2", "Vipr2", "Tyw5", "Sh3bgrl", "Tagln2", "Tmem102", "Mag",
	"Siglecf", "Mir503", "Ctsc", "Il18", "C430049B03Rik", "Klk9", "Csf1r", "Gm5547", "Prps1l1",
	"Opn4", "Arhgap11a", "Cd9", "Snx33", "Olig1", "Bfsp2", "Mir145b", "Bcl3", "Cd300lf",
	"Sap30bp", "Trpv4", "1700024I08Rik", "2810429I04Rik", "Vsir", "Asap3", "Ccr2", "Ssb", "Agtr1b",
	"0610039K10Rik", "Fgf7", "Ggct", "Slc7a7", "Mir344f", "Gjc2", "Fam71a", "Tlr1", "Dach2", "Heph",
	"Prkcq", "Cdh19", "S100a1", "Csf3r", "Gng10", "Zfp474", "Ctsh", "Fam19a4", "1700021F07Rik",
	"Stk17b", "Zfp36", "A2ml1", "Pla2g15", "Gm33104", "Dynlt3", "Nrap", "Enpp1", "Hhex", "A530053G22Rik",
	"Fa2h", "Pdgfra", "Rhoh", "Xkr5", "Prkd3", "Plekhg3", "Pdlim5", "Mal", "Plp1", "Car2", "Cyp2j9",
	"Tns3", "Cfi", "Recql5", "E230029C05Rik", "Kdm1b", "Lgals2", "Tspan2", "Mir17hg", "Cp", "Has2",
	"Hhip", "Cyp2j11", "Mmp2", "Pon1", "Gm5126", "Timm8a2", "Tifab", "4930563H07Rik", "Echs1",
	"Gm10863", "Lncbate1", "Mir7001", "Tspan2os", "Carhsp1", "Slc25a27", "Lhfpl2", "Sema3d", "Zfp658",
	"Nat8f4", "Cpt1a", "Ctsd", "Olig2", "2210409D07Rik", "Tent5a", "Tgfb1", "AI467606", "Efemp1",
	"4833418N02Rik", "Cers2", "Dll1", "2700070H01Rik", "Slc5a11", "Cryaa", "Kctd11", "Plekhg1",
	"Frmd4b", "Cela1", "Dnah17", "Lgmn", "Tmem163", "Urm1", "2310010J17Rik", "Fli1", "Nfe2l2",
	"Lrig3", "6430503K07Rik", "Galnt4", "BC024139", "Rab37", "Tgm2", "Ccr6")

#C3 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C3 <-c("Opalin", "Gjc3", "Mag", "2310005A03Rik", "Ppp1r14a", "Mal", "Glt6d1", "Mir219a-2",
	"Mir7001", "Fbxw15", "Acsm3", "Gm10863", "Gng11", "5031410I06Rik", "Ermn", "Aspa", "Gm5862",
	"Mir6917", "Pkd2l1", "Fga", "A230009B12Rik", "Fa2h", "Sox10", "Hapln2", "Gm10471", "Defb30",
	"Folh1", "Gm3402", "Mir5135", "Gm1979", "Erbb3", "Trf", "Rnase1", "Ankub1", "Plp1", "Sgk2",
	"Speer4a", "Fbxw19", "Mir219b", "Kcnrg", "Rnf220", "Cldn11", "Olfr1331", "Tmem63a", "C030029H02Rik",
	"Defa32", "A530053G22Rik", "2700070H01Rik", "4930544G11Rik", "Car14", "Pp2d1", "Cd82", "Zfp488",
	"Gjc2", "Gm13293", "Plxnb3", "Tshb", "Cyp2j12", "Mir7239", "Gpr17", "Tmprss12", "5430425K12Rik",
	"Lgals2", "Zfp36l3", "Smco3", "Opn4", "Gngt1", "B230206H07Rik", "Galnt6", "Scp2d1", "Ptgds",
	"Efhd1", "Pde8a", "Fgfr2", "Olfr52", "Thsd1", "4930405J17Rik", "Mir7230", "9330199G10Rik",
	"Nkain1", "Olfml1", "B230344G16Rik", "Rpl10l", "Myrf", "Tmem125", "Thumpd1", "Vmn2r77", "Ugt8a",
	"Nkx2-2", "Plin4", "Bpifa6", "Nmral1", "Prr18", "Vmn1r176", "Anln", "Carhsp1", "Lpo", "Cnp", "Lims2",
	"1500015L24Rik", "Tfpi2", "Cdc42ep1", "Sp7", "Tas2r105", "Tspan2", "Cmtm5", "Snx33", "Gpr37",
	"Vmn1r2", "Fgf7", "Olig1", "Gm21671", "1300017J02Rik", "Bcas1os2", "Art4", "S1pr5", "BC024139",
	"Chst4", "AI463170", "5430431A17Rik", "Efnb3", "Olfr50", "Prima1", "Lect2", "Ggt6", "Plekhg3",
	"Defb47", "Trim13", "Cltrn", "Tubb4a", "4930502E18Rik", "Sec14l5", "Gm12530", "Fam71e2", "Mir8115",
	"Trpv3", "Eda2r", "Mucl1", "Ace2", "4930449I24Rik", "Nkx2-2os", "Ggct", "Ppfibp2", "Prkd3", "Pola2",
	"Bcas1", "Ttyh2", "Pla2g16", "Olfr1333", "Rffl", "Elovl7", "Azgp1", "Cd9", "7630403G23Rik", "Plekhf1",
	"1700044K03Rik", "Kcnk13", "Krt40", "Unc5b", "Smim5", "Cryaa", "Cdh19", "Mir6909", "Vmn2r78", "Aaas",
	"Olfr1023", "Olfr1039", "Trim59", "Il25", "1700081H04Rik", "Rnase11", "Acod1", "Klk6", "Hspb2", "Rhog",
	"Mir7243", "Hhip", "Tmem258", "Lrrc63", "S100a16", "A930017K11Rik", "Flnc", "Sapcd1", "Olfr1043",
	"Abca8a", "Gab1", "Olfr503", "Tspan2os", "Olfr204", "6430503K07Rik", "Nkx6-2", "Gm9895", "Prkcq",
	"Dcpp1", "4732490B19Rik", "Gm10804", "2810429I04Rik", "Rpl32l", "Gal3st1", "Olfr1042")

#C4 markers: FDR <= 0.05 & Log2FC >=1.15, all 75 sorted by Log2FC
known_markers_C4 <-c("Olfr1040", "Olfr906", "Ube2dnl2", "Mir188", "Olfr1184", "Olfr715b", "Olfr1214", "Olfr503",
	"Olfr1181", "Olfr490", "Hbb-bs", "Olfr486", "Mir362", "Gm3402", "Olfr472", "Gm1979", "Rnase11",
	"Vmn1r204", "Olfr1024", "Mir501", "4930449I24Rik", "Olfr1086", "Olfr1232", "Gm3286", "2310005A03Rik",
	"Olfr479", "Vmn2r51", "Olfr487", "Neil3", "D5Ertd577e", "Tmem205", "Vmn2r45", "Gng11", "Vmn2r98",
	"Mucl1", "Cyp2j12", "Pilrb2", "Samt3", "Cyp2j11", "Mir5135", "Olfr1118", "Olfr186", "Olfr183",
	"Sult2a5", "Ermn", "Cldn11", "Prss32", "Trf", "A530053G22Rik", "Rnase1", "Zscan4d", "Folh1", "Gm6370",
	"Serpinb6b", "5430425K12Rik", "Ptgds", "1500015L24Rik", "Gjc3", "Cyp2j8", "7630403G23Rik", "Lgals2",
	"Pet2", "Ppp1r14a", "Sult2a3", "Fa2h", "Fam71e2", "Klra4", "Vmn2r72", "Vmn2r121", "Olfr1200", "Fcrl6",
	"Gm7978", "Olfr129", "Vmn2r78", "Rpl10l")

#C5 markers: FDR <= 0.01 & Log2FC >=1.25, all 102 sorted by Log2FC
known_markers_C5 <-c("Spaar", "Gm13293", "Ppp1r14a", "Zic1", "Pald1", "2310005E17Rik", "Rbpms", "Ror1",
	"Slc6a13", "Zic4", "2700069I18Rik", "Atf3", "Trpv3", "Foxd1", "Card10", "Snx33", "Fzd6", "Xdh",
	"1700083H02Rik", "Gm35978", "Ebf1", "Eya2", "Htra3", "Pdlim2", "Gypc", "Vangl1", "Nr2f2", "Afap1l2",
	"Ano1", "Unc5b", "Agmo", "Snai1", "4930581F22Rik", "Egfl7", "Slc16a12", "Mc2r", "1810010D01Rik",
	"Id3", "Clic4", "Lrrk1", "Cgnl1", "Flt1", "Zfhx4", "Notch1", "Tmc7", "Heg1", "Cdc42ep1", 
	"4933440J02Rik", "Rnf7", "Sik1", "Rftn1", "Adgra2", "Kirrel", "Rnf220", "Cd9", "Ttyh2", "Tgfbi",
	"Frmd4b", "Dixdc1", "Amotl2", "Cdh19", "Lbh", "Mylk4", "Sorbs3", "Plce1", "Sept4", "Adamtsl4",
	"Prkcq", "Mal", "Tanc1", "Adamtsl1", "8430430B14Rik", "Fam178b", "Nkain1", "Cryl1", "Sema3d",
	"Gm10863", "Gm15713", "Mob3b", "Klf3", "Plekhg3", "Cald1", "Tns1", "Qk", "Pld1", "Pola2", "Epas1",
	"Tfeb", "Abca4", "Cdc42", "Foxn3", "Hsf2bp", "Rbms1", "Gab1", "Inpp5a", "A630001O12Rik", "Neat1",
	"A230009B12Rik", "Ets1", "Colec12", "Fgfr2", "4933433H22Rik")

#C6 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C6 <-c("Hba-x", "Olfr315", "Dcpp1", "Smok3a", "Klk9", "4930504O13Rik", "Dcpp2", "Pdgfra",
	"Olfr556", "Klk8", "Ifna4", "Olfr312", "C1ql1", "Sun5", "Gpr17", "Mir7025", "BC048671", "Cacng4",
	"Gm7903", "Clba1", "Olfr30", "Mif4gd", "Doxl2", "Lims2", "Gsx1", "Neu4", "B3gnt7", "Tat", "Mir5621",
	"Rbpjl", "Pnlip", "Olig2", "Olfr1348", "Gm4632", "AA545190", "Fam89a", "Prl2a1", "Myt1", "Slfn1",
	"Matn4", "Prss21", "Sox10", "Irx6", "4930563H07Rik", "Has2os", "Prl3d2", "Ntn1", "Kank1", "Gjc3",
	"Acan", "H2al3", "Traf4", "Prl7b1", "Tmem255b", "Has2", "Gm26633", "Gphb5", "Gm2176", "Mir7684",
	"A930009A15Rik", "Olig1", "Zfp488", "Kif18b", "Prl8a1", "Olfr341", "Mir6406", "Tmem176a", "Mir759",
	"Olfr1537", "Caskin2", "Fnd3c2", "Cspg4", "Vwa2", "Kcnj16", "Cdrt4", "Gpnmb", "Inava", "Mir542",
	"2310002F09Rik", "Slc36a3", "Chst5", "Gm5091", "Krt5", "Afap1l2", "Rlbp1", "Tcl1b4", "Tmem100", "Klk7",
	"S100a16", "Cyp2d37-ps", "Olfr328", "Mkrn3", "Gm10863", "Mir450-1", "Olfr1453", "Xkr5", "Tagln2",
	"Cmtm5", "Mir450-2", "Snord68", "Kif12", "Oas1h", "Mir7686", "5031434C07Rik", "Krtap21-1", "Adora2b",
	"Rpl13", "Tmem176b", "Il25", "Mmp15", "Tgfa", "Bcas1os2", "Gm13290", "Fzd9", "Epn2", "Scrg1", "Mir6996",
	"Gzmc", "Sdr42e1", "1700047L14Rik", "2610528J11Rik", "Gpr37l1", "Prl7c1", "Susd5", "Megf11", "Svs1",
	"P2ry1", "Olfr324", "Gm15326", "Cyp2d10", "1600015I10Rik", "Mir351", "Olfr12", "Lrit1", "Gm3402",
	"Stk32a", "Prss50", "Azgp1", "Mir145b", "Mir503", "Olfr157", "Gm3143", "Mrgprd", "Tstd1", "Sox3", 
	"Tac1", "1700057H15Rik", "4930486F22Rik", "1700044C05Rik", "C430049B03Rik", "Notch1", "Vmn1r50",
	"Cyp2d13", "Sh3bp4", "Lcat", "Pcdhb2", "Gm6961", "Lad1", "Plac1", "Fam240a", "Meox1", "Fgfrl1",
	"Ccl7", "Sema5b", "Myo7b", "Foxl1", "Spg7", "Lhfpl3", "E2f8", "Mir7074", "Halr1", "Ear14", 
	"6430503K07Rik", "C1ql2", "Ccnd1", "Sapcd1", "Slc2a10", "Slc12a4", "Mir322", "Fam19a4", "Unc93a2",
	"Ncmap", "Dcaf12l2", "Ptprz1", "Pmel", "Vipr2", "Col5a3", "Nkx2-2os", "Unc13c", "Plut", "Rnf43",
	"Camp", "Mannr", "Cdo1", "Pabpc6", "4930444F02Rik", "9230117E06Rik", "Bcas1", "H2-M10.3", "Ednrb")

#C7 markers: FDR <= 0.01 & Log2FC >=1.25, all 119 sorted by Log2FC
known_markers_C7 <-c("Fzd9", "Acsf2", "Slc30a10", "Wfdc1", "Chad", "Zic4", "Tnfsf12", "Mir26b", "Mir3966",
	"Jag1", "Slc7a2", "Twist1", "Phka1", "Heyl", "Foxc1", "Nkd2", "Flt1", "Gm36595", "Colec12", "Cd9",
	"Snai1", "Ppara", "Traf4", "Tmem51", "Foxd1", "Lef1", "Kirrel", "B130024G19Rik", "Slc38a3", "Gulp1",
	"1700018B08Rik", "Slc6a13", "Slc9a3r1", "Sorbs3", "Lrrk1", "1700060C20Rik", "Adgra2", "2310005E17Rik",
	"Zic1", "Nr2f2", "Atp1a2", "Bach1", "Vstm4", "Bmp6", "Ampd3", "Svil", "Gm5089", "Oaf", "Cdc42ep4",
	"Prdm16", "9530026P05Rik", "Grifin", "Htra3", "Prodh", "Cgnl1", "Abca8b", "Uaca", "Ets1", "D630003M21Rik",
	"Itih5", "Smad6", "Afg1l", "5830428M24Rik", "Rbpms", "Wwtr1", "Abca4", "Rassf2", "Cst3", "Eva1c",
	"Ppp1r12c", "Tns2", "Tbx18", "Gm35978", "Rapgef3", "Maml2", "Ptn", "Zic2", "1700054A03Rik", "Igsf8",
	"Notch1", "Mfge8", "Irak2", "Zcchc24", "Sash1", "Tle3", "Rcsd1", "Nox3", "Nrarp", "2610035F20Rik",
	"Klf15", "8430430B14Rik", "Emx2os", "Arhgap29", "Zic5", "Zeb1", "Gjb6", "9530091C08Rik", "F3", "Zfp36l2",
	"Plpp3", "Lbh", "Zfp648", "Amotl2", "Eps8", "4933433H22Rik", "Cnn3", "Tcf7l1", "Slc1a3", "Lrmda", "Sparc",
	"Plcl1", "Msi2", "Zbtb20", "Rbms1", "Chst11", "Gm17597", "Mmp14", "Col1a2", "Tef")

#C8 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C8 <-c("Vmn1r121", "5930403L14Rik", "Gm7861", "4930550C17Rik", "Prss45", "Cyp4f15", "Itih3",
	"Vmn1r78", "Mlc1", "Mir6236", "Ankrd66", "Prdm16", "Neurog1", "Gjb6", "Fzd2", "Lcat", "Olfr159",
	"Mir135a-2", "Scara3", "Lcn3", "Tectb", "Fgfr3", "Gm3219", "Olfr287", "Prss43", "Gm11627", 
	"2900052N01Rik", "Sdc4", "Obp2b", "Sox1", "Lyg1", "Klhdc7a", "F630040K05Rik", "Endou", "D630033O11Rik",
	"Tmie", "3110082J24Rik", "Als2cl", "Gm13872", "Apoc1", "Aldoc", "Slc7a10", "Vmn2r13", "Frmpd1os",
	"Cbs", "Gm29508", "Gm32511", "Arhgef19", "Actrt2", "Tmem89", "Dio2", "Sox1ot", "Alms1-ps2", "Mertk",
	"Gm9125", "Garem2", "Olfr675", "Gm5089", "Tppp2", "Mir6375", "Cxcl14", "Arhgef16", "Ndrg2", "Slc38a3",
	"F3", "BC016548", "Hoxb6", "Zic5", "Defb13", "Grifin", "Aqp4", "Hils1", "5033404E19Rik", "Crisp3",
	"Grin2c", "Atp1b2", "Frmpd1", "Slc12a4", "Plpp3", "Slc15a2", "Gm20757", "Pla2g7", "Gm11681", "Lfng",
	"A630077J23Rik", "Tril", "Hoxb5", "Notch1", "Fam107a", "Gpr37l1", "Olfr663", "Mir182", "Mov10l1",
	"0610039H22Rik", "Rida", "Bcan", "Apoe", "Gm35978", "Bco2", "Gm6083", "Slc7a2", "Tmem82", "Gjb2",
	"Pigs", "Calr4", "Sytl3", "Emx2os", "Cmtm1", "1700108J01Rik", "Ccdc70", "4930432M17Rik", "Ephb4",
	"Dmrta2", "Slc9a3r1", "Ranbp3l", "Olfr769", "Il18", "Eva1a", "Spatc1", "Nat8f6", "H2al1b", "Vmn2r22",
	"C030037D09Rik", "Trim34b", "Msi2", "Etnppl", "Rapgef3", "Mt2", "Slc1a2", "Entpd2", "Otx1", "Gm17660",
	"Nat8", "Ntsr2", "Gp2", "Emx2", "4930563H07Rik", "G630093K05Rik", "Olfr9", "Pax6", "Celsr1", "Gm2109",
	"LOC105245869", "Nr2e1", "Nat8f5", "Ppara", "Mir6715", "Cyp4f14", "Gpr179", "Hapln3", "Lgr6", "Prok1",
	"1700047L14Rik", "Prss46", "Ppp1r3g", "Mir7074", "Olfr854", "Sall3", "Pax8", "Mir6996", "S1pr1", "Tex24",
	"Mgat4e", "Foxn4", "1700096J18Rik", "Gli3", "Dbx2", "Oaf", "Pnp2", "Pou3f4", "Mir7684", "Nat8f2", "Sox21",
	"Cyp2d22", "Phkg1", "Serpinf2", "Slc25a34", "9230104L09Rik", "Rpe65", "Acsm5", "Smok3a", "Prr5", "Prodh",
	"Plekho2", "Fgfrl1", "1700012H19Rik", "Defa35", "Nrarp", "F930015N05Rik", "Vmn1r2", "Ghrh", "Mfge", "Gas1",
	"4933402N03Rik", "Htra1", "1700051A21Rik", "5430427M07Rik", "Hspb8", "Msi1", "B230311B06Rik")

#C9 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C9 <-c("Dcxr", "Mir7040", "Mir6973a", "C5ar2", "Adap2", "Dok1", "Tnni2", "Ltb", "H2-M3", "C1qb",
	"6330407A03Rik", "Tmem119", "Pilra", "Pbxip1", "Mir6236", "AF067061", "Samsn1", "Ccl3", "Nckap1l",
	"Crybb1", "Pon1", "Myo1f", "Ccr1", "Ctss", "Cx3cr1", "Lyn", "1810013A23Rik", "Adam26a", "Csf1r", "Fcgrt",
	"C5ar1", "Mir135a-2", "Apobec1", "Ldoc1", "Slc9a3r1", "Ccr2", "Gm12596", "Syndig1l", "Mir7086", "Csf3r",
	"Lzic", "Loxl3", "Gm11747", "Trim30a", "Lilrb4a", "Mir6983", "Bub1b", "Pon3", "Tgfbr1", "Fzd5", "Xcr1",
	"Gng10", "Nid2", "Padi2", "Fli1", "F11r", "Unc93b1", "Bin2", "A130077B15Rik", "Hfe", "Spata1", "Ccr6",
	"Acacb", "Akna", "Itpripl2", "Il10ra", "Ikzf1", "9530027J09Rik", "Mpeg1", "Siglech", "Arhgap45", "Irf5",
	"Txnip", "Zfp69", "Tectb", "Gm6607", "P2ry12", "Gm8013", "C030005K06Rik", "Nlrp3", "Cd300a", "Klhl31",
	"E230029C05Rik", "Tbc1d10a", "Gm16998", "Fcgr4", "Adora3", "Siglecf", "Rgs10", "Tlr3", "Rmi2", "Ggta1",
	"8030423J24Rik", "Nt5dc2", "Mertk", "Anxa3", "1700072O05Rik", "Ctsc", "Il6", "Hk2", "Tnfrsf14", "Gpatch2l",
	"Foxl1", "Il10rb", "Gm35978", "Rnase13", "B4galt1", "4933433H22Rik", "Tgif1", "Gm1966", "Slc1a3", "Ednrb",
	"Derl1", "4933431E20Rik", "Il21r", "Zc3h12a", "Pnp2", "Endou", "Mir6993", "Msh5", "1700017B05Rik", "Irf8",
	"Cebpzos", "4933406J10Rik", "3110082J24Rik", "Stab1", "Zfp36l1", "Vasp", "Runx1", "Gm2109", "Grifin",
	"Il13ra1", "Oaf", "Actrt2", "Gm21057", "Susd3", "B430010I23Rik", "Notch1", "Lfng", "Gpr183", "Ptgs1", "Cst3",
	"Golm1", "Prkch", "AI661453", "Adgrg1", "Timm8a2", "Dapp1", "Ifnar2", "Gm15441", "Nckap5l", "P2ry6", "Cdh26",
	"Prkag3", "4933433G08Rik", "Bcas1", "Dennd2c", "Gli3", "Car6", "Zscan4a", "Cd33", "Aga", "Mir6996", "Garem2",
	"Cryba4", "Mindy3", "Tppp2", "Clic1", "Gm35135", "Pik3ap1", "Sall1", "Tal1", "Prps1l1", "Neurl3", "Hexb",
	"Elmod3", "8430430B14Rik", "Tmigd3", "Ifitm10", "Mmp14", "Il6st", "Gm11696", "Cdkn1b", "Nfatc1", "Slco2b1",
	"Id3", "Ndrg2", "Selenop", "Hmox1", "Mtus1", "F630040K05Rik", "Bbox1", "Mt2", "Med29", "Aldh2", "Lrrk1",
	"Slc38a3", "Crlf3", "Oplah", "Bmf")

#C10 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C10 <-c("Mir470", "Vsig4", "Ccl12", "Olfr503", "Fcrla", "Olfr655", "Olfr433", "Pilra", "Ctla2a", 
	"Btnl5-ps", "Gm14548", "Olfr1085", "Slfn4", "Cd209d", "Olfr955", "BC147527", "Mir6983", "E230016K23Rik",
	"Adap2", "Gngt2", "Csf1r", "Shisa5", "Cd300c2", "Clec1b", "Btnl2", "Clec4a3", "Pira1", "Cd33", "Ccl3",
	"Lyz2", "Siglece", "C5ar2", "Ccr1", "C1qc", "Olfr894", "Mefv", "Psg23", "Ccl9", "B020014A21Rik", "Osm",
	"Gm4841", "Mir6961", "Cd86", "Cx3cr1", "Siglecf", "Mir450-2", "Mir450-1", "Tmem119", "Ccr5", "Ccl6",
	"Naip2", "Ccr6", "1700013H16Rik", "Il10ra", "Olfr116", "Arhgap45", "Olfr109", "Gm12185", "Slc7a7", 
	"Vmn1r170", "Hcar2", "Ceacam12", "Selplg", "Lilra5", "C1qb", "Irf5", "Nlrp1b", "Olfr394", "Vsir", "Cd5l",
	"Olfr1242", "Csf3r", "Siglech", "Mir142", "6030468B19Rik", "Fcgr1", "Tnfrsf1b", "Olfr110", "Milr1",
	"Mir7040", "Tnni2", "Dyrk4", "Bin2", "0610039K10Rik", "H2-Ob", "Olfr698", "Klrc1", "Cd300lf", "Ccr2",
	"Mir223", "Klri2", "Tlr6", "Abcg3", "Cd300a", "Ccl4", "Ifi205", "Olfr1163", "Lst1", "Tagap", "Mir7086",
	"Rbm47", "Trim30b", "Pik3ap1", "Vmn1r57", "Slfn3", "Klhl6", "Slfn8", "Tlr11", "Rgs10", "9130015L21Rik",
	"Adora3", "Olfr460", "Itgam", "Clec2i", "Clec4a4", "Lilra6", "Il21r", "Hexb", "Myo1f", "Serpinf1",
	"Tnfaip8l2", "Lgals9", "Tnfrsf17", "Ceacam18", "Olfr156", "Ms4a6b", "C1qa", "Pilrb1", "Havcr2", "Ctla2b",
	"Vmn1r184", "Mir542", "Bcl2a1a", "Olfr293", "Ccr1l1", "H2-Ea-ps", "Pfpl", "1110028F18Rik", "Wfdc17", 
	"Fcgr2b", "Klri1", "Slc2a5", "Gm33104", "Olfr485", "Gm5431", "Lyn", "Olfr869", "Adgre1", "Hhex", "Gdf3",
	"4931402G19Rik", "Tifab", "C5ar1", "Irf8", "Tubb1", "A530088E08Rik", "Ltb", "Prss43", "F11r", "Clec4a1",
	"Susd3", "Clec5a", "Afm", "H2-Ab1", "Olfr694", "Ikzf1", "Rfpl4", "Dennd2c", "Folr2", "Mir6973a", "Cyba",
	"Rps15a-ps4", "Olfr520", "Ildr1", "Ctss", "Gm1966", "Foxr1", "Ccr3", "Gm26641", "Fli1", "Hck", "Platr7",
	"Tnfsf14", "C3ar1", "Olfr417", "Apobec1", "Trim5", "Casp4", "Naip6", "6330407A03Rik", "P2ry13", "Gbgt1",
	"Gbp8", "Fcer1g", "Gpr34", "Gm5122", "Adam3", "Tal1", "Mpeg1", "Olfr292")

#C11 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C11 <-c("Olfr735", "Pilra", "Olfr955", "Mir881", "Zfp36l3", "Defa5", "Vmn1r122", "Olfr971", "Siglecf",
	"Ccr2", "Ccr1", "C1qc", "Ctla2b", "E230016K23Rik", "Olfr109", "Cx3cr1", "A530088E08Rik", "C5ar2", "Tmem119",
	"Cd300c2", "Cd33", "Fcgr1", "Sprr2d", "Olfr73", "Gngt2", "Ccl9", "Slamf8", "Olfr893", "Olfr110", "Ccl6",
	"H2-M2", "Gtsf2", "Siglece", "Bcl2a1a", "Ccr1l1", "Selplg", "Cysltr1", "Olfr894", "Clec4a3", "Olfr630",
	"Bin2", "Dyrk4", "Il10ra", "Fcrla", "Vsir", "Ly6a", "Tnni2", "Defa34", "Mir6961", "Mir142", "Olfr101",
	"Olfr1156", "Csf1r", "0610039K10Rik", "Olfr1053", "Ctss", "Mir6983", "BC147527", "Hhex", "Adap2", "Olfr784",
	"1700013H16Rik", "Clec4a2", "Mir223", "Ccl3", "Tnfrsf1b", "Tifab", "F9", "Ccr5", "Bcl2a1b", "Clec5a",
	"1700009J07Rik", "Gm4841", "Slc7a7", "E230025N22Rik", "Ccr6", "Abcg3", "Olfr111", "Rhoh", "Mir7086", "Cd5l",
	"Adora3", "Olfr203", "Enam", "Mmp12", "Csf3r", "Vsig4", "Itgam", "AF067061", "Siglech", "Siglecg", "Clec2i",
	"Klri1", "Trim30b", "F11r", "C1qa", "Milr1", "Clec9a", "Klrd1", "Scgb2b20", "C1qb", "Lilr4b", "Mir3471-2",
	"Clec4a1", "Mcemp1", "Mmp13", "Gpr34", "Lilra5", "1700024F13Rik", "Gm5431", "Naip2", "Lyz2", "Hexb", "Ccr3",
	"Rgs10", "Il21r", "Olfr112", "Olfr936", "Gm826", "Tnfrsf17", "Klrc3", "Myo1f", "Pik3ap1", "Olfr695", "Afm",
	"Upk1b", "Fcrls", "Irf8", "4931402G19Rik", "1600010M07Rik", "Trim12a", "Clec1b", "Gm5127", "Lalba", "Itgax",
	"Olfr469", "Gm6607", "Irf5", "Susd3", "Mir200c", "Ncf1", "Gm5728", "Olfr694", "Pira6", "Gm1966", "Clec2h",
	"Trim30a", "Olfr433", "Samsn1", "Vmn1r101", "Dusp27", "Rbm47", "Ctsh", "Gbp9", "Olfr1223", "Olfr1040",	
	"Olfr889", "Gykl1", "P2ry13", "Ccl12", "Hcar2", "Lyn", "2210414B05Rik", "Serinc3", "Reg3a", "Fli1", "Tal1",
	"E230029C05Rik", "Klrc1", "Olfr156", "Olfr107", "Slfn8", "Hpgd", "Gm10665", "Cd53", "Olfm5", "Olfr414",
	"Cdh23", "Xcr1", "Defa2", "Olfr461", "Apobec1", "Ildr1", "Gcnt1", "Havcr2", "Cd300a", "Vmn1r184", "Pirb",
	"Olfr1221", "Cd52", "Trim15", "Abi3", "Il7r", "Scgb2b19", "1700066B19Rik", "Clec2f", "Olfr432", "Trim12c",
	"Olfr494", "Ncaph")

#C12 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C12 <-c("Gm16063", "Saa2", "Vmn1r83", "Prss43", "Olfr1100", "Prp2", "Olfr1361", "Saa1", "Spaar",
	"Higd1b", "Olfr911-ps1", "Olfr888", "1700055C04Rik", "Olfr807", "AW549542", "Mir145a", "Olfr889",
	"Olfr1167", "Foxs1", "Olfr1157", "Olfr723", "Adap2", "Olfr1143", "5330439B14Rik", "Aspn", "Olfr1034",
	"Olfr485", "Abcc9", "Adora2a", "Tbx3", "Gm38416", "Cysltr2", "Eva1b", "Accs", "Gm428", "Usp17lb",
	"Defb25", "Rem1", "Olfr1311", "Tbx2", "Tas2r125", "Rgs5", "Orm3", "Olfr1242", "Olfr679", "Ecm2", "Rarres2",
	"Gprc5c", "Olfr486", "Atp2a3", "Slc6a20a", "Hrct1", "Kcnj8", "Tpbpb", "Fzd6", "Slc6a12", "Amelx",
	"Vmn1r75", "Olfr476", "Atp13a5", "Zic4", "Rbpms", "Anpep", "1700113B09Rik", "Ano1", "A4galt", "Vmn1r62",
	"Slc30a10", "Mir450-2", "Mir450-1", "Cgnl1", "Ctla2a", "Ednra", "Angptl8", "Mir542", "Rab3d", 
	"2610027K06Rik", "Rbpms2", "Slc38a11", "Hapln3", "Olfr446", "Olfr494", "Vmn1r169", "Olfr710", "Wfdc6a",
	"Ggt5", "AI463170", "Fpr-rs4", "Myl3", "Pald1", "Mylk4", "Olfr484", "Tdgf1", "Gja4", "Slc16a12", "Tmem45a",
	"Mir6384", "Cald1", "Ebf1", "Pcolce", "Pdgfrb", "Omd", "Mir6386", "Rasl12", "Gm7538", "Mir7083", "Nodal",
	"Zic1", "Bmp5", "1700023C21Rik", "Itpripl2", "Olfr1307", "Mir143", "Plscr1", "Bcas3os1", "Fgf3", "Olfr143",
	"Heyl", "Vmn1r177", "Olfr1283", "Gm2061", "Jaml", "Ptpn22", "Notch3", "Mup12", "Olfr887", "Carmn", "Tlr12",
	"Tbxa2r", "Oasl2", "Mrvi1", "E030025P04Rik", "Ephx3", "Kbtbd13", "Aadacl4", "Olfr131", "Slc6a13", "Plce1",
	"Gm16336", "Gjc1", "Mylk2", "4930515L03Rik", "Vmn2r53", "A630077J23Rik", "Nuak2", "Vmn1r24", "Olfr498",
	"Olfr347", "Cox4i2", "Vmn2r10", "Snai1", "Ifitm1", "Cdh16", "Plscr2", "5033404E19Rik", "1700092C02Rik",
	"Olfr1000", "Olfr132", "Olfr294", "1700012B07Rik", "Epha2", "Vmn2r15", "Olfr622", "Flt1", "Rnf135", "Olfr509",
	"Slc2a4", "D630033O11Rik", "Bcam", "Ntn1", "Gm15638", "Mir1954", "Dock6", "Mir7057", "Vstm4", "Uaca", "Cabp4",
	"Ifitm3", "Olfr910", "Olfr845", "Colec12", "4933407G14Rik", "Alx3", "Mir5622", "Olfr1182", "Olfr1280",
	"Sh3d21", "Pear1", "Cxcl10", "Olfr1349", "Gm15217", "2310005E17Rik", "Sntb1", "Cyb5r3", "Gm33104", 
	"E130008D07Rik", "Zdhhc19", "Vmn1r71", "Gm11827", "Pla1a")

#C13 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C13 <-c("Mir503", "Mir351", "Olfr914", "Olfr857", "Mir322", "Vmn1r58", "C430049B03Rik", "Zic4",
	"Pcolce", "Orm2", "Aspn", "Saa2", "Saa1", "Serpina1d", "2310005E17Rik", "AW549542", "Anpep", "Vmn2r50",
	"Vpreb3", "Zic1", "Tfap2b", "Gm16063", "Gbp9", "Tbx18", "Bcam", "Spaar", "Mmp1b", "Serpina1b", "Skint1",
	"Olfr1193", "A630095E13Rik", "Slc6a20a", "Slc6a13", "C730027H18Rik", "1700119H24Rik", "Nfatc4", "Hcar1",
	"Plscr2", "Higd1b", "Zic2", "Gm26641", "Olfr345", "Gm5779", "Orm3", "Omd", "F830016B08Rik", "Mir145a",
	"Sult2a2", "Adap2", "Slc6a12", "Mrgprb2", "Gpr182", "Ecm2", "2610035F20Rik", "Gm11827", "Skint4", "Il13ra2",
	"Mir483", "Vnn1", "Eva1b", "2200002D01Rik", "Mrgprb3", "Olfr77", "Gm8989", "Lamc3", "Nupr1", "Rem1", 
	"Serping1", "Rbp1", "Mir6386", "Vamp8", "Gm3458", "Gm20556", "Myzap", "Fcgrt", "1700113B09Rik", "Fam187b",
	"Aldh1a2", "Slc22a6", "Vmn2r51", "Foxd1", "Gm16287", "Vmn2r114", "Foxd2os", "Igf2os", "Foxl1", "Sult2a1",
	"Rbpms", "Foxs1", "Samd9l", "Gm15217", "Ceacam3", "Cgnl1", "Ifitm1", "A4galt", "Uaca", "Serpina1a",
	"Nlrp5-ps", "Pdgfrb", "Klrb1-ps1", "Vmn1r90", "Tmem45a", "Mir7057", "Mesp2", "Mndal", "Cfh", "Slc1a5",
	"Cdh1", "Slco1a4", "Spo11", "Fgfbp1", "Gprc5c", "Nr2f2", "Kcnj8", "Slc16a12", "Clec2d", "Tagln2", 
	"4930579J19Rik", "Ifi211", "Hsd17b2", "Gm5592", "Bmp5", "4933407G14Rik", "Cd209f", "Osr1", "St6galnac2",
	"Foxc1", "Calhm6", "Mir6396", "Gm19277", "Six2", "Foxd2", "Zic5", "Defb25", "Fgl2", "Prdm6", "Abcc9", "Fmo1",
	"Eya2", "2310001H17Rik", "Tbx3", "Hrct1", "Fn1", "Vmn2r52", "Nodal", "Col1a2", "Ephb4", "Mill2", "Jaml",
	"Bmp7", "2810433D01Rik", "Mc2r", "Cp", "Plscr1", "Adora2a", "B130024G19Rik", "Mir148a", "Lef1", "Ggt5", 
	"Zscan4d", "Ifitm2", "Mc5r", "Mylk4", "Rcn3", "Hspa12b", "Esp38", "Vstm4", "Zc3hav1", "Ctgf", "Msgn1",
	"Vmn1r15", "Ajuba", "Cysltr2", "Gm12830", "Isyna1", "Pear1", "Vmn1r66", "Gbp10", "Bmp4", "Edn3", "Igf2",
	"Gm29683", "Tas2r143", "5330439B14Rik", "Tmem30c", "Mir7050", "Iigp1", "Flt1", "Alx3", "Fzd6", "Gja4",
	"4930404A05Rik", "Mir143", "Col3a1", "Lrrc32", "Snai1", "Tpm4", "Amy1", "Slc9b2", "Colec12")

#C14 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C14 <-c("Olfr1263", "Olfr611", "Olfr612", "Olfr303", "Olfr394", "Olfr1289", "Olfr411", "Olfr971",
	"Slc22a6", "Rhox2f", "Olfr26", "Defa31", "Olfr508", "Mir6386", "Olfr646", "Defa5", "Olfr50", "Oog2",
	"Tas2r102", "Olfr74", "Olfr113", "Olfr1490", "Olfr1232", "Spo11", "Tas2r113", "Obp1b", "Olfr1145",
	"Olfr486", "Olfr298", "Olfr1173", "Olfr799", "Pcdhb17", "Vmn2r28", "Defa22", "Anpep", "Olfr745",
	"Olfr1461", "Vmn2r106", "Olfr1182", "Olfr59", "Slc6a13", "Gm4884", "Olfr473", "Aldh1a2", "Olfr1462",
	"Defa17", "Olfr308", "Igf2os", "Olfr1148", "Olfr678", "Slco1a4", "Mir7057", "H2-Q1", "2310005E17Rik",
	"Vmn1r13", "Gm10487", "Olfr397", "Olfr1446", "Scgb2b15", "AA792892", "Lce1a1", "Olfr403", "Olfr472",
	"Vmn1r170", "Gm3458", "Rbp1", "Vmn2r15", "Olfr1451", "Zic4", "Olfr1085", "Mir483", "Vmn2r105", 
	"C730027H18Rik", "Vmn1r89", "Olfr1140", "Krtap5-5", "Olfr1312", "Gm16063", "Vmn1r171", "Olfr137",
	"Olfr1472", "Olfr1502", "Olfr733", "Clec4n", "Olfr1036", "Olfr138", "Defb2", "Gm14295", "Olfr715",
	"Olfr1454", "Vmn2r33", "Slc6a20a", "Vmn1r231", "Olfr259", "Tas2r107", "Gbp4", "Olfr1453", "Zfp616",
	"Gm14305", "Pcolce", "C87977", "4930594M22Rik", "Olfr1162", "BB287469", "Plscr2", "Slc25a2", "Olfr1157",
	"Olfr1495", "Olfr1494", "Olfr1153", "Omd", "Olfr913", "Gm13119", "Olfr910", "Olfr124", "Vmn1r62",
	"Olfr868", "Olfr139", "A4galt", "Gm14325", "Slfn5os", "Slc22a29", "Gm13088", "A630077J23Rik", "Vmn2r16",
	"Olfr889", "Edn3", "Uaca", "Olfr167", "Olfr746", "Cd3d", "Trim5", "1700018B08Rik", "Gpr182", "Mrgprb1",
	"Samt4", "Vmn2r80", "Olfr747", "Foxd1", "Cyp2c67", "Olfr955", "BC021767", "A430089I19Rik", "Olfr1501",
	"Olfr133", "Mir126b", "Ocln", "Lama1", "Igf2", "Olfr1317", "Olfr814", "Nfatc4", "Olfr293", "Vmn2r110",
	"Serpina1d", "Zic5", "Ctla2a", "Eya2", "Six2", "Olfr837", "Olfr1053", "Olfr1208", "Olfr358", "Pramef20",
	"Olfr476", "Olfr395", "Slc22a2", "Olfr884", "Olfr384", "Olfr952", "Cyp2c53-ps", "Olfr159", "Cldn5",
	"Olfr819", "Gm15217", "Olfr166", "Slc22a14", "Olfr165", "Cyp2c39", "Serpinb9g", "Iigp1", "Olfr474",
	"F830016B08Rik", "Zic1", "4930579J19Rik", "Sytl3", "Anxa2", "Ccl7", "Slfn5", "Olfr43", "Gm12886",
	"Gbp10", "Tmem45a", "Olfr1170", "Slc22a27", "Olfr1247", "Bcam", "Cgnl1", "1700055C04Rik", "Olfr854")

#C15 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C15 <-c("Gm14548", "H2-M2", "Olfr1028", "Pira1", "Olfr1307", "Pilra", "Mir223", "Clec1b",
	"Bcl2a1a", "Olfr693", "A530088E08Rik", "Hhex", "Olfr432", "Ccr1l1", "Olfr414", "Cd300c2", "C1qa",
	"BC064078", "Siglec", "Gbp9", "Dyrk4", "Clec9a", "1700009J07Rik", "Gm21119", "Olfr111", "Ccl12",
	"Enam", "Siglech", "Ccr5", "Adap2", "Pira2", "Gdf3", "Clec4a3", "Clec5a", "Olfr109", "Fcrla",
	"Clec4a2", "C5ar2", "Hpgd", "Gngt2", "C1qc", "Cd5l", "Olfr112", "Cd52", "Csf1r", "Klri2", "Clec2i",
	"Cd86", "Clec4a1", "Mir6983", "Cd33", "Il10ra", "Bin2", "Ncf1", "Mmp12", "Tmem119", "Cd300a",
	"0610039K10Rik", "Adora3", "Dennd2c", "Ccr6", "Slc7a7", "Tnfrsf1b", "Tifab", "Mir6961", "Khdc3",
	"1700024F13Rik", "Cx3cr1", "Vsir", "Lalba", "Jchain", "Lilrb4a", "Tlr13", "Tnni2", "Gm826", "Trim30a",
	"Rgs10", "C1qb", "BC147527", "4930563M21Rik", "Gm1966", "Ly9", "Adam3", "Tnfrsf17", "Lilr4b", 
	"Gpr160", "Sash3", "Fcgr2b", "Arhgap45", "Gm4841", "Gm156", "Mrln", "Cks1brt", "Mir7086", "Gbp8",
	"Mir200c", "Ltb", "Ccr1", "Fli1", "Serpinb12", "Ccl6", "Ccr2", "Samsn1", "E230029C05Rik",
	"Olfr110", "1700100I10Rik", "Mir6973a", "AF067061", "Gcnt1", "Fcrls", "Il21r", "Olfr108", "Mir142",
	"Adgre1", "Olfr107", "E230016K23Rik", "Ucp2", "6330407A03Rik", "Milr1", "Tmem273", "Lta", "Rbm47",
	"Ccl4", "Pik3ap1", "Dapp1", "Ctsc", "Selplg", "B2m", "Cpb1", "Hpgds", "Ccl3", "Skint4", "Ncf4",
	"Mir141", "Ikzf1", "4930512M02Rik", "Abi3", "Cd53", "Gpr35", "Skint3", "C5ar1", "E230025N22Rik",
	"Itgam", "Siglecf", "Ctss", "Gpa33", "Klri1", "Pon3", "Gm33619", "Tmsb4x", "Stab1", "Fam71a",
	"Susd3", "Fgl2", "Fcgr4", "Uts2b", "Inka1", "Hexb", "Abcg3", "Gpr183", "Prkcd", "F11r", "Klk12",
	"Sqor", "1700021F07Rik", "Lrmp", "Cmtm7", "Irf5", "Batf3", "Ctsh", "C3ar1", "5430434I15Rik",
	"Slco2b1", "Cass4", "Lyn", "2310069B03Rik", "Gm21057", "1700122C19Rik", "Ctsz", "4930552P12Rik",
	"Afm", "Cebpzos", "Apobec2", "Slamf8", "Gm10578", "Mrgpra3", "Tlr3", "Clec4g", "Rgs1", "Msh5",
	"Dusp27", "Lncbate1", "A130077B15Rik", "Cdrt4os1", "Siglecg", "Fam162b", "F9", "Hps3", "Bcl2a1b",
	"4933406I18Rik")

#C16 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C16 <-c("Krtap4-8", "Gm22650", "Snora35", "Cyp4a14", "Krtap4-9", "Cyp4a30b", "Olfr322", "AY702102",
	"Mir764", "Olfr798", "Krt82", "Krt90", "Tshz2", "Ces1d", "A630075F10Rik", "Mir669g", "Mir7072", "Card14",
	"Cd3d", "Gm12169", "Pom121l12", "2610528A11Rik", "Adgrg3", "Bmp10", "Kank4", "Ifi208", "Myzap", "Tbata",
	"Krt87", "Mir466g", "Lypd1", "1700017J07Rik", "4930487D11Rik", "Olfr71", "4933406M09Rik", "Mir669m-2",
	"Trim58", "Hist1h1e", "Cyp2e1", "Drd5", "4930578G10Rik", "Mir466n", "Mir669o", "Kank4os", "Galp",
	"1700121N20Rik", "Mir6385", "Krtap13", "Nxph3", "Ccl21a", "Gm8630", "4930500F04Rik", "Gm10768", "Slc17a8",
	"Il1f6", "Gm12794", "Oas1f", "Ceacam16", "Glyatl3", "Krt25", "Abca12", "Samd7", "Ces1c", "1700012A03Rik",
	"Trhr", "Vsig1", "Bcas3os2", "Cyp4a10", "Lrriq4", "Zp1", "Aox4", "Impg2", "Abcc12", "Aplnr", "Ptch2",
	"Nlrc5", "Htr4", "Fbxw24", "Mir1298", "Asb18", "Nr1h4", "Ankrd22", "Htr2c", "Mir598", "Krt26", "Gm2696",
	"Fbxw27", "Mki67", "Mir6969", "Snord11", "Mir375", "Trpc3", "Slc15a5", "Mettl7b", "2210407C18Rik",
	"Mir6923", "Slc7a3", "Abca15", "Slc47a2", "Spta1", "Fgf21", "Cyp4a12b", "Gm30505", "4930579J19Rik",
	"Il1bos", "Hsh2d", "Mir448", "Tmc3", "Gm38416", "Cd244a", "Aqp1", "Adamdec1", "Elfn1", "Gm21221", "Gja10",
	"Mirt1", "4921529L05Rik", "Antxrl", "Cpne4", "Gm5460", "4930455G09Rik", "Krt24", "Olfr530", "Mir7677",
	"Fezf2", "D030068K23Rik", "Pla2g2d", "Adad1", "Syce1", "Gm19668", "1700113H08Rik", "Ctgf", "Ly6k", "Muc5b",
	"Mir6365", "Mst1r", "A630019I02Rik", "Actl11", "Vwc2l", "Cacng5", "Khdc1c", "Rp1l1", "9330162B11Rik",
	"Trem3", "Acaa1b", "Spon1", "4933430M04Rik", "Dmrtc1a", "8430436N08Rik", "Psmg3", "Olfr967", "4931429P17Rik",
	"Drd1", "Fcrl5", "Olfr434", "Tmco5", "Hs3st4", "Ctf2", "Gm12354", "Ifi207", "Nrsn2", "Kif26a", "Mir466l",
	"Sntb1", "Chrna5", "Il1b", "Serpinb6c", "Sla2", "Ffar2", "BC049352", "Drc7", "Layn", "Tmem159", "Nlrp6",
	"Nek11", "Fhl5", "Gm19402", "Gas2l3", "D630023F18Rik", "Bcl11b", "Tex13b", "Trim50", "Gimap1", "Mfsd2b",
	"Slc5a9", "Mir6406", "Podnl1", "Arg2", "Ces1e", "Gm6042", "Npy4r", "Mir200b", "Abca14", "A930041C12Rik",
	"Grik3", "Notum", "Pcsk5", "Grp", "Ptger3", "Mir6974")

#C17 markers: FDR <= 0.05 & Log2FC >=1.15, all 15 sorted by Log2FC
known_markers_C17 <-c("Nkg7", "Olfr328", "Ppbp", "Eppin", "Cntnap3", "Car1", "Omt2a", "Hba-a2", "Gad2", "Dlx6os2",
	"Fam240b", "Ankk1", "2210407C18Rik", "Mir7072", "Crhbp")

#C18 markers: FDR <= 0.01 & Log2FC >=1.25, all 127 sorted by Log2FC
known_markers_C18 <-c("Wfdc15a", "Gm11938", "Tmco2", "Svs6", "Vmn1r233", "Wfdc15b", "Fam183b", "Krtap4-1",
	"Olfr1392", "Olfr808", "Olfr1340", "Agxt2", "Olfr1391", "Tmem215", "Gm11937", "Olfr316", "Krtap3-2",
	"Nrl", "Krtap1-4", "Rhox4a", "Mir1952", "Cdsn", "Serpina1f", "Olfr1427", "4933400B14Rik", "Exph5",
	"Dact2", "4933428C19Rik", "Olfr809", "Rtn4rl2", "Ubxn10", "Ms4a6c", "Gm8369", "Serpina1e", "Gm4719",
	"Lrit2", "Nptxr", "Olfr464", "Slc30a3", "Gm13710", "4930452N14Rik", "Krtap4-2", "Olfr314", "Trim54",
	"Olfr1428", "Spag11b", "Spata31d1a", "Mir8106", "Svs2", "Olfr10", "Calb1", "Olfr1441", "Sntg2", 
	"Ptgs2", "Rspo1", "Olfr789", "Olfr16", "Gm8267", "Defb34", "Ocm", "Lce1m", "Eef1e1", "Olfr877", "Crhr1",
	"Themis", "Vwa5b1", "Krt12", "Slfn10-ps", "Rtn4r", "4921517D22Rik", "Vsig2", "Pla2g2c", "Krtap1-5",
	"Olfr145", "Gm11559", "Clca3a2", "Olfr56", "Mir216c", "Dgat2l6", "Olfr317", "Cerkl", "Gm11568", "Mir6239",
	"Cd163", "A430090L17Rik", "Arsj", "Prss54", "Sh2d7", "Unc5d", "Grm2", "Mir5107", "Kcnv1", "Itpka", 
	"1700009C05Rik", "S100a7a", "Sh3rf2", "Gm6654", "Psors1c2", "Mpl", "4930532M18Rik", "Olfr1383", "Samd3",
	"Mir802", "Mir215", "1700006F04Rik", "Glt8d2", "Ovol2", "Tmem213", "Slc24a4", "Uchl4", "Neurod6", "Smoc2",
	"Car4", "Trpc5os", "Retnlb", "Lingo1", "Krtap4-7", "Snai3", "Sstr3", "Olfr120", "Krtap3-1", "Cux2",
	"Ackr2", "Apol7b", "Panx3", "Kcnh4", "Gzmn")

#C19 markers: FDR <= 0.01 & Log2FC >=1.25, all 126 sorted by Log2FC
known_markers_C19 <-c("Tshz2", "Hist1h1e", "Svs5", "Olfr685", "Gm12354", "Mir448", "Krt79", "A630075F10Rik",
	"Gm9112", "Gm5166", "Mir1264", "Htr2c", "Lipo1", "Trem3", "Mir6354", "Khdc1c", "Ankdd1b", "Gm13580",
	"Olfr881", "Lipo2", "Mup21", "Aox4", "Gm8267", "Olfr890", "Adcy10", "Smim6", "Gm32141", "Ear2",
	"AY702102", "Neil2", "Teddm1a", "Lipo4", "Itprid1", "Olfr235", "Mir8113", "Platr21", "BC061237",
	"Slc12a1", "Cldn34-ps", "Slpi", "4921517D22Rik", "Ankrd22", "2410012E07Rik", "Olfr519", "Ip6k3",
	"Pbk", "St8sia6", "Pglyrp4", "Htr7", "Fam170b", "Cst6", "4930556N09Rik", "Scnn1b", "1700017J07Rik",
	"Pip", "Krt73", "Mir1952", "Tnxb", "Rbpjl", "Mir1298", "Htr4", "Slc17a6", "Aard", "Igdcc3", "Lamp3",
	"Ovch2", "Gm2848", "Gm4719", "Cyp8b1", "Cd27", "Grm2", "Vwa5a", "Evc2", "Olfr930", "Pigr", "Gal3st2",
	"Ush1c", "Ucma", "Has3", "Olfr1396", "Wfdc5", "Mir6356", "Wdr93", "Olfr874", "S100a3", "Slc12a7",
	"Platr14", "Neurod6", "Abca12", "Panx2", "Evc", "Ang6", "Fthl17a", "Lipo3", "Treml2", "Mir6919",
	"Gpr161", "Gm40893", "1700028P14Rik", "Mir6923", "Iqcf6", "Scgb2b1", "Rhov", "Olfr878", "Nell1os",
	"Naa11", "4933400B14Rik", "Kcnj11", "4930550L24Rik", "Lrrc15", "Gm10466", "Tnp1", "Tll1", "Mir3970",
	"Ces1f", "Tapbpl", "4930563F08Rik", "Gm6040", "E130112N10Rik", "Chrna6", "4930546C10Rik", "Chrnb3",
	"B430306N03Rik", "Clec7a", "Kcns1", "Ch25h")

#C20 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C20 <-c("Skint1", "Olfr684", "Tpsb2", "Lce3c", "Olfr683", "Serpina3b", "Lce3e", "4933412E24Rik",
	"Olfr665", "Olfr168", "Tpsab1", "Olfr907", "Skint10", "Olfr1200", "Mir547", "4933408N05Rik",
	"Vmn1r84", "Mup17", "Krtap19-4", "Olfr99", "Olfr1189", "Olfr887", "Lce3d", "Vmn1r11", "4933406K04Rik",
	"Olfr740", "Sprr2k", "Olfr25", "Sprr2j-ps", "Gm12171", "Olfr906", "Lce3f", "Olfr785", "Mup3", "Vmn1r200",
	"Mup7", "Tpsg1", "Dppa1", "Olfr1437", "1700121N20Rik", "Olfr885", "Olfr68", "Skint7", "Olfr1246", 
	"Olfr1204", "Rhox4a", "Olfr786", "Prop1", "Olfr599", "Olfr1180", "Gm3434", "Oosp1", "Mir7007", "Olfr1512",
	"Olfr685", "Olfr1231", "Tas2r109", "Olfr170", "Olfr774", "Olfr169", "Gm5916", "Prpmp5", "Olfr1448",
	"Prss29", "Olfr1202", "Olfr1181", "Tarm1", "Rprl1", "Vat1l", "Mir201", "Olfr978", "Kcng1", "Btg1-ps2",
	"Vmn1r-ps103", "Olfr1436", "1700063O14Rik", "Vmn1r44", "Krt73", "Adgrg3", "Mir467h", "Tcerg1l", 
	"1700029M20Rik", "Skint8", "Msgn1", "Olfr1208", "Olfr1467", "Olfr376", "Mir290a", "Olfr1205", "Acrv1",
	"Zp3r", "Gm5414", "Layn", "Oas3", "Klrc3", "Mc4r", "Ms4a5", "Slc43a3", "Gm11985", "Fam19a1", "Olfr64",
	"1700123L14Rik", "Npr3", "Klra5", "Olfr773", "Lincmd1", "Olfr63", "Olfr970", "Vmn1r197", "4930422M22Rik",
	"Teddm3", "Mup21", "Stac", "Cldn23", "Gm19784", "Olfr1469", "4933407E24Rik", "Ptpn22", "Adam21", "Cks1brt",
	"Olfr527", "Mup11", "Cst8", "Olfr1450", "Mir3473e", "Nat3", "Ebf4", "Olfr1387", "Olfr1226", "Mbd3l2",
	"Cpa3", "Lipn", "Gm5144", "Platr20", "Olfr945", "4930448C13Rik", "Crym", "Spdye4a", "Mup19", "4933417A18Rik",
	"Depp1", "4933425L06Rik", "Olfr814", "C7", "Ctcflos", "Dnah9", "Olfr689", "Olfr273", "A1bg", "Cldn34d",
	"Mir7069", "Bcl11b", "Raver2", "Strc", "Ppp1r14d", "Il1f6", "Npvf", "1700034K08Rik", "Olfr516", "Ms4a3",
	"Adam5", "Olfr1434", "Olfr775", "Olfr796", "Klk1b3", "Gm5615", "Pzp", "Lce3b", "Olfr371", "Slc5a9",
	"L3mbtl4", "3300002P13Rik", "Kcnj15", "Vmn1r43", "Cpb1", "4921539E11Rik", "Rab3c", "Fam124b", "Olfr1447",
	"Gm3428", "Omt2b", "Gapt", "Olfr1432", "Olfr1329", "Erg", "Mir6918", "Olfr1031", "Fras1", "Lrrc19",
	"C530044C16Rik", "Serpina3c", "Gm5512", "Ankrd34b", "Spata45", "Olfr60", "4930407I19Rik", "Toporsl",
	"Tmem163", "Olfr976", "Olfr1395")

#C21 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C21 <-c("C920009B18Rik", "Klk1b3", "Klk1b4", "Hs3st4", "Cabp5", "Olfr1497", "Gm8630", "Raet1d",
	"Olfr319", "Raet1c", "Gm10825", "Raet1a", "Olfr1499", "Sval2", "Prss1", "Bsph1", "Cyp2c70", "Olfr577",
	"D730048I06Rik", "Clec4g", "Gm5771", "Cd209a", "Fam90a1a", "Anks4b", "Chrna5", "Defb46", "9230110F15Rik",
	"Olfr1487", "Defb4", "Ccr3", "Foxp2", "Syt6", "Mettl7b", "Tas2r126", "Gm19668", "Gm17727", "Aldh3a1",
	"Capns2", "Mcpt-ps1", "Ces2d-ps", "Olfr178", "Il18r1", "Samd7", "Krt17", "A630019I02Rik", "Gm4850",
	"1600027J07Rik", "Gm3238", "Krt82", "Defb6", "Tbpl2", "Sash3", "Pigr", "Pih1d3", "Krt16", "Tex19.2",
	"Krt90", "Olfr458", "Gm3428", "Olfr177", "Tas2r135", "Mir6923", "Pcp4", "Tchhl1", "1810020O05Rik",
	"Ccr1l1", "Mir1941", "A630023P12Rik", "Klk1b5", "Raet1e", "Olfr641", "Islr", "Myh3", "Olfr320",
	"3100003L05Rik", "Slc1a7", "Dpp4", "Klk1b22", "Gm12863", "Ceacam9", "Gm2696", "Aldh1b1", "Emc3", 
	"Arhgap25", "Ceacam15", "Gal3st2c", "Gal3st2", "Lrrc15", "Fancd2", "C7", "BC053393", "Olfr974",
	"Btbd35f7", "Olfr1413", "Gm10389", "Krt24", "H60b", "Mup17", "Zfpm2", "LOC102631757", "2210010C04Rik",
	"Mndal", "1700092K14Rik", "Slc37a1", "Gm19276", "Xkr7", "Gal3st2b", "Adam21", "Gm2848", "Gm5859",
	"Gm10334", "Gm26760", "Rai14", "Tmem267", "Cyp4a30b", "Ly6d", "Thsd7b", "1700028J19Rik", "Olfr535",
	"B020004J07Rik", "Olfr576", "Ces1c", "Olfr541", "Prss2", "Krt86", "Cyp4a14", "Treml1", "Hnrnpa1l2-ps2",
	"Pate4", "Podn", "Mir21c", "Sval1", "Saxo1", "Nlrp6", "Olfr59", "Adam6b", "Olfr533", "4930556M19Rik",
	"Defb50", "Gm19402", "Pla2g4c", "Grik3", "A130030D18Rik", "Pnliprp1", "Serinc2", "Igsf5", "Olfr1107",
	"Aqp3", "D930048N14Rik", "Tbata", "Fezf2", "Gfral", "Mir592", "Vsig1", "Chrna4", "Slc16a10",
	"4930423M02Rik", "1700011B04Rik", "Mir7069", "1700007P06Rik", "Ptk7", "Il1rl1", "Gm17644", "Diras2",
	"Grb7", "Slc34a2", "Slc35f3", "Cyr61", "Espnl", "Mir203", "Hcrtr2", "Bcl11b", "Mdh1b", "Olfr47", "F10",
	"4933430N04Rik", "Pnliprp2", "Spata46", "Il1f5", "Map3k7cl", "Daw1", "N4bp3", "Gm9899", "1700013D24Rik",
	"Igfbp4", "Defb14", "Fam163b", "Serpina3b", "Galnt9", "4930553E22Rik", "Krtap16-1", "Defb5", "Ces2f",
	"1700084F23Rik", "Olfr194", "Col5a1", "Gm12295", "Crym", "Gm7694", "Krt222")

#C22 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C22 <-c("Lhx6", "Mir486", "Crhbp", "Mir3107", "Fam240b", "Slc32a1", "Krt31", "Krt34", "Olfr286",
	"Gad2", "Olfr287", "Gm14204", "Dlx1as", "Gm15816", "1700120E14Rik", "Pvalb", "Krtap15", "Gad1", "Gm5833",
	"Rab3b", "Cort", "Dlx1", "Vmn1r217", "Dlx6os2", "4930477N07Rik", "Btbd11", "Pcare", "Lrrc38", "Cacna2d2",
	"H1fnt", "Nxph1", "Olfr33", "Gm12408", "Mir3967", "Nipsnap3a", "Minar1", "1700025F24Rik", "Sox6",
	"Pcp4l1", "Gjd2", "Kcnmb2", "Krtap14", "Zfp641", "1700031M16Rik", "Olfr270", "Grhl3", "Mir1a-1", 
	"Mir30c-1", "Dlx6os1", "Mir30f", "Gad1os", "Rbp4", "Rpp25", "1700010J16Rik", "Akain1", "Dlx6", "Dlx5",
	"Rerg", "Tmem132c", "Kcns3", "Ffar4", "4933432K03Rik", "Oprd1", "Mir7234", "Gm4371", "Mir429", "Paqr5",
	"Slc5a4a", "Pde6c", "C130080G10Rik", "Cntnap3", "Ifna7", "Erbb4", "S100z", "Cib4", "Ch25h", "Rnaset2b",
	"Gm12339", "Gm13648", "Kcnip1", "Spp2", "Ngf", "Krt32", "Ahsg", "Filip1", "Krt35", "Adamts15", "Nkx6-3",
	"Sytl5", "A530058N18Rik", "Mir1a-1hg", "Slfn4", "Olfr355", "Mir673", "Hist1h2ba", "Th", "Olfr281",
	"Olfr282", "Rbms3", "Tmem233", "Krt13", "Cox6a2", "Elfn1", "Olfr30", "Tmem132cos", "Ces2g", "Grip1",
	"Afap1", "Acan", "Hist1h2aa", "Gm12159", "Fetub", "Gpr176", "Olfr332", "Slc27a2", "Srpx2", "Krt36",
	"8030443G20Rik", "4930426D05Rik", "Ntn4", "Ank1", "Myo5b", "Rln1", "Hhla1", "Vwc2", "4930590A17Rik",
	"1700054M17Rik", "Cntnap4", "Trim67", "4930573O16Rik", "Mrpl43", "Gm20611", "Iltifb", "Glrp1", "Mir7212",
	"Grip1os1", "1700016K05Rik", "Ces2f", "Atp5g1", "Klhl14", "BC018473", "Mir3473e", "Gstm2", "Rxfp3",
	"Lrrc52", "Ptprm", "Grip1os2", "Them5", "Nppb", "Mir466j", "Lyz1", "Gm15881", "Inpp5j", "Vmn1r85",
	"6530411M01Rik", "Kdelr1", "Gm5941", "Acoxl", "Als2cr11", "Ccdc102a", "9230112J17Rik", "Actl7b",
	"4930552N02Rik", "Col19a1", "Bend4", "Grip2", "Mir7027", "Dynap", "Sp9", "Gm8580", "Gm14137", "Cr2",
	"Sst", "Reln", "Grik1", "Vat1", "Npy", "Pwp1", "Ascl2", "Gm2721", "Moxd1", "Mafb", "Ifit3", "Gsdme",
	"Spsb4", "Kcnc2", "Gm572", "Tox3", "Gm4871", "Lrrd1", "Olfr559", "Mir300", "Adad2", "Defb21", "Mir376a",
	"Izumo1r", "Morn5", "Mir3099", "Guca2b", "Olfr320")

#C23 markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_C23 <-c("Olfr285", "Adarb2", "Htr3a", "Vip", "Olfr284", "Glrp1", "Slc32a1", "Krtap9-1", "Gad1",
	"Gm14204", "Gad2", "Dlx1", "Htr3b", "Olfr1330", "Dlx1as", "Dlx6os1", "Prox1", "Mir3967", "Mir7080",
	"Krt71", "Asb17os", "Kit", "Dlx6os2", "Gm19589", "Gad1os", "2310034C09Rik", "Krt39", "Pnoc", "Krt6a",
	"Ngf", "Igf1", "Taar7f", "Mir30f", "Ppbp", "Krtap13-1", "6530411M01Rik", "Tmem225", "Dlx6", "Metap1d",
	"Mir30c-1", "Tnfaip8l3", "Themis3", "Btbd11", "Pldi", "Gjd2", "Cpne7", "Krt4", "Gm1140", "4933432K03Rik",
	"Slc44a5", "Erbb4", "Tmem132e", "Msantd1", "Elfn1", "D930016D06Rik", "Mir6408", "Chat", "Gm3704", "Npas1",
	"Mir135a-1", "Bpifa5", "Reln", "Bpifa1", "7420701I03Rik", "Mir701", "1700003G13Rik", "Chrna7", "Krtap10-4",
	"Gm10787", "Bpifb2", "Phf11d", "Kcnip1", "Olfr1383", "Gm29683", "4930444F02Rik", "Kcnmb2", "Gm11985",
	"Gja10", "Rgs12", "Olfr555", "Olfr986", "Uts2r", "Clrn3", "Dlx5", "Gdf2", "Slfn3", "Npy", "Tcerg1l",
	"Zfp536", "Gm10845", "Spp2", "Nnmt", "Akp3", "Mir7-2", "Drd2", "Gm10782", "Hgfac", "Grpr", "9430019J16Rik",
	"Egfr", "Rbms3", "Chrna5", "Sp8", "Rgs16", "Ffar4", "Krt5", "Gm11567", "Mir7025", "A530058N18Rik",
	"LOC102632463", "B430212C06Rik", "Slfn4", "Fxyd6", "Grip1", "Arhgap40", "Sult5a1", "Spert", "Cers3",
	"4930448I06Rik", "Grip1os1", "Gm6623", "Erich6b", "Prok2", "C130080G10Rik", "Th", "Nr2e1", "Ankk1", "Necab2",
	"Rspo4", "C8b", "Tmem235", "Olfr1335", "Cacng5", "4930463O16Rik", "Dpysl5", "Krt13", "Myo3b", "Gapt",
	"Slc22a22", "1700025F24Rik", "1700121N20Rik", "Il22ra1", "Pde11a", "Gm20741", "Npas3", "Mir149", "Nrip3",
	"Calb2", "Tex19.1", "Krt19", "Tac2", "Col26a1", "Tmprss11e", "Mir328", "1700054M17Rik", "9630028B13Rik",
	"Krt79", "C8a", "Vmn1r46", "Cldn1", "St6galnac1", "Glyctk", "Mir3099", "2310039L15Rik", "Pf4", "Pcp4l1",
	"Gm10318", "Gm16364", "Cpa5", "Olfr1384", "Olfr638", "Rdh1", "Mir551b", "4930447N08Rik", "Sema3c", "Rpp25",
	"Cort", "Gm4371", "Rab3c", "Pbk", "Cybrd1", "Gdf10", "Cyp2d40", "Glb1l2", "Eppin", "Krtap15", "2310057N15Rik",
	"Krt15", "4930402F11Rik", "Wfdc6a", "Cacng4", "C530008M17Rik", "Gm3285", "Siah3", "Asic4", "Npy2r", "Nr2f2",
	"Gm7607", "Fut4-ps1", "Adra1a")

#C24 markers: FDR <= 0.01 & Log2FC >=1.25, all 87 sorted by Log2FC
known_markers_C24 <-c("Olfr476", "Vmn2r13", "Vnn1", "Vmn1r78", "Defb10", "Vmn1r65", "Usp17le", "Scgb1b24",
	"Sytl3", "Vmn1r79", "Mir6236", "Olfr982", "Olfr663", "Bsph2", "Scgb2b24", "BC049730", "Cd177",
	"Psg26", "Vmn1r64", "Mir3079", "C87414", "Olfr1320", "Olfr9", "Ceacam5", "5930403L14Rik", "Mlc1",
	"Vmn2r43", "AY761185", "Pax4", "4930550C17Rik", "Prdm16", "Prss45", "Gfap", "Defb12", "Neurog1",
	"Prok1", "Hils1", "Gjb6", "Itih3", "Ankrd66", "Cyp2d10", "Endou", "Arx", "Tmie", "Als2cl", "Mir135a-2",
	"Tppp2", "Ntsr2", "Vmn2r22", "Tmem176a", "Scara3", "Vmn1r80", "2900052N01Rik", "Gm9125", "Tmem89",
	"Gm2109", "Vmn2r101", "Smpd5", "4930539C22Rik", "Mir101c", "Mov10l1", "F630040K05Rik", "Upk3bl", "Cbs",
	"Olfr286", "Ndrg2", "Vmn1r75", "Crisp3", "C030037D09Rik", "Aqp4", "1700012H19Rik", "Slc7a2", "Ranbp3l",
	"Arhgef16", "Psg-ps1", "Slc9a3r1", "D630033O11Rik", "AA619741", "Fgfr3", "Crisp1", "5033404E19Rik",
	"Lcat", "Tectb", "Hoxb7", "Lcn3", "Gm13872", "0610039H22Rik")

#C25 markers: FDR <= 0.05 & Log2FC >=1.15, all 99 sorted by Log2FC
known_markers_C25 <-c("Mir450-2", "Mir450-1", "Olfr1336", "Mir542", "Etv2", "Gpr17", "Slc22a27", "Has2os",
	"Mif4gd", "Dcpp1", "Doxl2", "Pnlip", "Ms4a4d", "Pdgfra", "Olfr1507", "Map7d3", "Dcpp2",
	"Ppp1r2-ps3", "Mir503", "C1ql1", "Mir351", "Plac1", "Lims2", "Has2", "Trim25", "Nipal1", "Vwa2",
	"C430049B03Rik", "Oas1h", "Fam89a", "Wincr1", "Tmem176a", "Mir7025", "Gm4632", "Tcl1", "Mir7678",
	"Sapcd1", "Susd5", "A930003O13Rik", "Tmem176b", "Antxrl", "Pramel1", "Gm5091", "Neu4", "Caskin2",
	"Gsx1", "Cacng4", "B3gnt5", "4922502D21Rik", "Tac1", "Inava", "Mir6406", "AA545190", "Sun5", "Hao1",
	"Spg7", "Irx6", "Mir5621", "Slco1a5", "Zfp488", "5031434C07Rik", "4930486F22Rik", "Arhgef16", "Myt1",
	"Rab38", "Lad1", "Lhfpl3", "Fam240a", "Chst5", "Pla2g2f", "Sox10", "Slc5a8", "Plpp4", "D730050B12Rik",
	"Xkr5", "Cd24a", "AY702103", "Ntn1", "Bpifc", "Rpl13", "Slc38a11", "Plut", "Krtap21-1", "Sstr1",
	"Mnd1-ps", "Sapcd2", "Cdo1", "Gm15908", "Gzmb", "Prkg2", "A930009A15Rik", "Ttr", "Matn4", "Kif18b",
	"Cspg4", "Snord68", "Ltbr", "Calb2", "Mtap")



######################################

#Subset the markerGS object to just these known markers
se_idx <- which(rowData(markerGS)$name %in% known_markers_C4)  

#Get an index of these from the summarized experiment
subset_markerGS <- markerGS[se_idx,]  

##Plot it out:

pdf(file="C4_TxComp_Heatmap_2024-03-21.pdf", width=15, height=30)
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
  labelMarkers = known_markers_C4,  #Label the markers of interest we're plotting here
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
  labelMarkers = known_markers_C4,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  # Change this to TRUE
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

# Print z-scores from plotMarkerHeatmap
write.csv(z_scores, "C4_byTx_zscores_2024-03-21.csv", row.names = TRUE)

##########################################
##########################################
##########################################
##########################################

