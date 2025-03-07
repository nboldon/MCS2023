
######################################################
######################################################


Setup an interactive session
salloc --account=eon -t 0-08:00:00 --mem=256G --nodes=5 --ntasks-per-node=16

#Updated conda env 12-2023
module spider load miniconda3/23.1.0
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
#setwd("/project/eon/nboldon/MCS2023/ClusterID")
setwd("/Volumes/DataBox/MCS2023/Tx_Comp")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
#projMCS7 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS7)
table(projMCS7$Clusters)


###############################################################
###############################################################
###############################################################


# Define cell types
cellTypes <- list(
  glut = projMCS7[projMCS7$Clusters %in% c("C18", "C19", "C21"), ],
  glutPrecursor = projMCS7[projMCS7$Clusters %in% c("C15", "C16", "C17", "C20", "C25"), ],
  gaba = projMCS7[projMCS7$Clusters %in% c("C22", "C23"), ],
  microglia = projMCS7[projMCS7$Clusters %in% c("C10", "C11"), ],
  oligo = projMCS7[projMCS7$Clusters %in% c("C2", "C3"), ],
  oligoPrecursor = projMCS7[projMCS7$Clusters %in% c("C5", "C6"), ],
  glutOligo = projMCS7[projMCS7$Clusters %in% "C4", ],
  endo = projMCS7[projMCS7$Clusters %in% c("C12", "C13", "C14"), ],
  astro = projMCS7[projMCS7$Clusters %in% "C8", ],
  astroPrecursor = projMCS7[projMCS7$Clusters %in% c("C1", "C7", "C9"), ],
  glutAstro = projMCS7[projMCS7$Clusters %in% "C24", ]
)



############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################


################################################################
########## Get top genes for each cell type  #####################
################################################################


# Create a named vector with NA values for each cell
cellTypeVector <- rep(NA, length(projMCS7$cellNames))
names(cellTypeVector) <- projMCS7$cellNames

# Assign cell types
cellTypeVector[projMCS7$Clusters %in% c("C18", "C19", "C21")] <- "glut"
cellTypeVector[projMCS7$Clusters %in% c("C15", "C16", "C17", "C20", "C25")] <- "glutPrecursor"
cellTypeVector[projMCS7$Clusters %in% c("C22", "C23")] <- "gaba"
cellTypeVector[projMCS7$Clusters %in% c("C10", "C11")] <- "microglia"
cellTypeVector[projMCS7$Clusters %in% c("C2", "C3")] <- "oligo"
cellTypeVector[projMCS7$Clusters %in% c("C5", "C6")] <- "oligoPrecursor"
cellTypeVector[projMCS7$Clusters %in% "C4"] <- "glutOligo"
cellTypeVector[projMCS7$Clusters %in% c("C12", "C13", "C14")] <- "endo"
cellTypeVector[projMCS7$Clusters %in% "C8"] <- "astro"
cellTypeVector[projMCS7$Clusters %in% c("C1", "C7", "C9")] <- "astroPrecursor"
cellTypeVector[projMCS7$Clusters %in% "C24"] <- "glutAstro"

# Add this vector to the cellColData
projMCS7$cellTypes <- cellTypeVector


# Get marker features for all clusters in the ArchR project
markersGenes <- getMarkerFeatures(
  ArchRProj = projMCS7,
  useMatrix = "GeneScoreMatrix",
  groupBy = "cellTypes",  # Reference the new column
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


# Get significant marker list
markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

# Get the names of all cell types
celltype_names <- names(markerList)

# Loop through each cellType and save its marker list to a separate CSV file
for (cellType in celltype_names) {
  # Create a file name based on the cellType
  outputFileName <- paste0(cellType, "_GeneMarker_List_", Sys.Date(), ".csv")
  
  # Save the marker list for the current cellType to a CSV file
  write.csv(markerList[[cellType]], file = outputFileName)
  
  # Print message to track progress
  print(paste0("Saved marker list for cellType: ", cellType))
}


######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################


#glut markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_glut <- c("Mir684-1", "Olfr1499", "Olfr145", "Spata31d1a", "Fam183b", "Pate4",
                        "Gm6654", "Olfr535", "Olfr1441", "Spata46", "Olfr877", "Ptgs2", "Olfr1340",
                        "Serpina1f", "1700013D24Rik", "Klk1b3", "Gm11549", "Gm11938", "Olfr533", "Serpina11",
                        "Nrl", "Rprml", "Uts2", "D730048I06Rik", "Mir7055", "Exph5", "Nxf7", "Dgat2l6",
                        "Olfr209", "Gm17751", "4930452N14Rik", "Slc36a1os", "Ubxn10", "Neurod6",
                        "1700009C05Rik", "Pigr", "Nrgn", "Defb14", "Vmn1r199", "Cd37", "Krt17",
                        "Rin1", "Svs6", "Tead2", "Pla2g2c", "Rtn4r", "Ptcra", "Panx3", "Gfy", "Snai3",
                        "Defb34", "Ccdc155", "Lrrc74b", "Kcnv1", "Slc17a7", "2410137M14Rik", "Olfr1440",
                        "Nrn1", "Fap", "Vsig2", "Gsg1l2", "Sstr3", "1700072I22Rik", "Ptgs2os", "Mir128-1",
                        "AU023762", "4930532M18Rik", "Arr3", "Vsig1", "4921517D22Rik", "Rcvrn", "Cd6",
                        "4930405D11Rik", "Utp4", "Acot5", "Myf6", "P2rx6", "Gm11937", "Mir6957", "Flg2",
                        "B430306N03Rik", "Lingo1", "Glp2r", "Pth2", "Olfr878", "Bpifb3", "Pth2r", 
                        "4933407I05Rik", "Olfr1391", "Glt8d2", "Gm765", "Ackr2", "Prrg3", "Slamf1", "Ipcef1", 
                        "Slc6a16", "Ccm2l", "Islr", "Prdx6b", "A830019L24Rik", "Gm5771", "Serpina1e", "Gpr150", 
                        "Camk4", "Olfr1392", "Slc7a4", "Dmc1", "F630042J09Rik", "1700007P06Rik", "AI847159", 
                        "Olfr876", "Olfr1427", "Lrit2", "Urah", "A230070E04Rik", "Itpka", "Grap2", "Wfdc15a", 
                        "4930447F24Rik", "4933400B14Rik", "Mapk12", "Nptx1", "P2ry4", "Ddn", "2900079G21Rik", 
                        "Macc1", "Ocm", "Grm2", "Dlgap3", "Gpr26", "Tmco2", "Dnal4", "Fmo6", "Gm6642", "Agxt2", 
                        "Kdf1", "Cdsn", "Slc8a2", "Cdr1", "Myh2", "Spag11b", "A830010M20Rik", "Lrrc57", 
                        "Mir133a-1", "Fgf5", "Gm41341", "Adgrf4", "Chn1os3", "Fam163b", "Odaph", "Mir1900", 
                        "Ankrd33b", "Gm3279", "Mir6898", "A930017M01Rik", "Ddrgk1", "Tnfrsf25", "Prss40", 
                        "Gpr68", "Olfr577", "Serpina3h", "Mir215", "Krtap4-2", "Vipr1", "4930557J02Rik", 
                        "9830166K06Rik", "Ntf5", "Spn", "Pla2g2a", "4933428C19Rik", "Gm7694", "Clec14a", 
                        "2310075C17Rik", "Meltf", "Bloodlinc", "Mir6366", "Mybpc3", "Pcsk1", "Mir466k", 
                        "Mir7231", "Ccdc34", "Mir129-2", "Gm13710", "Tbxa2r", "Fcnb", "Ddo", "Nefm", "Ntsr1", 
                        "Acr", "A430093F15Rik", "Rnf148", "Myot", "Olfm1", "Itprid1", "Ces1f", "Ascl5", 
                        "Gm13030", "Pfkl", "Serinc2", "H2-M5")

#glutPrecursor markers: FDR <= 0.01 & Log2FC >=1.25, all 150 sorted by Log2FC
known_markers_glutPrecursor <- c("Olfr1200", "Fbxw16", "Serpina3b", "Krt90", "Skint1", "Gm6042", "Olfr786", 
                                 "Oas1e", "Vmn1r44", "Vmn1r200", "4933417A18Rik", "Olfr25", "Skint8", "Olfr1329", 
                                 "1700121N20Rik", "Tpsb2", "Krt82", "Mettl7b", "Oas3", "Krt4", "Skint7", "4933408N05Rik", 
                                 "4933406K04Rik", "Adgrg3", "Mir7007", "Olfr684", "Mup7", "AY702102", "Lrrc19", 
                                 "4933412E24Rik", "Olfr685", "Tpsab1", "Acaa1b", "Abca15", "Prop1", "1700029M20Rik", 
                                 "Olfr1467", "Cd244a", "Vmn1r-ps103", "Ffar3", "Prss29", "Tshz2", "Clec7a", "D630023F18Rik", 
                                 "Ebf4", "Htr4", "Ces1d", "Ms4a5", "Tcerg1l", "4933425L06Rik", "A630075F10Rik", "Olfr1448", 
                                 "Itln1", "Anks4b", "Olfr881", "Epsti1", "Nat3", "Vat1l", "Ctcflos", "Pcare", "Serpina3j", 
                                 "Phgdh", "Kcng4", "Gpr37", "Zfhx4", "Fa2h", "1700003F12Rik", "6430503K07Rik", "Fxyd3", 
                                 "Zcchc24", "Gab1", "Endou", "Cnp", "Mir6539", "Pde8a", "Olig2", "1190005I06Rik", "Txnip", 
                                 "Gm10804", "Itgb5", "Smad6", "Sall1", "Rab37", "A730020M07Rik", "Bcas1", "Gm15441", "Gm5627", 
                                 "Nkain1", "Ttyh2", "Amotl2", "Nat8f2", "Tgfa", "Angptl4", "Plpp3", "Eva1a", "Clic4", 
                                 "Plekho2", "4933433H22Rik", "Plekhf1", "Slc38a3", "Plekhg3", "Ptgds", "Calr4", "Rnase13", 
                                 "Nrarp", "Snx33", "8430430B14Rik", "Oaf", "Carhsp1", "St3gal4", "Plekhg1", "Zfp36l2", 
                                 "Mir219b", "Slc1a3", "Nat8", "D630033O11Rik", "Defb34", "Sox1", "3110099E03Rik", "Notch1", 
                                 "Pdk4", "Ermn", "Sox10", "Dio2", "Gm35978", "5430425K12Rik", "Prkd3", "Zic5", "Bco2", "Cd82", 
                                 "A230009B12Rik", "Kcnj10", "Foxs1", "Olfr761", "AI463170", "Cdc42ep1", "Olig1", "Fgfr2", 
                                 "Ankrd66", "5930403L14Rik", "Btnl1", "Vmn2r52", "Mir6996", "Folh1", "Vmn2r95", "Lgals2", 
                                 "Aspa", "4932411N23Rik", "Gng11", "Glt6d1")

#glutOligo markers: FDR <= 0.01 & Log2FC >=1.25, all 30 sorted by Log2FC
known_markers_glutOligo <- c("Pilrb2", "Oog3", "Olfr1118", "Krtap21-1", "Mucl1", "Mir133b", "Olfr830", "Chst4", 
                             "Tfpi2", "Gng11", "2310005A03Rik", "Ppp1r14a", "Ace2", "Cyp2j12", "Gm1979", "5430425K12Rik", 
                             "Folh1", "4930430A15Rik", "Olfr1116", "Sapcd1", "Scp2d1", "Nfe2l3", "Tmem100", "4930550C17Rik", 
                             "Rarres2", "Tbx3", "Mertk", "Gm16063", "Fzd2", "Cyp2f2")

#glutAstro markers: FDR <= 0.01 & Log2FC >=1.25, all 76 sorted by Log2FC
known_markers_glutAstro <- c("Vmn1r80", "Olfr131", "Vmn2r-ps129", "Vmn1r78", "Gm3219", "Krtap5-4", "Mlc1", 
                             "Olfr663", "Vmn1r79", "Prss45", "Hoxb6", "Arhgef16", "Gm9125", "Gfap", "BC016548", 
                             "Tmem52", "Usp17le", "Ceacam5", "Tmie", "Mir6236", "Als2cl", "Cyp2b19", "Gm2109", 
                             "Castor1", "Mov10l1", "Apoc1", "5430427M07Rik", "Mbl1", "5930403L14Rik", "Fzd2", 
                             "Atf7ip2", "Fgfr3", "Neurog1", "Lyg1", "Taar9", "9630013K17Rik", "Cyp2a22", "Sox21", 
                             "Prdm16", "Sox1", "Itih3", "Gm11681", "4930550C17Rik", "Chrdl1", "Tmem89", "Chad", 
                             "D630033O11Rik", "Cxcl14", "Mmrn2", "Gm4632", "LOC105245869", "Sp5", "Fam107a", "Myoz2", 
                             "Frmpd1os", "Endou", "Mirlet7c-2", "Vmn2r4", "Arx", "Hapln3", "Vmn2r103", "F630040K05Rik", 
                             "Gpr179", "Dmrta2", "Tril", "Atp1b2", "Tekt5", "Bcan", "2900052N01Rik", "Ffar2", "Trpv3", 
                             "Fbxw15", "Fa2h", "Opalin", "Trf", "Mag")

#gaba markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_gaba <- c("Krt5", "Mir7234", "Adad2", "Fbxw18", "Dlx6os2", "Slc32a1", "Olfr54", "Gm14204", "Gad1", 
                        "Dlx1", "Dlx6os1", "Mir3967", "Gad2", "Hist1h2ba", "Vmn1r217", "Olfr285", "Dlx1as", 
                        "Nipsnap3a", "Mir30f", "Dlx6", "Lhx6", "Crhbp", "Glrp1", "Slc5a4a", "Btbd11", "Krt15", 
                        "Htr3b", "Spp2", "Krt39", "Krt31", "Fam240b", "Gad1os", "Sytl5", "Dlx5", "Kcns3", "Mir376a", 
                        "1700025F24Rik", "Srpx2", "Mir30c-1", "Mir654", "Mir300", "1700120E14Rik", "Hoxd8", "Gm12408", 
                        "Twnk", "Ngf", "Pnoc", "Olfr282", "Ffar4", "Lrrc38", "Gm19589", "Erbb4", "S100z", "Cox6a2", 
                        "Csf2", "Cib4", "Kcnmb2", "Mir3107", "Mir376b", "Adamts15", "Cort", "Npas1", "Pcp4l1", 
                        "Kcnip1", "Krtap15", "Iltifb", "Npy", "4931431B13Rik", "Mir376c", "Grip1", "Gm4371", 
                        "4833422M21Rik", "Actr5", "Gm15421", "4930547E08Rik", "Rpp25", "4933432K03Rik", "Mettl21c", 
                        "Tmem132cos", "Grip1os1", "4930477N07Rik", "Nxph1", "Nr1h5", "Wfdc2", "Tmem233", "Rerg", 
                        "Reln", "Lincred1", "Hhla1", "Adarb2", "Kit", "Pvalb", "Vip", "Mir486", "Elfn1", "Zscan21", 
                        "Igf1", "Sh2d2a", "Rbms3", "Dynap", "Pldi", "Olfr1413", "Cacna2d2", "Th", "7420701I03Rik", 
                        "Grhl3", "Erich6b", "Minar1", "Tmem132c", "Metap1d", "Drd2", "Pde6c", "Rab3b", "Trim67", 
                        "1700016K05Rik", "Rbp4", "Gjd2", "Urad", "Zfp641", "Dcaf17", "Olfr284", "Akain1", "Actl7b", 
                        "Krt32", "Nuggc", "Kirrel", "Olfr270", "4930552N02Rik", "Gm20611", "Grin2d", "Nrip3", 
                        "Star", "Paqr5", "Sox6", "E130018O15Rik", "Pde5a", "Spert", "Ahsg", "Gm10845", "Gm15816", 
                        "Vat1", "1700031M16Rik", "Crtap", "Afap1", "Ncmap", "4933416E03Rik", "Grip1os2", "Sst", 
                        "Cracr2a", "6530411M01Rik", "1810063I02Rik", "Mir701", "C130080G10Rik", "Myadml2", "Gpr176", 
                        "Zim3", "Kdelr1", "Gphb5", "Ptprm", "Rnaset2b", "Aoc2", "Tdrd6", "Cntnap4", "Bend4", "Myo3a", 
                        "Tac2", "Gm5941", "BC030500", "Gm12339", "Htr3a", "Spag4", "Mir6408", "1700010J16Rik", 
                        "Kif26b", "Cers3", "Hykk", "Gm5833", "Speer4b", "Hnf4aos", "Tlr6", "Lrrd1", "1700054M17Rik", 
                        "Gsdme", "Myo5b", "Nppb", "Mir7027", "Acoxl", "Mir539", "Rnd2", "4930426D05Rik", "Ankrd55", 
                        "Vwc2", "Dpysl5", "4930447N08Rik", "Pycr1", "Lyz1", "8030423F21Rik", "Rnf207", "Mir3473e", "Drd1")

#oligo markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_oligo <- c("Pkd2l1", "Olfr357", "4930449I24Rik", "Gm5862", "Gm10471", "Fga", "Rnase11", "Olfr243", 
                         "Gm3402", "Gm6370", "Olfr624", "5031410I06Rik", "Olfr141", "Spink14", "Olfr116", "Olfr123", 
                         "Azgp1", "Krtap19-1", "Mir6917", "Mir7001", "Fa2h", "Vmn1r175", "2310005A03Rik", "Gm13023", 
                         "Krtap14", "Fbxw15", "Ppp1r14a", "Mir219a-2", "Sapcd1", "Erbb3", "Gm10863", "Vmn2r24", 
                         "Plp1", "Bpifc", "Gm5938", "Tshb", "Mal", "Olfr1331", "Folh1", "Pp2d1", "Vmn2r78", "Sox10", 
                         "Hoxd1", "Fbxw19", "Acsm3", "Vmn1r58", "A230009B12Rik", "Gm1979", "Olfr796", "Mir219b", 
                         "C1rb", "Ermn", "Rida", "Mir24-1", "Speer4a", "Trf", "Olfml1", "Gm15638", "Gjc3", "Car14", 
                         "Vmn2r33", "Glt6d1", "Gm14781", "Mir27b", "Aspa", "Smco3", "Vmn2r89", "Olfr292", "Mag", 
                         "Mir5135", "Opalin", "Rnase1", "Ankub1", "Gngt1", "Cd82", "Olfr1045", "Vmn2r77", "Gm5592", 
                         "4930405J17Rik", "Eri2", "Opn4", "Cts8-ps", "Cnp", "Hapln2", "Mrgprb2", "Scp2d1", "Olfr1491", 
                         "Gkn1", "Efhd1", "Ptgds", "Rpl10l", "Olfr1170", "Cltrn", "Ms4a4c", "Gng11", "Mir23b", "Lpo", 
                         "Olfr1333", "Gm13293", "Il25", "7630403G23Rik", "Gpr17", "Tfpi2", "Thsd1", "Pramel7", "Defb47", 
                         "B230206H07Rik", "A530053G22Rik", "Esp4", "Gm21671", "Rnf220", "Cldn11", "Tmem63a", "Pla2g16", 
                         "Sgk2", "Plxnb3", "Cyp2j11", "Cmtm5", "Olfr1093", "Clec4b1", "Gm5065", "I730030J21Rik", 
                         "5430425K12Rik", "Tktl2", "C030029H02Rik", "Rdh12", "Kcnrg", "Ugt8a", "Mir7239", "Cyp2j12", 
                         "Tmem125", "1500015L24Rik", "1700044K03Rik", "Vmn2r32", "Eda2r", "Vmn2r92", "Cyp27a1", "Ceacam5", 
                         "Pde8a", "Abca8a", "Defb30", "Scgb2b24", "Zscan4a", "Haglr", "Adam3", "Adam18", "Vmn1r76", 
                         "Fam71e2", "Vmn1r71", "Fgfr2", "Vmn1r60", "Cyp3a59", "Bpifa6", "Tspan2", "Rhox13", "Plin4", 
                         "9530014B07Rik", "Rnase6", "Nkx2-2", "Ifnk", "Cdc42ep1", "Galnt6", "Cyp2j9", "Olfr1082", 
                         "Il27", "Lrrc63", "Serpinb1a", "4930529C04Rik", "Carhsp1", "E330034G19Rik", "Vmn2r61", 
                         "Lims2", "Sult2a4", "Cdh19", "Selenop", "Vmn2r-ps60", "4930555G21Rik", "R3hcc1", "Cyp4a12b", 
                         "Trim13", "Prr18", "Psg-ps1", "4930544G11Rik", "AI463170", "Tppp3", "Thumpd1", "Nkx2-2os", 
                         "Nmral1", "Gjc2", "Gpr37", "Clca4a", "Mir6909", "Olfr52", "Art4", "Plekhf1", "2210409D07Rik", 
                         "Neil3", "Plpp2", "Efnb3", "Tnnt2")

#oligoPrecursor markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_oligoPrecursor <- c("Serpina3b", "Olfr157", "Olfr1051", "Gm38416", "Dcpp1", "Olfr103", "1110028F18Rik", 
                                  "Dcpp2", "Prss21", "Pdx1", "Orm2", "Tmem176a", "Doxl2", "Pdgfra", "C1ql1", "Gpr17", "Taf7", 
                                  "Cacng4", "Rbpjl", "H2-M10.3", "Matn4", "Mir759", "1600015I10Rik", "Tmem176b", "Gjc3", 
                                  "Adora2b", "Sun5", "E330017A01Rik", "Lad1", "Neu4", "Tmem255b", "Mif4gd", "Prss32", "Traf4", 
                                  "Myt1", "Mc5r", "9230009I02Rik", "Vmn2r71", "Gm4632", "Cyp4f37", "Fam89a", "Mir7025", "Gphb5", 
                                  "Pnlip", "Prl7b1", "Pgk2", "Lims2", "B3gnt7", "Gsx1", "4930505G20Rik", "Prl3c1", "Plut", 
                                  "Has2", "Sdr42e1", "A4gnt", "Olig2", "Mir7678", "Slc36a3", "Clba1", "Ly6g", "Ppp1r2-ps9", 
                                  "Mir3063", "Rgr", "Irx6", "A930009A15Rik", "Olfr1331", "Kank1", "4930563H07Rik", "Tgfa", 
                                  "1700021F07Rik", "Polr2f", "Rlbp1", "S100a3", "Xkr5", "S100a14", "Cspg4", "Has2os", "S100a16", 
                                  "Nkd2", "Slc22a18", "Chst5", "Esp31", "Ear10", "Slc2a10", "Cdrt4", "Zfp488", "C1ql2", "Gm813", 
                                  "Caskin2", "Gpr37l1", "Sycp1-ps1", "Lypd3", "Pmel", "Dcaf12l2", "Gpnmb", "Inava", "Ly6f", 
                                  "1110028F11Rik", "Mir5621", "Gm1553", "Lrit1", "Psg22", "Kif18b", "Acan", "Tnfsf9", "Atp6v1e2", 
                                  "Ear14", "Afap1l2", "6430503K07Rik", "Gm30052", "Kcnj16", "Mir6338", "Mybpc1", "Mir6996", 
                                  "Eprn", "Susd5", "Myo7b", "Gm13043", "Igf2bp1", "9230117E06Rik", "Spsb4", "Dok1", "Vwa2", 
                                  "Scrg1", "AA545190", "Fam69c", "Tat", "4930486F22Rik", "Cacng1", "Gm6961", "1700047L14Rik", 
                                  "Sox10", "Krt2", "Sox3", "Ntn1", "Fgfrl1", "Snord68", "Il25", "Sapcd2", "Bpifc", "Sbpl", 
                                  "Rpl13", "E2f8", "Sh3bp4", "Olig1", "Sox21", "Mmp2", "Gm15401", "Pilrb2", "Zfp957", "Spg7", 
                                  "Hsd17b2", "Gm13286", "Nipal1", "Zcchc24", "Habp2", "Klk9", "Cd209f", "Epn2", "Cyp2d10", 
                                  "Nat8f5", "4930444F02Rik", "F2r", "Col5a3", "Matn1", "P2rx2", "Meox1", "Ncmap", "Mir6406", 
                                  "1700044C05Rik", "Angpt4", "Dmrt3", "Oas1h", "Gm10863", "Klk5", "Pcmtd2", "Cmtm5", 
                                  "LOC105245869", "Gm2176", "Sema5b", "Lpo", "Stk32a", "Mir466f-4", "Xylt1", "Rnf43", "Mmp15", 
                                  "Halr1", "Gm15326", "Gm11110", "Ddi1", "Nkx2-2os", "Megf11", "Bcas1os2", "Krtap6-2", "Bcas1", 
                                  "Cdo1", "Plce1", "Mir7074", "Slc5a8", "P2ry1")

#microglia markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_microglia <- c("Gm4841", "Olfr433", "Mir6961", "Olfr693", "Clec4a3", "A530088E08Rik", "E230016K23Rik", 
                             "Il10ra", "Olfr1282", "Mir680-1", "Olfr1283", "Trim21", "Ccr1", "Cd300c2", "Olfr912", 
                             "Olfr1215", "C1qc", "Olfr108", "Ctla2b", "Bin2", "Cd33", "Nlrp1b", "Cx3cr1", "Olfr485", 
                             "Olfr520", "C5ar2", "Was", "BC147527", "Clec2i", "E230025N22Rik", "Olfr1232", "Olfr910", 
                             "Scgb1b7", "Itpripl1", "Susd3", "Tas2r108", "Olfr834", "Olfr461", "Csf1r", "Olfr110", "Ccr5", 
                             "Tmem119", "1600010M07Rik", "Ccl9", "Cd53", "Olfr1537", "Trim30b", "Vmn2r116", "Selplg", 
                             "Clec5a", "Olfr460", "Bcl2a1a", "Ccl6", "Olfr259", "Vmn1r171", "Clec9a", "Vsir", "Olfr1095", 
                             "Olfr1029", "Siglece", "Scgb1b29", "Ly6a", "Olfr1037", "Ifi203", "Cd52", "Serpinf1", "Olfr395", 
                             "Pilra", "2210414B05Rik", "Trim30a", "Gngt2", "Olfr1221", "Slc7a7", "Hhex", "4930515G16Rik", 
                             "Mir3079", "Mir3471-2", "Tubb1", "Cd300lf", "Mlxipl", "Mir6983", "Olfr107", "Dyrk4", "Mir7057", 
                             "Cd300ld3", "Cd5l", "AF067061", "Siglecf", "Olfr156", "Mir142", "Ikzf1", "Tlr6", "Fcgr1", 
                             "Ccr6", "Mir7086", "Olfr944", "Ccl12", "Gm5122", "Gif", "Olfr698", "Gm33619", "Olfr414", 
                             "Klri1", "Cd300ld4", "Cd300a", "Olfr1137", "H2-Ob", "Defa24", "Gtsf2", "Mmp12", "Cks1brt", 
                             "Il7r", "C1qb", "Mir6386", "Pfpl", "Clec2h", "Naip2", "Gm1966", "Siglecg", "Zfp658", 
                             "4930486L24Rik", "Vmn1r184", "Hpgd", "Ccr1l1", "Lair1", "4931402G19Rik", "Mrgprb1", "Gm11758", 
                             "Ccr2", "Mir223", "Ccl3", "Gm12185", "Ctss", "Gm826", "Pira6", "Gykl1", "Mmp3", "Adora3", 
                             "Tnfrsf17", "1700009J07Rik", "Rhoh", "Olfr1252", "Gm29811", "Gm5662", "Myo1f", "Prss37", 
                             "Olfr103", "Tifab", "Gpr160", "Lalba", "Olfr1193", "Lyn", "0610039K10Rik", "Gm5431", 
                             "Olfr109", "Mmp13", "Il6", "Olfr112", "Klrd1", "Mrgprx2", "Hcar2", "Milr1", "Ifnb1", "Fcgrt", 
                             "Gimap6", "Ncf1", "Gm5127", "Irf8", "Olfr111", "Fli1", "Ceacam14", "Csf3r", "Serpinb12", 
                             "Itgam", "Tlr11", "Lyz2", "Rps15a-ps4", "Ildr1", "Hexb", "Olfr963", "Olfr1028", "F11r", 
                             "Actl9", "Tgfb1", "Fcrls", "Ms4a6b", "Fcrla", "Irf5", "Tnni2", "Adap2", "Olfr694", "AI661453", 
                             "Samsn1", "1700024F13Rik", "Klri2", "Clec4a1", "Klhl6", "Bach2os", "Casp4", "Crybb1")

#endo-vasc markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_endo <- c("Olfr480", "Olfr1143", "Olfr508", "Saa1", "Olfr704", "Olfr723", "2310034C09Rik", "Olfr1042", 
                        "Vmn1r60", "Vmn2r12", "Gm13103", "Vmn1r179", "Olfr623", "Olfr1195", "Vmn1r83", "Gm5938", 
                        "Vmn2r28", "Olfr551", "Fpr-rs4", "Spata31d1b", "Aspn", "Olfr137", "Olfr727", "Olfr23", 
                        "Olfr393", "Olfr373", "Olfr1066", "Olfr724", "Zic4", "Mir7057", "Spaar", "Vmn2r105", 
                        "Vmn2r33", "BC021767", "Anpep", "Olfr139", "Vmn1r58", "Vmn2r16", "Olfr402", "Nlrp4a", 
                        "Vmn2r81", "Bcam", "Rbp1", "Omd", "C87977", "Olfr649", "Vmn1r172", "Vmn1r59", "Rem1", 
                        "Vmn1r178", "Higd1b", "Vmn1r230", "2310005E17Rik", "Gm2381", "Olfr782", "Sycp1-ps1", 
                        "Fam180a", "Scgb2b12", "Gm5779", "Zic1", "Vmn1r184", "Olfr792", "Mndal", "Mill2", "Amelx", 
                        "Vmn1r64", "Vmn2r32", "Jaml", "Cgnl1", "Eva1b", "Gprc5c", "Vmn1r74", "Olfr293", "Slc22a6", 
                        "Ifi211", "Olfr1220", "Cyp3a41a", "Ctla2b", "Sytl3", "Olfr985", "Slc6a20a", "Adap2", 
                        "Olfr301", "Olfr822", "Wfdc6a", "Olfr1461", "Vmn1r65", "Vmn2r14", "Slc6a13", "Olfr794", 
                        "Prf1", "Vmn2r108", "Gm3458", "Olfr821", "Gpr182", "Vmn1r176", "4933404G15Rik", "Olfr199", 
                        "Slfn5os", "Gm6756", "Gm36633", "Slc16a12", "Mrgprb2", "Olfr27", "Vmn2r13", "Flt1", 
                        "Cyp2c67", "Olfr19", "4930598F16Rik", "Saa2", "Vmn2r60", "Tfap2b", "Mrgprb8", "Mir503", 
                        "9230112D13Rik", "Olfr862", "Olfr310", "Serpina1c", "Rbpms", "Olfr1084", "Olfr366", "Fzd6", 
                        "Vmn1r50", "Vmn2r98", "Defb25", "Mrgprb1", "Olfr1252", "Olfr853", "Serpina1b", "Scgb1b19", 
                        "Uaca", "Slc30a10", "Tnfsf10", "Mir145a", "Slco1a4", "Mir3073a", "Ecm2", "Kcnj8", "Ocln", 
                        "Olfr304", "Tbx3", "Esp4", "Olfr891", "Orm2", "Atp2a3", "Psg22", "Vmn2r17", "Tmem45a", 
                        "Gm11827", "9430007A20Rik", "Mir483", "Olfr697", "Snai1", "Vmn2r115", "Sult2a7", "Olfr473", 
                        "Csn1s1", "Fcgrt", "Mir351", "Mir3109", "Vnn1", "Pramel4", "Nfatc4", "Olfr1253", 
                        "C730027H18Rik", "Aldh1a2", "Foxd1", "Mir322", "Igf2os", "Olfr938", "Olfr165", "Olfr1484", 
                        "Slc6a12", "Gm2083", "Gm27162", "Olfr780", "Olfr1234", "Plscr1", "Olfr1302", "Rab3d", 
                        "Spo11", "Scgb2b7", "Vmn2r15", "Cald1", "Olfr710", "Vmn2r106", "Asgr2", "Olfr912", "Pramef6", 
                        "Mrgprb4", "Vmn2r6", "Abcc9", "Nlrp9a", "AW549542", "Foxs1", "Rgs5", "Lama1", "Heyl", 
                        "Vmn1r73", "Arhgdib")

#astro markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_astro <- c("Gm4961", "Mir182", "Defb13", "Vmn1r64", "Olfr342", "Vmn1r65", "5930403L14Rik", "Phf11b", 
                         "Sytl3", "Gm3219", "Usp17la", "Ankrd66", "Vnn1", "Olfr281", "Scgb2b7", "Rida", "Ceacam3", 
                         "Sdsl", "Mir7222", "Fzd2", "Drc7", "Vmn1r71", "Tectb", "Vmn1r78", "Prss45", "Tma7", "Nlrp2", 
                         "Gm29508", "Prss43", "Endou", "Scara3", "4930550C17Rik", "F13b", "Zscan4c", "Mir6236", 
                         "Gjb2", "Neurog1", "Vmn2r107", "Fgfr3", "Sdc4", "Lcat", "Cyp4f15", "Prss50", "Pax8", "Sox1", 
                         "Itih3", "D630033O11Rik", "Mlc1", "Acsm5", "Olfr287", "Slc22a28", "Vmn2r13", "Tmem89", 
                         "Prdm16", "Olfr155", "Mir5098", "Olfr663", "Stat5a", "Spink11", "2900052N01Rik", "Sox1ot", 
                         "Olfr1349", "Ndrg2", "F630040K05Rik", "Slc7a10", "Gjb6", "Tdgf1", "Tppp2", "Klhdc7a", 
                         "Ranbp3l", "Gm32511", "Ceacam-ps1", "Syt15", "Gm11627", "Obox3", "Gm6083", "AA387883", 
                         "Ccdc51", "Sult2a4", "3110082J24Rik", "Cbs", "Als2cl", "Prss46", "Slc38a3", "Gm5089", 
                         "Gm20757", "Olfr770", "Garem2", "Prok1", "Cyp4f14", "Grin2c", "Rnase13", "G630093K05Rik", 
                         "Aqp4", "4930563H07Rik", "Gm13872", "Atp1b2", "Itih1", "Hils1", "Frmpd1", "Calr4", "Lfng", 
                         "Cmtm1", "Slc15a2", "Pax6", "Dio2", "Slc12a4", "Oplah", "Gm9125", "Phox2a", "Pla2g7", 
                         "Mir135a-2", "Castor1", "Fam107a", "C030037D09Rik", "Cxcl14", "Aldoc", "Frmpd1os", "Tril", 
                         "Vmn2r23", "Nr2e1", "4930432M17Rik", "Pou3f4", "Rapgef3", "Ccdc70", "Fzd9", "Mrgprb3", 
                         "Oaf", "Efcab3", "S1pr1", "Etnppl", "Mir7074", "Scrg1", "Arhgef19", "Atf7ip2", "Mir6993", 
                         "Vmn2r106", "1700012H19Rik", "Cyp2b19", "F3", "Slc25a18", "Gm11681", "5033404E19Rik", 
                         "Plxnb1", "Csf2", "Emx2os", "1700029P11Rik", "Hspb8", "Emx2", "Dmrta2", "Arhgef16", 
                         "Olfr1320", "Plpp3", "Tmie", "Foxn4", "Gpr37l1", "Plekhd1os", "Smpd2", "Zic5", "Hoxc12", 
                         "Cyp2j5", "Smpd5", "Gramd3", "9230102O04Rik", "0610039H22Rik", "BC061195", "Sox13", 
                         "Fam181a", "Eva1a", "Zscan4b", "Chad", "9630013K17Rik", "Acsf2", "Lgr6", "Hoxc9", "Slc9a3r1", 
                         "Actrt2", "Msi2", "Vmn1r61", "Grifin", "Gm17660", "Vmn2r45", "Usp17le", "Spatc1", 
                         "9230104L09Rik", "Dmp1", "Gm9268", "Ppp1r36", "5430427M07Rik", "Plin1", "Gpr179", "Col4a5", 
                         "Vmn2r61", "D930032P07Rik", "Bmx", "Hoxb5", "Nat8f6", "Mt2", "Tex12", "Obp2b")

#astroPrecursor markers: FDR <= 0.01 & Log2FC >=1.25, top 200 sorted by Log2FC
known_markers_astroPrecursor <- c("H2-Ea-ps", "Olfr699", "Olfr849", "Olfr298", "Olfr1256", "Olfr836", "Olfr835", 
                                  "Oog2", "Vmn1r10", "Olfr126", "Scgb1b12", "Olfr1136", "Olfr847", "Vmn1r171", "Iapp", 
                                  "Btnl2", "Olfr495", "Gm13088", "Klra13-ps", "Olfr811", "Psg26", "Vmn2r13", "Olfr1130", 
                                  "Vmn2r28", "Obox3", "Esp23", "Olfr866", "2310079G19Rik", "Olfr731", "Scgb2b26", "Olfr77", 
                                  "Cyp2a4", "Olfr497", "Ceacam3", "Vmn2r61", "Gm4884", "Olfr127", "Zim2", "Olfr362", "Vmn1r78", 
                                  "Acsm5", "Vmn1r172", "Prss45", "Olfr698", "Vmn2r12", "Lbx2", "Mrgprb4", "Vmn2r63", "Crisp3", 
                                  "1110028F18Rik", "Vmn1r61", "Adam26a", "Zscan4f", "Usp17ld", "Olfr159", "Clec4a4", "Mrgpra4", 
                                  "Mir133b", "Psg28", "Klk6", "Olfr730", "Ccl3", "Klk9", "BC016548", "Psg18", "Fga", "Nkx2-4", 
                                  "Vmn2r44", "Mageb1", "Vmn1r76", "Vmn1r75", "Scgb1b3", "Dyrk4", "Mir6236", "Klk8", "Snrpe", 
                                  "A630077J23Rik", "Tas2r121", "Skint4", "Vmn1r184", "Ccdc70", "Krtap16-3", "Ceacam5", "Defa24", 
                                  "Kif4-ps", "Nlrp4a", "Arhgef16", "Cox6b2", "Vmn1r64", "E230016K23Rik", "Vmn1r60", "Olfr502", 
                                  "Ugt2b37", "Vmn1r74", "Mrgpra3", "Gm8579", "Sult2a1", "Olfr134", "Ctsm", "Tmem82", "Vmn1r65", 
                                  "Mkrn3", "Vmn2r95", "Olfr293", "Gm9268", "Trim12a", "Scgb2b12", "Lyg1", "1700003P14Rik", 
                                  "Lpar4", "5430427O19Rik", "Scgb2b7", "Olfr487", "Mrc1", "5033404E19Rik", "BC064078", "Vmn2r96", 
                                  "Gm1979", "2810047C21Rik1", "F630028O10Rik", "Gm2061", "Sytl3", "Dio2", "Cts8", "Siglecf", 
                                  "Vmn2r23", "Mucl1", "Vmn2r20", "Slc7a10", "Vmn1r173", "Zscan4a", "Esp3", "Klhdc7a", "Cmtm1", 
                                  "Clec9a", "Sult2a7", "Gm16063", "Gm12887", "Ankrd66", "Nlrp9a", "Vmn2r45", "Trim34a", 
                                  "Mir135a-2", "Fam71e2", "Olfr1095", "Adam6a", "Vmn2r31", "Pmel", "Ranbp3l", "Prkag3", 
                                  "Serpina1a", "Cd209d", "Cts7", "Ecm2", "Hapln3", "Vmn1r77", "", "Cyp4f15", "Hotairm2", "Oaf", 
                                  "Igsf1", "Mir511", "Gm20757", "Tex24", "Paqr6", "Slc30a10", "Rida", "Olfr843", "Rbp1", 
                                  "5430427M07Rik", "Nat8", "4930550C17Rik", "Fgfr3", "Vmn1r72", "Nat8f2", "Zic5", "5930403L14Rik", 
                                  "Neurog1", "Cyp2j5", "Mir33", "Clec4a3", "Zfp36l2", "Rpl10l", "Vmn2r98", "Vmn2r22", "C87414", 
                                  "Slc25a34", "Vmn2r74", "5330439B14Rik", "Gm16287", "D230030E09Rik", "Scgb2b11", "Lgals1", 
                                  "Gm38423", "Esp8", "Scgb1b29", "Calr4", "Olfr829", "Heyl", "Tmprss12", "Fam107a")




######################################
######################################
######################################
######################################
######################################
######################################


## Subset projMCS7 by cell type

glut_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C18", "C19", "C21"), ]
glutPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C15", "C16", "C17", "C20", "C25"), ]
gaba_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C22", "C23"), ]
microglia_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C10", "C11"), ]
oligo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C2", "C3"), ]
oligoPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C5", "C6"), ]
glutOligo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C4"]
endo_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C12", "C13", "C14"), ]
astro_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C8"]
astroPrecursor_ArchRSubset <- projMCS7[projMCS7$Clusters %in% c("C1", "C7", "C9"), ]
glutAstro_ArchRSubset <- projMCS7[projMCS7$Clusters %in% "C24"]


############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################



## Complete following steps for each subset tx group:


# Specify which treatment group each sample is in:
# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+


treatment <- astro_ArchRSubset$Sample

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
astro_ArchRSubset$treatment <- treatment
# Check that this worked - if not, make sure the previous line was run successfully
head(astro_ArchRSubset$treatment)


##########################################
##########################################
##  Get marker genes for each cluster:
##########################################
##########################################


## Get marker features
markerGS_astro <- getMarkerFeatures(
  ArchRProj = astro_ArchRSubset,
  useMatrix = "GeneScoreMatrix",
  groupBy = "treatment",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

#Get the list of markers
markerList_astro <- getMarkers(markerGS_astro, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")


######################################
######################################
######################################
######################################
######################################
######################################



#Subset the markerGS object to just these known markers
se_idx <- which(rowData(markerGS_astro)$name %in% known_markers_astro)

#Get an index of these from the summarized experiment
markerGS_astro <- markerGS_astro[se_idx,]


##Plot it out:

library(viridis)

pdf(file = "astro_TxComp_Heatmap_2025-01-16.pdf", width = 15, height = 30)

plotMarkerHeatmap(
  seMarker = markerGS_astro,
  cutOff = "FDR <= 1 & Log2FC >= 0.01",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  #groupBy = "treatment",  # Group by treatment and order them based on the factor levels
  labelMarkers = known_markers_astro,  # Extracted marker names
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  pal = viridis(256),
  logFile = createLogFile("plotMarkerHeatmap")
)

dev.off()




######################################


# To print the Z-scores from the above heatmap

z_scores <- plotMarkerHeatmap(
  seMarker = markerGS_astro,
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
  labelMarkers = known_markers_astro,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  # Change this to TRUE
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

# Print z-scores from plotMarkerHeatmap
write.csv(z_scores, "astro_byTx_zscores_2025-01-16.csv", row.names = TRUE)


##########################################
##########################################
##########################################



## Plot descending zscores by the 2N treatment group 


# Read the data
mat <- read.csv("astro_byTx_zscores_2025-01-16.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert the data to a matrix
mat <- as.matrix(mat)

# Check the structure of your data
str(mat)

# Sort the data based on the t1 values in descending order
sorted_mat <- mat[order(mat[, "t1"], decreasing = TRUE), ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(sorted_mat)

color_palette <- colorRampPalette(c("#1F9E89FF", "white", "#440154FF"))(100)
# Create a color mapping for treatment groups 
annotation_colors <- list(Treatment = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF"))



# Open a JPEG device with calculated dimensions
jpeg(file = "/Volumes/DataBox/MCS2023/Tx_Comp/astro_byTx_sorted_heatmap_by_t1_2025-01-17.jpg", 
     width = 1200, 
     height = 5000, 
     res = 150)

pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Astro Sorted by 2N",
  color = color_palette,
  annotation_colors = annotation_colors,
  fontsize = 10, 
  cellwidth = 10, 
  cellheight = 10,
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  legend = FALSE            # Disable the legend
)

# Close the JPEG device to save the plot
dev.off()


##########################################
##########################################
##########################################


## Plot descending zscores by the Ts treatment group 


# Read the data
mat <- read.csv("astro_byTx_zscores_2025-01-16.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert the data to a matrix
mat <- as.matrix(mat)

# Check the structure of your data
str(mat)

# Sort the data based on the t3 values in descending order
sorted_mat <- mat[order(mat[, "t3"], decreasing = TRUE), ]

# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(c("2N", "2N+", "Ts", "Ts+"), levels = treatment_labels))
rownames(annotation_df) <- colnames(sorted_mat)


# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("#1F9E89FF", "white", "#440154FF"))(50)

# Create a color mapping for treatment groups 
annotation_colors <- list(Treatment = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF"))


# Open a JPEG device for saving the plot
jpeg(file = "/Volumes/DataBox/MCS2023/Tx_Comp/astro_byTx_sorted_heatmap_by_t3_2025-01-17.jpg", width = 1200, height = 5000, res = 150)

# Create the pheatmap without annotations and legends
pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Astro Sorted by Ts",
  color = color_palette,
  annotation_colors = annotation_colors,
  fontsize = 10, 
  cellwidth = 10, 
  cellheight = 10,
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  legend = FALSE            # Disable the legend
)

# Close the JPEG device to save the plot
dev.off()


##########################################
##########################################
##########################################


## Plot descending zscores by the 2N & 2N+ treatment groups


# Read the data and set the first column as row names
mat <- read.csv("astro_byTx_zscores_2025-01-16.csv", header = TRUE, row.names = 1)

# Convert the data frame to a numeric matrix
mat <- as.matrix(mat)

# Ensure the matrix is numeric
mat <- apply(mat, 2, as.numeric)  # Converts all columns to numeric

# Check the structure of the data to confirm it is numeric
str(mat)


# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(treatment_labels, levels = treatment_labels))
rownames(annotation_df) <- colnames(mat)

# Calculate the mean of t1 and t2 for sorting
mean_values <- rowMeans(mat[, c("t1", "t2")], na.rm = TRUE)

# Sort the data based on the mean of t1 and t2 values in descending order
sorted_mat <- mat[order(mean_values, decreasing = TRUE), ]

# Verify the sorted matrix
head(sorted_mat)



# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("#1F9E89FF", "white", "#440154FF"))(50)

# Create a color mapping for treatment groups 
annotation_colors <- list(Treatment = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF"))

# Open a JPEG device for saving the plot
jpeg(file = "/Volumes/DataBox/MCS2023/Tx_Comp/astro_byTx_sorted_heatmap_by_t1-t2_2025-01-17.jpg", width = 1200, height = 5000, res = 150)

# Create the pheatmap without annotations and legends
pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Astro Sorted by 2N & 2N+",
  color = color_palette,
  annotation_colors = annotation_colors,
  fontsize = 10, 
  cellwidth = 10, 
  cellheight = 10,
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  legend = FALSE              # Disable the legend
)

# Close the JPEG device to save the plot
dev.off()


##########################################
##########################################
##########################################


## Plot descending zscores by the Ts & Ts+ treatment group 


# Read the data
mat <- read.csv("astro_byTx_zscores_2025-01-16.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert the data to a matrix
mat <- as.matrix(mat)

# Check the structure of your data
str(mat)


# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(treatment_labels, levels = treatment_labels))
rownames(annotation_df) <- colnames(mat)

# Calculate the mean of t1 and t2 for sorting
mean_values <- rowMeans(mat[, c("t3", "t4")], na.rm = TRUE)

# Sort the data based on the mean of t1 and t2 values in descending order
sorted_mat <- mat[order(mean_values, decreasing = TRUE), ]

# Verify the sorted matrix
head(sorted_mat)


# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("#1F9E89FF", "white", "#440154FF"))(50)

# Create a color mapping for treatment groups 
annotation_colors <- list(Treatment = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF"))

# Open a JPEG device for saving the plot
jpeg(file = "/Volumes/DataBox/MCS2023/Tx_Comp/astro_byTx_sorted_heatmap_by_t3-t4_2025-01-17.jpg", width = 1200, height = 5000, res = 150)

# Create the pheatmap without annotations and legends
pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Astro Sorted by Ts & Ts+",
  color = color_palette,
  annotation_colors = annotation_colors,
  fontsize = 10, 
  cellwidth = 10, 
  cellheight = 10,
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  legend = FALSE              # Disable the legend
)

# Close the JPEG device to save the plot
dev.off()


##########################################
##########################################
##########################################


## Plot descending zscores by the 2N+ & Ts+ treatment group 


# Read the data
mat <- read.csv("astro_byTx_zscores_2025-01-16.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert the data to a matrix
mat <- as.matrix(mat)

# Check the structure of your data
str(mat)


# Define the updated treatment labels as a character vector
treatment_labels <- c("2N", "2N+", "Ts", "Ts+")  # New labels for treatments

# Create a data frame for column annotations with updated labels
annotation_df <- data.frame(Treatment = factor(treatment_labels, levels = treatment_labels))
rownames(annotation_df) <- colnames(mat)

# Calculate the mean of t2 and t4 for sorting
mean_values <- rowMeans(mat[, c("t2", "t4")], na.rm = TRUE)

# Sort the data based on the mean of t1 and t2 values in descending order
sorted_mat <- mat[order(mean_values, decreasing = TRUE), ]

# Verify the sorted matrix
head(sorted_mat)


# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("#1F9E89FF", "white", "#440154FF"))(50)

# Create a color mapping for treatment groups 
annotation_colors <- list(Treatment = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF"))

# Open a JPEG device for saving the plot
jpeg(file = "/Volumes/DataBox/MCS2023/Tx_Comp/astro_byTx_sorted_heatmap_by_t2-t4_2025-01-17.jpg", width = 1200, height = 5000, res = 150)

# Create the pheatmap without annotations and legends
pheatmap(
  sorted_mat,
  scale = "none",  # No scaling to maintain the original values
  cluster_cols = FALSE,  # Disable column clustering
  cluster_rows = FALSE,  # Disable row clustering to show sorted order
  show_colnames = TRUE,  # Show column names
  main = "Astro Sorted by 2N+ & Ts+",
  color = color_palette,
  annotation_colors = annotation_colors,
  fontsize = 10, 
  cellwidth = 10, 
  cellheight = 10,
  display_numbers = FALSE,  # Disable displaying numbers in the heatmap
  legend = FALSE              # Disable the legend
)

# Close the JPEG device to save the plot
dev.off()




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################



# Calculation of means between treatment groups 
## NOT pairwise comparisons



## Calculate significant genes between 2N & 2N+ vs. Ts & Ts+
# significant_(abs>=1.25)


# Load the data
data <- read.csv("astro_byTx_zscores_2025-01-16.csv")
names(data)
str(data)


# Calculate means for each group with descriptive column names
colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")  # rename first column to 'Gene'
data$"mean_phase1_(t1+t2)/2" <- (data$t1 + data$t2) / 2
data$"mean_phase2_(t3+t4)/2" <- (data$t3 + data$t4) / 2

# Calculate difference between means
data$"difference_phase1-phase2" <- data$"mean_phase1_(t1+t2)/2" - data$"mean_phase2_(t3+t4)/2"

# Add column to check if absolute difference between means is >= 1.25
data$"significant_(abs>=1.25)" <- abs(data$"difference_phase1-phase2") >= 1.25

# Save the updated data to a new CSV file
write.csv(data, "astro_sig_zscores_t1-t2vst3-t4_2025-01-17_analyzed.csv", row.names = FALSE)

# Check for significant differences
if (sum(data$"significant_(abs>=1.25)") == 0) {
  print(FALSE)  # Print FALSE if no significant differences are found
} else {
  # Print the subset of data with significant differences
  significant_data <- data[data$"significant_(abs>=1.25)" == TRUE, ]
  print(significant_data)
}



##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################



# Calculation of means between treatment groups 
## NOT pairwise comparisons



# List of cell types to loop over (replace with actual cell types)
cellTypes <- c("glut", "gaba", "microglia", "oligo", "endo", "astro")

## Not yet analyzed: "glutPrecursor", "oligoPrecursor", "glutOligo", "astroPrecursor", "glutAstro"


##########################################
##########################################
##########################################




## Loop for calculating differences between groups (t1+t2 vs t3+t4) by cellType:
# significant_(abs>=1.25)


# Loop through each cell type
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type (replace 'cellType' in filename)
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Calculate means for each group with descriptive column names
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")  # Rename first column to 'Gene'
  data$"mean_phase1_(t1+t2)/2" <- (data$t1 + data$t2) / 2
  data$"mean_phase2_(t3+t4)/2" <- (data$t3 + data$t4) / 2
  
  # Calculate difference between means
  data$"difference_phase1-phase2" <- data$"mean_phase1_(t1+t2)/2" - data$"mean_phase2_(t3+t4)/2"
  
  # Add column to check if absolute difference between means is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_phase1-phase2") >= 1.25
  
  # Save the updated data to a new CSV file
  output_filename <- sprintf("%s_t1-t2_t3-t4_zscores_2025-01-18_analyzed.csv", cellType)
  write.csv(data, output_filename, row.names = FALSE)
  
  # Check for significant differences
  if (sum(data$"significant_(abs>=1.25)") == 0) {
    print(sprintf("%s: No significant differences", cellType))  # Print if no significant differences are found
  } else {
    # Extract and print the subset of data with significant differences
    significant_data <- data[data$"significant_(abs>=1.25)" == TRUE, ]
    print(sprintf("%s: Significant differences found", cellType))
    print(significant_data)
    
    # Optionally, save the significant data to a separate file
    significant_output_filename <- sprintf("%s_significant_t1-t2_t3-t4_zscores_2025-01-18.csv", cellType)
    write.csv(significant_data, significant_output_filename, row.names = FALSE)
  }
}


##########################################
##########################################
##########################################



## Loop for calculating differences between groups (t2+t4 vs t1+t3) by cellType:



# Loop through each cell type
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type (replace 'cellType' in filename)
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Calculate means for each group with descriptive column names
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")  # Rename first column to 'Gene'
  data$"mean_group2_(t2+t4)/2" <- (data$t2 + data$t4) / 2
  data$"mean_group1_(t1+t3)/2" <- (data$t1 + data$t3) / 2
  
  # Calculate difference between means
  data$"difference_group2-group1" <- data$"mean_group2_(t2+t4)/2" - data$"mean_group1_(t1+t3)/2"
  
  # Add column to check if absolute difference between means is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_group2-group1") >= 1.25
  
  # Save the updated data to a new CSV file
  output_filename <- sprintf("%s_t2-t4_t1-t3_zscores_2025-01-18_analyzed.csv", cellType)
  write.csv(data, output_filename, row.names = FALSE)
  
  # Check for significant differences
  if (sum(data$"significant_(abs>=1.25)") == 0) {
    print(sprintf("%s: No significant differences", cellType))  # Print if no significant differences are found
  } else {
    # Extract and print the subset of data with significant differences
    significant_data <- data[data$"significant_(abs>=1.25)" == TRUE, ]
    print(sprintf("%s: Significant differences found", cellType))
    print(significant_data)
    
    # Optionally, save the significant data to a separate file
    significant_output_filename <- sprintf("%s_significant_t2-t4_t1-t3_zscores_2025-01-18.csv", cellType)
    write.csv(significant_data, significant_output_filename, row.names = FALSE)
  }
}


##########################################
##########################################
##########################################



## Loop for calculating differences between groups (t3 vs t4) by cellType:


# Loop through each cell type
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type (replace 'cellType' in filename)
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the difference between t3 and t4
  data$"difference_t3-t4" <- data$t3 - data$t4
  
  # Add a column to check if the absolute difference between t3 and t4 is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_t3-t4") >= 1.25
  
  # Save the updated data to a new CSV file with the updated filename
  output_filename <- sprintf("%s_t3_t4_zscores_2025-01-18_analyzed.csv", cellType)
  write.csv(data, output_filename, row.names = FALSE)
  
  # Check for significant differences
  if (sum(data$"significant_(abs>=1.25)") == 0) {
    print(sprintf("%s: No significant differences", cellType))  # Print if no significant differences are found
  } else {
    # Extract and print the subset of data with significant differences
    significant_data <- data[data$"significant_(abs>=1.25)" == TRUE, ]
    print(sprintf("%s: Significant differences found", cellType))
    print(significant_data)
    
    # Optionally, save the significant data to a separate file
    significant_output_filename <- sprintf("%s_significant_t3_t4_zscores_2025-01-18.csv", cellType)
    write.csv(significant_data, significant_output_filename, row.names = FALSE)
  }
}


##########################################
##########################################
##########################################



## Loop for calculating differences between groups (t3 vs t1, t2, t4) by cellType:


# Loop through each cell type
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the differences between t3 and the other groups
  data$"difference_t3_t1" <- data$t3 - data$t1  # t3 - t1
  data$"difference_t3_t2" <- data$t3 - data$t2  # t3 - t2
  data$"difference_t3_t4" <- data$t3 - data$t4  # t3 - t4
  
  # Add columns to check if the absolute difference is >= 1.25
  data$"significant_t3_t1_(abs>=1.25)" <- abs(data$"difference_t3_t1") >= 1.25
  data$"significant_t3_t2_(abs>=1.25)" <- abs(data$"difference_t3_t2") >= 1.25
  data$"significant_t3_t4_(abs>=1.25)" <- abs(data$"difference_t3_t4") >= 1.25
  
  # Save the updated data to a new CSV file
  output_filename <- sprintf("%s_t3-t1_t2_t4_zscores_2025-01-18_analyzed.csv", cellType)
  write.csv(data, output_filename, row.names = FALSE)
  
  # Check for significant differences for each comparison
  if (sum(data$"significant_t3_t1_(abs>=1.25)") == 0 &&
      sum(data$"significant_t3_t2_(abs>=1.25)") == 0 &&
      sum(data$"significant_t3_t4_(abs>=1.25)") == 0) {
    print(sprintf("%s: No significant differences for any comparison", cellType))  # Print if no significant differences are found
  } else {
    # Extract and print the subset of data with significant differences
    significant_data <- data[data$"significant_t3_t1_(abs>=1.25)" == TRUE |
                               data$"significant_t3_t2_(abs>=1.25)" == TRUE |
                               data$"significant_t3_t4_(abs>=1.25)" == TRUE, ]
    print(sprintf("%s: Significant differences found", cellType))
    print(significant_data)
    
    # Optionally, save the significant data to a separate file
    significant_output_filename <- sprintf("%s_significant_t3-t1_t2_t4_zscores_2025-01-18.csv", cellType)
    write.csv(significant_data, significant_output_filename, row.names = FALSE)
  }
}




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################



## Combine spreadsheets - Ts vs Ts+


# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the difference between t3 and t4
  data$"difference_t3_t4" <- data$t3 - data$t4
  
  # Add column to check if absolute difference between t3 and t4 is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_t3_t4") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_t3_t4")]
  
  # Rename the 'difference_t3_t4' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_Ts-Ts+", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences <- merge(combined_differences, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences, "Combined_CellType_significant_t3_t4_comparison_2025-01-18.csv", row.names = FALSE)


##########################################
##########################################
##########################################



#Combine spreadsheets - 2N+/Ts+ vs 2N/Ts

# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate means for group 1 and group 2
  data$"mean_group2_(t2+t4)/2" <- (data$t2 + data$t4) / 2
  data$"mean_group1_(t1+t3)/2" <- (data$t1 + data$t3) / 2
  
  # Calculate difference between means
  data$"difference_group2_group1" <- data$"mean_group2_(t2+t4)/2" - data$"mean_group1_(t1+t3)/2"
  
  # Add column to check if absolute difference between means is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_group2_group1") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_group2_group1")]
  
  # Rename the 'difference_group2_group1' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_2N+/Ts+_2N/Ts", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences <- merge(combined_differences, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences, "Combined_CellType_significant_t2-t4_t1-t3_comparison_2025-01-18.csv", row.names = FALSE)


##########################################
##########################################
##########################################



## Combine spreadsheets - 2N/2N+ vs Ts/Ts+


# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate means for group 1 and group 2
  data$"mean_group2_(t1+t2)/2" <- (data$t1 + data$t2) / 2
  data$"mean_group1_(t3+t4)/2" <- (data$t3 + data$t4) / 2
  
  # Calculate difference between means
  data$"difference_group2_group1" <- data$"mean_group2_(t1+t2)/2" - data$"mean_group1_(t3+t4)/2"
  
  # Add column to check if absolute difference between means is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_group2_group1") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_group2_group1")]
  
  # Rename the 'difference_group2_group1' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_2N_2N+_Ts_Ts+", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences <- merge(combined_differences, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences, "Combined_CellType_significant_t1-t2_t3-t4_comparison_2025-01-18.csv", row.names = FALSE)


##########################################
##########################################
##########################################


## t3 vs all - run each comparison seperately


## t3 vs t1

# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences_t3_vs_t1 <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the difference between t3 and t1
  data$"difference_t3_t1" <- data$t3 - data$t1
  
  # Add column to check if absolute difference between t3 and t1 is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_t3_t1") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_t3_t1")]
  
  # Rename the 'difference_t3_t1' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_T3_vs_T1", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences_t3_vs_t1) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences_t3_vs_t1 <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences_t3_vs_t1 <- merge(combined_differences_t3_vs_t1, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences_t3_vs_t1, "Combined_CellType_significant_t3_vs_t1_comparison_2025-01-18.csv", row.names = FALSE)


## t3 vs t2

# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences_t3_vs_t2 <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the difference between t3 and t2
  data$"difference_t3_t2" <- data$t3 - data$t2
  
  # Add column to check if absolute difference between t3 and t2 is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_t3_t2") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_t3_t2")]
  
  # Rename the 'difference_t3_t2' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_T3_vs_T2", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences_t3_vs_t2) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences_t3_vs_t2 <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences_t3_vs_t2 <- merge(combined_differences_t3_vs_t2, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences_t3_vs_t2, "Combined_CellType_significant_t3_vs_t2_comparison_2025-01-18.csv", row.names = FALSE)


## t3 vs t4

# Initialize an empty data frame to store all the combined differences across cellTypes
combined_differences_t3_vs_t4 <- data.frame()

# Loop through all cell types
for (cellType in cellTypes) {
  
  # Construct the filename for each cell type
  filename <- sprintf("%s_byTx_zscores_2025-01-16.csv", cellType)
  
  # Load the data for the current cell type
  data <- read.csv(filename)
  
  # Rename the first column to 'Gene' and keep t1, t2, t3, t4 as is
  colnames(data) <- c("Gene", "t1", "t2", "t3", "t4")
  
  # Calculate the difference between t3 and t4
  data$"difference_t3_t4" <- data$t3 - data$t4
  
  # Add column to check if absolute difference between t3 and t4 is >= 1.25
  data$"significant_(abs>=1.25)" <- abs(data$"difference_t3_t4") >= 1.25
  
  # Filter for only significant results
  significant_data <- data[data$"significant_(abs>=1.25)", ]
  
  # Select relevant columns for cellType differences (only significant ones)
  cellType_differences <- significant_data[, c("Gene", "difference_t3_t4")]
  
  # Rename the 'difference_t3_t4' column to reflect the cellType
  colnames(cellType_differences)[2] <- sprintf("%s_DA_T3_vs_T4", cellType)
  
  # Merge the current cellType's differences into the combined data frame
  if (nrow(combined_differences_t3_vs_t4) == 0) {
    # If it's the first iteration, just assign it directly
    combined_differences_t3_vs_t4 <- cellType_differences
  } else {
    # Merge based on the 'Gene' column
    combined_differences_t3_vs_t4 <- merge(combined_differences_t3_vs_t4, cellType_differences, by = "Gene", all = TRUE)
  }
}

# Export the combined data with only significant differences into a new spreadsheet
write.csv(combined_differences_t3_vs_t4, "Combined_CellType_significant_t3_vs_t4_comparison_2025-01-18.csv", row.names = FALSE)



## Combine t3 comparisons into 1 spreadsheet for downstream analysis


# Step 1: Read the data from each of the 3 spreadsheets (adjust file paths if needed)
data_t3_vs_t1 <- read.csv("Combined_CellType_significant_t3_vs_t1_comparison_2025-01-18.csv")
data_t3_vs_t2 <- read.csv("Combined_CellType_significant_t3_vs_t2_comparison_2025-01-18.csv")
data_t3_vs_t4 <- read.csv("Combined_CellType_significant_t3_vs_t4_comparison_2025-01-18.csv")

# Step 2: Merge the data frames based on the 'Gene' column
combined_data <- merge(data_t3_vs_t1, data_t3_vs_t2, by = "Gene", all = TRUE)
combined_data <- merge(combined_data, data_t3_vs_t4, by = "Gene", all = TRUE)

# Step 3: Check the resulting combined data
head(combined_data)

# Step 4: Save the combined data to a new CSV file
write.csv(combined_data, "Combined_CellType_significant_t3_vs_all_2025-01-19.csv", row.names = FALSE)




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################



## Tally total clusters that have the share DA gene markers between tx comps

# Repeat the following code for each combined cell type significant .csv file


# Load necessary library
library(dplyr)

# Load your data
file_path <- "Combined_CellType_significant_t1-t2_t3-t4_comparison_2025-01-18.csv"  
data <- read.csv(file_path)

# Count significant values for each gene
data_with_counts <- data %>%
  mutate(Significant_Count = rowSums(!is.na(select(., -Gene))))  # Count non-NA values for each gene

# Count the total significant entries for each column (excluding the Gene column)
total_significant_counts <- colSums(!is.na(select(data, -Gene)))

# Convert the counts to a data frame
total_counts_df <- as.data.frame(t(total_significant_counts))  # Transpose to make it a row
total_counts_df$Gene <- "Total"  # Add a Gene label for the total row

# Combine the original data with the total counts
final_data <- bind_rows(data_with_counts, total_counts_df)

# Save the updated data frame to a new CSV file
output_file_path <- "Combined_CellType_Total-sig_t1-t2_t3-t4_comparison_2025-01-19.csv"  # Update this to your desired output path
write.csv(final_data, output_file_path, row.names = FALSE)

# Print the results
print(final_data)


################################################


## Tally up the total number of genes in the "Significant_Count" column

# Run for each set of treatment group comparisons


# Load necessary libraries
library(dplyr)

# Load your data from the CSV file
file_path <- "Combined_CellType_Total-sig_t3_t4_comparison_2025-01-19.csv"  # Update this to your file path
data <- read.csv(file_path)

# Count occurrences of 1's, 2's, and 3's
tally <- data %>%
  count(Significant_Count) %>%
  filter(Significant_Count %in% c(1, 2, 3))

# Print the tally
print(tally)


## Results from above tallies:

Combined_CellType_Total-sig_t1-t2_t3-t4_comparison_2025-01-19.csv
Significant_Count   n
                 1 320
                 2   3
                 
Combined_CellType_Total-sig_t2-t4_t1-t3_comparison_2025-01-19.csv
Significant_Count   n
                 1 303
                 2   3

Combined_CellType_Total-sig_t3_t4_comparison_2025-01-19.csv
Significant_Count   n
                 1 548
                 2   8        
                 
Combined_CellType_Total-sig_t3_vs_all_comparison_2025-01-19.csv
Significant_Count   n
                 1 610
                 2 318
                 3 126

