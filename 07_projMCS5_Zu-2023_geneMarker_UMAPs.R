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
library(dplyr)
library(tidyr)
library(cowplot)
library(ggtree)
library(aplot)
library(patchwork)
library(ggplot2)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/ClusterID/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

getAvailableMatrices(projMCS6)

table(projMCS6$Clusters)

############################################
############################################

projMCS6 <- addHarmony(
  ArchRProj = projMCS6,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Use MAGIC to impute gene scores by smoothing signal across nearby cells
projMCS6 <- addImputeWeights(projMCS6)
getImputeWeights(projMCS6)

############################################
############################################

# Assigning clusters with gene scores

projMCS6 <- addUMAP(ArchRProj = projMCS6, reducedDims = "Harmony", force = TRUE)

############################################

# Overlay marker gene scores in 2D UMAP embedding

########

# Zu et al., 2023
# Isocortex gene markers

gaba_markerGenes <- c(
"2610028E06Rik", "4930467D21Rik", "A830036E02Rik", "Acvr1c", "Adamts12", "Adamtsl5",
  "Adgrg6", "Ahrr", "Akr1c18", "Ano3", "Anxa2", "Ap1s3", "Aqp5", "Arhgap31", "Atp6ap1l",
  "Barx2", "C1qtnf7", "Cadps2", "Calb1", "Calb2", "Cbln4",
  "Cd24a", "Cd34", "Cd44", "Cdh23", "Cdh7", "Cdk6", "Chat", "Chrna2", "Chrnb3",
  "Chst7", "Cntn6", "Cntnap5a", "Cobll1", "Col14a1", "Col15a1", "Col19a1", "Col24a1",
  "Col25a1", "Corin", "Cort", "Cox6a2", "Cp", "Cplx3", "Cpne8", "Crabp1", "Creb5", "Crh",
  "Crhr2", "Crim1", "Crlf1", "Crybg1", "Csgalnact1", "Cxcl14", "Cyp1b1", "D030045P18Rik",
  "Dach2", "Dmrt2", "Dnah5", "Dock5", "Drd1", "Dysf", "Ecel1",
  "Edaradd", "Edn1", "Edn3", "Egfr", "Egln3", "Egr2", "Erbb4", "Eya1", "Eya4", "F830208F22Rik",
  "Fam129a", "Fam163a", "Fbn2", "Fgfr3", "Fibin", "Fign", "Fstl4", "Fxyd6", "Gad1", "Gad2",
  "Galntl6", "Gas2l3", "Gdpd2", "Gldn", "Glis3", "Glp1r", "Gm5087", 
  "Gpc3", "Gpc4", "Gpc5", "Gpc6", "Gpr149", "Gpx3", "Grm8", "Gulp1", "Has2os", "Hpse",
  "Hs6st2", "Hs6st3", "Hsd11b1", "Htr2c", "Htr3a", "Htr4", "Ifit1", "Igfbp6", "Il7", "Irs4",
  "Itga4", "Itga8", "Itih5", "Kcns3", "Kit", "Klhl14", "Krt12", "Krt73", "Lama3", "Lamp5",
  "Lgals1", "Lhx6", "Lonrf3", "Lpl", "Lrrc38", "Mc5r", "Met", "Mfge8", "Mgat4c", "Moxd1",
  "Mybpc1", "Myh13", "Myh4", "Myh8", "Myl1", "Myo5b", "Ndnf", "Necab1", "Nell1", "Ngf", "Npas1",
  "Npnt", "Nptx2", "Ntn1", "Ntng1", "Nts", "Nxph1", "Nxph2", "Oprd1", "Ostn", "Otogl", "Oxtr",
  "Paqr5", "Pax6", "Pcdh18", "Pcsk5", "Pde11a", "Pde5a", "Pdgfra", "Pdyn", "Pdzrn3", "Pla2g5",
  "Platr14", "Pnoc", "Postn", "Pou6f2", "Ppp1r1c", "Prdm1", "Prdm8", "Prkg2", "Prokr2", "Prss12",
  "Prss23", "Ptger3", "Ptgs2", "Pth2r", "Pthlh", "Ptprk", "Ptprt", "Pvalb", "Qrfpr", "Rasl11a",
  "Rassf2", "Rbp4", "Reln", "Ret", "Rgs5", "Rgs6", "Rph3al", "Rprml", "Rspo2", "Rtl4", "Rxfp1",
  "Rxfp3", "Scube2", "Sema3e", "Sertm1", "Sfrp2", "Sfta3-ps", "Shisa9", "Shisal2b", "Slc24a4",
  "Slc5a7", "Slc9a2", "Slit1", "Sln", "Snx31", "Sorcs1", "Sostdc1", "Spp1", "Sst", "St8sia6",
  "Stac", "Stk32a", "Stk32b", "Stxbp6", "Syndig1l", "Syt10", "Syt2", "Syt6", "Tac1", "Tac2",
  "Tacr3", "Tcap", "Teddm3", "Tent5a", "Tgfbi", "Th", "Thbs2", "Tmem100", "Tmem26", "Tnfaip8l3",
  "Tnni3k", "Tnnt1", "Tph2", "Tpm2", "Unc13c", "Vip", "Vwc2", "Zfp804b"
)
#Error: FeatureName (4930544I03Rik,5033421B08Rik,5330429C05Rik,9630002D21Rik,A830012C17Rik,B020031H02Rik,
#C87487,Ccn3,D030055H07Rik,Dchs2,Gm14507,Gm16685,Gm26673,Gm28175,Gm30373,Gm32442,Gm39185,Gm40518,Gm43154,Gm49227) does not exist!

glut_markerGenes <- c(
  "0610040J01Rik", "1700023F02Rik", "2310002F09Rik", "4833415N18Rik",
  "4930407I19Rik", "4930467D21Rik", "4933406J09Rik", "5330416C01Rik", "9130024F11Rik",
  "9330158H04Rik", "9430014N10Rik", "A630023P12Rik", "A830009L08Rik", "A830029E22Rik",
  "A830036E02Rik", "Adam33", "Adamts1", "Adamts18", "Adamts2", "Adarb2", "Adgrd1", "Adgrg6", "Adra1a",
  "Aldh1a3", "Aldh1l1", "Alkal1", "Ankfn1", "Ano3", "Arhgap25", "Arhgap28", "Arhgap31", "Atp10a", "Atp2b4",
  "B430212C06Rik", "Barx2", "Batf3", "Baz1a", "BC006965", "Bcl11b", "Blnk", "Bmp5", 
  "Bmpr1b", "Brinp3", "C030013G03Rik", "C1ql2", "C1ql3", "C7", "Calb1", "Camk2d", "Car3", "Cartpt", "Cbln1",
  "Cbln2", "Cbln4", "Ccbe1", "Ccdc3", "Ccdc60", "Ccdc80", "Cck", "Ccnb1", "Cd24a",
  "Cd34", "Cd36", "Cdc25c", "Cdh10", "Cdh12", "Cdh13", "Cdh18", "Cdh20", "Cdh8", "Cdh9", "Cdhr1", "Cdk6",
  "Cemip", "Chn2", "Chrm3", "Chrna6", "Chst9", "Clec18a", "Clic5", "Cntn5", "Cntnap3", "Cntnap5a", "Col15a1",
  "Col23a1", "Col24a1", "Col5a1", "Col5a2", "Col6a1", "Colq", "Cort", "Cpa6", "Cplx3", "Cpm", "Cpne5", "Cpne7",
  "Cpne9", "Crh", "Crispld1", "Crispld2", "Cryab", "Crym", "Csrnp1", "Ctsc", "Cxcl12", "Cxcl14", "Cyp26b1",
  "Cyp39a1", "Ddit4l", "Diaph3", "Dio3", "Dkk2", "Dkkl1", "Dlgap2", "Dmrt2", "Dpp10", "Drd1", "Drd5", "Ednra",
  "Egfem1", "Egln3", "Egr2", "Egr3", "Endou", "Epas1", "Erg", "Fbxl7", "Fezf2", "Fgf10", "Fgf11", "Fign", "Flt3",
  "Fmn1", "Fosl2", "Foxp2", "Fst", "Fxyd6", "Gadl1", "Galnt14", "Gfra1", "Gipc2", "Glis1", "Glra3", "Gm10635",
  "Gm10754", "Gm11549", "Gm12295", "Gm20752", "Gm5820", "Gm6260", "Gm7271", "Gpc6")

#Error: FeatureName (1700047F07Rik,2600014E21Rik,4921539H07Rik,9630002D21Rik,B530045E10Rik,Ccn2,Ccn3,Ccn4,Gm11639,
#Gm11730,Gm12031,Gm12840,Gm12930,Gm13391,Gm13601,Gm13974,Gm15270,Gm17501,Gm19410,Gm20713,Gm30094,Gm30371,Gm30382,
#Gm31698,Gm32647,Gm34184,Gm34567,Gm35853,Gm39822,Gm40331,Gm40518,Gm41414,Gm42707,Gm45623,Gm47033,Gm48530,Gm49678,
#Gm49906,Gm50370,Gm6209,Gm9732) does not exist!

glut2_markerGenes <- c("Gpr139", "Gpr88", "Gria1", "Grik1", "Grik3", "Grin2a", "Grm8", "Grp", "Gxylt2", "Hapln1", "Hgf", "Hkdc1", "Hmga2",
  "Hpgd", "Hs3st4", "Hspb3", "Htra1", "Igfbp2", "Igfbp4", "Igfbp6", "Igfn1", "Il1rapl2", "Inhba", "Inhbb",
  "Inpp4b", "Iqgap2", "Itga9", "Itprid1", "Kcnh1", "Kcnh5", "Kcnh7", "Kcnj16", "Kcnt2", "Kctd8", "Klf5", "Klhl1",
  "Klhl13", "Krt80", "Krt9", "Kynu", "L3mbtl4", "Lama4", "Layn", "Lcp1", "Lepr", "Lgr5", "Lhfp", "Lhx2", "Lingo2",
  "Lipg", "Lpl", "Lrrc55", "Lrrtm4", "Ly6g6e", "Lypd1", "Macc1", "Maf", "Mafb", "Masp1", "Matn4",
  "Mc4r", "Mctp2", "Megf11", "Met", "Mgat4c", "Mgp", "Moxd1", "Myl4", "Myzap", "Ndnf", "Necab1", "Nell1", "Neurod6",
  "Ngf", "Nkx1-2", "Nos1", "Npnt", "Npsr1", "Nptx2", "Npy", "Npy2r", "Nr2f2", "Nr4a2", "Nrgn", "Nrn1", "Nrtn",
  "Ntn5", "Nxph3", "Olfr111", "Opn3", "Osr1", "Otof", "Otx1", "Ovol2", "P2rx5", "Pamr1", "Pappa", "Pappa2", "Pcdh15",
  "Pcdh19", "Pcp4", "Pde7b", "Pdia5", "Pdyn", "Pdzrn4", "Penk", "Pld5", "Plekha2", "Plekhd1", "Pln", "Plod2", "Plpp4",
  "Podxl", "Pou3f1", "Pou6f2", "Prdm8", "Prkch", "Prkcq", "Prkg2", "Prlr", "Prph", "Prss22", "Psd4", "Ptgfr", "Ptgs2",
  "Ptprd", "Ptpre", "Ptpru", "Pvalb", "Qrfpr", "Rab38", "Ramp3", "Rasl10a", "Reln", "Rgs16", "Rgs5", "Rgs8", "Rims3",
  "Robo3", "Rorb", "Rprm", "Rrad", "Rspo1", "Rxfp1", "Rxfp2", "S100b", "S1pr3", "Samd3", "Satb2", "Scml4", "Scn4b",
  "Scn5a", "Scn7a", "Scnn1a", "Scube2", "Sdk1", "Sema3d", "Sema3e", "Serpina3n", "Serpine2", "Sertm1", "Sgcd", "Sgcz",
  "Sh2d4b", "Sh3bp4", "Shisa6", "Shroom3", "Sla2", "Slc17a6", "Slc17a7", "Slc17a8", "Slc30a3", "Slc5a5", "Slit2", 
  "Smoc2", "Sostdc1", "Sowahb", "Spata13", "Spc25", "Spink8", "St8sia2", "Stac", "Stard13", "Stard8", "Stat4",
  "Stk17b", "Stxbp6", "Sulf1", "Susd5", "Sv2c", "Syndig1", "Synpr", "Syt17", "Syt2", "Syt6", "Tacr1", 
  "Tcap", "Tcerg1l", "Teddm3", "Tent5a", "Tgfbr2", "Tll2", "Tmem117", "Tmem125", "Tmem215", "Tmem40", "Tnfaip6", "Tnmd")

#Error: FeatureName (Ighm,Lncbate10,Lratd2,Sox5os4,Tafa1,Tafa2) does not exist! 

astro_markerGenes <- c(
"Cxcl5", "Dab1", "Ddn", "Fam163a", "Gja1", "Kcnq1ot1", "Kcnq3", "S1pr1"
)

oligo_markerGenes <- c(
"9630013A20Rik", "Bmp4", "Cck", "Cenpa", "Cenpf", "Foxg1", "Gp1bb", "Gpr17", "Mki67", "Mms22l", "Olig1", "Pclaf", "Pdgfra",
  "Pdzd2", "Plekhg1", "Ptger1", "Top2a", "Zfp36l1"
)

endoVasc_markerGenes <- c(
"Adgrl4", "Alpl", "Grm7", "Igf2", "Ly6c1", "Myl9", "Olfr558", "Rgs5", "Slc22a8", "Slco1a4", "Srpx2", "Tln2", "Vtn"
)

microglia_markerGenes <- c(
"C1qa", "Tmem119"
)

############################################

p <- plotEmbedding(
  ArchRProj = projMCS6, 
  colorBy = "GeneScoreMatrix", 
  name = microglia_markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.99),
  imputeWeights = getImputeWeights(projMCS6)
)

# To plot a specific gene, subset the plot list using the gene name

p$Olig1

# To plot all genes, use cowplot to arrange the different plots together.
# Each of these marker genes lights up the corresponding cell clusters. 

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 5),p2))

# Save an editable vectorized version of the plots

plotPDF(plotList = p, 
        name = "microglia_geneMarkers-Zu-2023_2024-07-10.pdf", 
        ArchRProj = projMCS6, 
        addDOC = FALSE, width = 5, height = 5)


##################################################
##################################################

# To plot UMAP by Clusters for comparisons

p <- plotEmbedding(
  ArchRProj = projMCS6,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

plotPDF(p, name = "projMCS6_UMAP-Clusters_2024-07-10.pdf", ArchRProj = projMCS6, addDOC = FALSE, width = 5, height = 5)
