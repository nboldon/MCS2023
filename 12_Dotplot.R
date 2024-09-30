
setwd("/Volumes/DataBox/ProjMCS6")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

##################################################

projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)
projMCS6
getAvailableMatrices(projMCS6)

##################################################

#Load libraries
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork) 

##################################################
##################################################

#To get a dataframe of the mean and relative accessibility of each gene for each cluster:

Cluster.markers <- getMarkerFeatures(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

assays(Cluster.markers)

#Take a look at how this is encoded
Cluster.markers@assays@data$Mean[1:10,1:3]
Cluster.markers@assays@data$MeanDiff[1:10,1:3]
Cluster.markers@assays@data$FDR>0

#Get out mean and relative accessibility
mean_acc <- Cluster.markers@assays@data$Mean
rel_acc  <- Cluster.markers@assays@data$MeanDiff
FDR_acc <- Cluster.markers@assays@data$FDR
Log2FC_acc <- Cluster.markers@assays@data$Log2FC

#Add rownames
rownames(mean_acc) <- rowData(Cluster.markers)$name
rownames(rel_acc) <- rowData(Cluster.markers)$name
rownames(FDR_acc) <- rowData(Cluster.markers)$name
rownames(Log2FC_acc) <- rowData(Cluster.markers)$name

#Add column names
colnames(mean_acc) <- colnames(Cluster.markers)
colnames(rel_acc) <- colnames(Cluster.markers)
colnames(FDR_acc) <- colnames(Cluster.markers)
colnames(Log2FC_acc) <- colnames(Cluster.markers)

#Melt it
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc))
melted_rel_acc <- reshape2::melt(as.matrix(rel_acc))
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc))
melted_Log2FC_acc <- reshape2::melt(as.matrix(Log2FC_acc))
colnames(melted_mean_acc) <- c("gene", "cluster", "mean")
colnames(melted_rel_acc) <- c("gene", "cluster", "rel")
colnames(melted_FDR_acc) <- c("gene", "cluster", "FDR")
colnames(melted_Log2FC_acc) <- c("gene", "cluster", "Log2FC")
melted_mean_acc$gene <- as.character(melted_mean_acc$gene)
melted_rel_acc$gene <- as.character(melted_rel_acc$gene)
melted_FDR_acc$gene <- as.character(melted_FDR_acc$gene)
melted_Log2FC_acc$gene <- as.character(melted_Log2FC_acc$gene)

#Merge em together
###adding melted_FDR_acc did not work
###Error: Error in fix.by(by.x, x) : 
###'by' must specify one or more columns as numbers, names or logical

acc_comb1 <- merge(melted_mean_acc, melted_rel_acc)
acc_comb2 <- merge (acc_comb1, melted_FDR_acc)
acc_comb3 <- merge(acc_comb2, melted_Log2FC_acc)
tail(acc_comb3)

################################################
################################################
################################################

#write.table(projMCS6, file = "/project/eon/nboldon/MCS2023/ClusterID/merged.tsv", 
#            append = FALSE, quote = TRUE, sep = "\t",
#            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#            col.names = TRUE, qmethod = c("escape", "double"),
#            fileEncoding = "")

###############################################
###############################################
###############################################

#gene_cluster <- read_tsv("C:\\Users\\Naomi\\Desktop\\TIN_CANN\\TIN_CANN.tsv")

#gene_cluster %>% sample_n(5)

###############################################
#https://davemcg.github.io/post/lets-plot-scrna-dotplots/

markers <- acc_comb3$gene %>% unique()

# Zu et al., 2023 marker genes 

# GABAergic (Inhibitory)
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

interesting_markerGenes <- c( "Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Gad1", "Aqp4", "Dio2", "Cx3cr1", "Olig1", "Opalin")

#########################################
#########################################
#########################################
#########################################
#########################################


## https://davemcg.github.io/post/lets-plot-scrna-dotplots/#how-do-i-make-a-dotplot
# gene_cluster %>% filter(Gene %in% markers) %>% 
# mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#   filter(count > 0, `% Expressing` > 1) %>% 
#   ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#   geom_point() 


#Remove dots where there is zero (or near zero) expression
acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(FDR) > 0.01) %>% 
  ggplot(aes(x=cluster, y = gene, color = FDR, size = mean)) + 
  geom_point() 

acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(Log2FC) > 0.5) %>% 
  ggplot(aes(x=cluster, y = gene, color = Log2FC, size = mean)) + 
  geom_point() 

acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(rel) > 1.0) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) + 
  geom_point() 

# gene_cluster %>% filter(Gene %in% markers) %>% 
# mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#  filter(count > 0, `% Expressing` > 1) %>% 
#  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#  geom_point() + 
#  scale_color_viridis_c(name = 'log2 (count + 1)') + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank()) 

#Better color, better theme, rotate x axis labels
acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(rel) > 1.0) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'MeanDiff') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(Log2FC) > 0.5) %>% 
  ggplot(aes(x=cluster, y = gene, color = Log2FC, size = mean)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'log2FC') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(FDR) > 0.5) %>% 
  ggplot(aes(x=cluster, y = gene, color = FDR, size = mean)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'FDR') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

## One gene is at much higher expression than anything else. 
## I would rather have anything at say log(count+1) ~ 4 be bright yellow. 
## We can do this by using scale_color_gradientn with limits set to c(0,4) 
## and have anything above 4 be “squished” down by oob = scales::squish
# gene_cluster %>% filter(Gene %in% markers) %>% 
# mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#  filter(count > 0, `% Expressing` > 1) %>% 
#  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#  geom_point() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank()) +
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

#Tweak color scaling
acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(rel) > 1.0) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), 
                        oob = scales::squish, name = 'MeanDiff')
                  

acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(FDR) > 0.5) %>% 
  ggplot(aes(x=cluster, y = gene, color = FDR, size = mean)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0.25, 1.00), 
                        oob = scales::squish, name = 'FDR')


acc_comb3 %>% filter(gene %in% interesting_markerGenes) %>% 
  filter(mean > 0, abs(Log2FC) > 1.0) %>% 
  ggplot(aes(x=cluster, y = gene, color = Log2FC, size = mean)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(-2,4), 
                        oob = scales::squish, name = 'Log2FC')

acc_comb3

#########################################

## Group the genes (cluster) by similar expression patterns and show the dendrogram. 

# Option: run hclust and reorder the genes by the results with mutate(Gene = factor(Gene, levels = YOURNEWGENEORDER). 
# This will visually (gene order) get the right result, but there’ll be no dendrogram.

## Create dendrogram:
# make data square to calculate euclidean distance
# mat <- gene_cluster %>% 
#  filter(Gene %in% markers) %>% 
#  select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
#  pivot_wider(names_from = cluster, values_from = count) %>% 
#  data.frame() # make df as tibbles -> matrix annoying
#row.names(mat) <- mat$Gene  # put gene in `row`
#mat <- mat[,-1] #drop gene column as now in rows
#clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
#ddgram <- as.dendrogram(clust) # create dendrogram
#ggtree_plot <- ggtree::ggtree(ddgram)
#ggtree_plot


# make data square to calculate euclidean distance
mat <- acc_comb3 %>% 
  filter(gene %in% oligo_markerGenes) %>%
  #filter(mean > 0, abs(rel) > 1.0) %>%
  select(-mean, -FDR, -Log2FC) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = rel) %>% 
  data.frame() # make df as tibbles -> matrix annoying

row.names(mat) <- mat$gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
mat <- as.matrix(mat)
clust <- hclust(dist(mat)) # hclust with distance matrix

ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

## Glue them together with cowplot
#Notice how rel_widths in plot_grid is used to tweak the relative width of each plot 
#and we are using align to attempt to line the plots up.
#dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
#  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
#  filter(count > 0, `% Expressing` > 1) %>% 
#  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#  geom_point() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank()) +
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
#plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

dotplot <- acc_comb3 %>% filter(gene %in% astro_markerGenes) %>% 
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

## Reorder the genes with the hclust ordering.
# And SQUEEZE the plots together with a cowplot trick of adding a fake plot in between and giving it a negative distance.
#dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
#  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
#         Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
#  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#  geom_point() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank()) +
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
#plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

dotplot <- acc_comb3 %>% filter(gene %in% astro_markerGenes) %>% 
  filter(mean > 0, abs(rel) > 1.0) %>%
  gene = factor(gene, levels = clust$labels[clust$order]) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')


## Move the gene names to the right side with scale_y_discrete(position = "right")
#dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
#  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
#        Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
#  filter(count > 0, `% Expressing` > 1) %>% 
#  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#  geom_point() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank()) +
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
#  scale_y_discrete(position = "right")
#plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

dotplot <- acc_comb3 %>% filter(Gene %in% astro_markerGenes) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
         Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
  scale_y_discrete(position = "right")





pdf(file="Ggtree_Zu-Glut2_2024-07-17.pdf", width=13, height=17)
plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')
dev.off()


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


## Define custom colors for groups

# Define your custom color palette
custom_colors <- c(
  "Group1" = "#1f77b4",  # Blue
  "Group2" = "#ff7f0e",  # Orange
  "Group3" = "#2ca02c",  # Green
  "Group4" = "#d62728",  # Red
  "Group5" = "#9467bd"   # Purple
  # Add more colors as needed
)

## Create custom x-axis labels with these colors

# Map clusters to new group labels
cluster_to_group <- c(
  "Cluster1" = "Group1",
  "Cluster2" = "Group1",
  "Cluster3" = "Group2",
  "Cluster4" = "Group3",
  "Cluster5" = "Group4",
  "Cluster6" = "Group5"
  # Add more mappings as needed
)

# Create a dataframe with x-axis labels and their corresponding colors
x_axis_labels <- data.frame(
  label = names(cluster_to_group),
  group = cluster_to_group,
  color = custom_colors[cluster_to_group],
  x = seq_along(cluster_to_group) - 0.5  # Adjust position if necessary
)

## Reorder and update x-axis labels

# Create a custom annotation function for x-axis labels
custom_x_axis <- function(labels_df, y_position) {
  grobs <- lapply(seq_len(nrow(labels_df)), function(i) {
    textGrob(
      labels_df$label[i],
      x = labels_df$x[i],
      y = y_position,
      just = "center",
      gp = gpar(col = labels_df$color[i], fontsize = 12)  # Adjust fontsize as needed
    )
  })
  grobTree(children = do.call(gList, grobs))
}

# Example dotplot with updated x-axis
dotplot <- filtered_data %>%
  ggplot(aes(x = cluster, y = gene, color = cluster, size = mean)) + 
  geom_point() + 
  scale_color_manual(values = custom_colors) + 
  theme_cowplot() + 
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + 
  labs(x = '', y = 'Gene')

# Add custom x-axis labels
dotplot <- dotplot + annotation_custom(
  custom_x_axis(x_axis_labels, y_position = -0.1)  # Adjust y_position if needed
)

# Example ggtree plot (replace with your actual ggtree plot)
ggtree_plot <- ggtree(clust) + 
  theme_tree() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# Combine the plots
combined_plot <- plot_grid(
  ggtree_plot,
  dotplot,
  nrow = 2,
  rel_heights = c(1, 2)
)

# Save the combined plot
ggsave("combined_plot_with_custom_x_axis.png", plot = combined_plot, width = 12, height = 10, dpi = 300)


#######################################################################################
#######################################################################################

# Generate colors using viridis
num_colors <- length(unique(filtered_data$cluster))
# Here, num_colors should match the number of distinct clusters in your plot. 
viridis_colors <- viridis(num_colors)
# viridis(num_colors) generates a palette with a range of colors.

# Map colors to labels
# Define your new cluster labels
new_labels <- unique(filtered_data$cluster)

# Map clusters to colors
color_mapping <- setNames(viridis_colors, new_labels)

# Apply colors to x-axis labels
# Create a dataframe for x-axis labels with colors
x_axis_labels <- data.frame(
  label = new_labels,
  color = viridis_colors,
  x = seq_along(new_labels) - 0.5  # Adjust position as needed
)

# Create a dataframe for x-axis labels with colors
x_axis_labels <- data.frame(
  label = new_labels,
  color = viridis_colors,
  x = seq_along(new_labels) - 0.5  # Adjust position as needed
)

# Define a custom annotation function for x-axis labels
custom_x_axis <- function(labels_df, y_position) {
  grobs <- lapply(seq_len(nrow(labels_df)), function(i) {
    textGrob(
      labels_df$label[i],
      x = unit(labels_df$x[i], "npc"),
      y = unit(y_position, "npc"),
      gp = gpar(col = labels_df$color[i])
    )
  })
  do.call(gList, grobs)
}

# Define your new cluster labels and colors
new_labels <- c("1", "2", "3", "4")  # Customize as needed
cluster_colors <- c("Cluster A" = "blue", "Cluster B" = "red", "Cluster C" = "green", "Cluster D" = "purple")  # Customize as needed

# Reorder dendrogram based on your preferences
order <- c("C1", "C3", "C2", "C4")  # Define your preferred order
mat_reordered <- mat[, order]  # Apply the new order

ggtree_plot <- ggtree(hclust(dist(mat_reordered))) +
  theme_tree() +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),  # Adjust margins
    axis.text.x = element_blank(),            # Remove x-axis labels
    axis.text.y = element_blank(),             # Remove y-axis labels
    geom_text2(aes(label=group_label), hjust=-.3)  # Replace `group_label` with your label column
  ) +
  geom_tiplab(size = 0)                      # Remove tip labels (if any)

# Create a dot plot with ordered clusters and custom labels
dotplot <- filtered_data %>%
  mutate(cluster = factor(cluster, levels = order, labels = new_labels),
         gene = factor(gene, levels = rownames(mat_reordered))) %>%
  ggplot(aes(x=cluster, y=gene, color=rel, size=mean)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'log2 (count + 1)') + 
  scale_x_discrete(labels = new_labels) +
  theme_cowplot() + 
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=8),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_text(color = "black")) +  # Adjust as needed
  labs(x = 'Cluster', y = 'Gene') +  # Add axis labels here
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2))


# Create and save the plot
ggsave("dotplot_with_custom_x_axis.png", plot = dotplot + 
         annotation_custom(
           grob = custom_x_axis(x_axis_labels, 0),
           xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
         ), width = 10, height = 6, dpi = 300)

# Combine the plots
combined_plot <- (ggtree_plot / dotplot) + plot_layout(nrow = 2, heights = c(1, 3))

# Save the combined plot
ggsave("combined_plot.png", plot = combined_plot, width = 8, height = 10, dpi = 300)

#######
#######

#Or save by:
  # Save the combined plot as a PNG file
  ggsave("combined_plot.png", plot = combined_plot, width = 10, height = 6, dpi = 300)

# Save the combined plot as a PDF file
ggsave("combined_plot.pdf", plot = combined_plot, width = 10, height = 6)

## OR 

# Open a JPEG device
jpeg("astro_dotPlot_2024-07-22.jpeg", width = 500, height = 200, quality = 90)

# Print the combined plot
print(combined_plot)

# Close the device
dev.off()

#########################################
#########################################
#########################################

