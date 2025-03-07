#Setup an interactive session
#salloc --account=eon -t 0-06:00 --mem=128G --nodes=2 --ntasks-per-node=16
#Updated conda env 12-2023
#module load miniconda3/23.1.0
#conda activate archr2023_12

#Load libraries
#R
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
library(tidyr)
library(cowplot)

#Additional setup
setwd("/Volumes/DataBox/MCS2023/Tx_Comp/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

projMCS5
#NumberOfCells: 104455
#medianTSS: 14.521
#medianFrags: 6653

getAvailableMatrices(projMCS5)
# "GeneScoreMatrix" "PeakMatrix"      "TileMatrix" 

############################################

# To compare samples amongst tx groups in a particular cluster:

t1ArchRSubset <-projMCS5[
        projMCS5$Sample %in% c("C302_", "C306_", "C309_", "C318_", "C323_", "C328_", "C332_", "C337_", "C339_",
        "C346_", "C351_", "C353_", "C360_"),]
# C323 is 2N (not Ts)

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

##############################################

t4markers <- getMarkerFeatures(
   ArchRProj = t4ArchRSubset,
   useMatrix = "GeneScoreMatrix",
   groupBy = "Clusters",
   testMethod = "wilcoxon",
   bias = c("TSSEnrichment",
   "log10(nFrags)"),
)

t4markers

##############################################

#To make comparisons and plots among specific treatments; use getMarkerFeatures to identify differentially accessible pe$
#Using groupBy = "Sample" returns empty individual cluster .csv files

#t1markers <- getMarkerFeatures(
#   ArchRProj = t1ArchRSubset,
#   useMatrix = "GeneScoreMatrix",
#   groupBy = "Sample",
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment",
#   "log10(nFrags)"),
#)

#t1markers

##############################################

table(projMCS5$Clusters)

#   C1   C10   C11   C12   C13   C14   C15   C16   C17   C18   C19    C2   C20 
#  257   975  2908   872   409  2136   560  1677  1916 44381   976   180  3426 
#  C21   C22   C23   C24   C25    C3    C4    C5    C6    C7    C8    C9 
# 10686  6750  3513   874   396 10313  1005    51  2027    61  8020    86 

table(t1ArchRSubset$Clusters)
#   C1   C10   C11   C12   C13   C14   C15   C16   C17   C18   C19    C2   C20 
#   49   252   732   223   100   519   148   388   441 10697   247    41   747 
#  C21   C22   C23   C24   C25    C3    C4    C5    C6    C7    C8    C9 
# 2437  1599   858   176    97  2428   239    13   518    14  1982    24 

table(t2ArchRSubset$Clusters)
#   C1   C10   C11   C12   C13   C14   C15   C16   C17   C18   C19    C2   C20 
#   75   264   907   263   123   689   147   479   515 12758   294    42   944 
#  C21   C22   C23   C24   C25    C3    C4    C5    C6    C7    C8    C9 
# 3027  1917   965   248   118  2954   254    15   574    25  2530    22 

table(t3ArchRSubset$Clusters)
#   C1   C10   C11   C12   C13   C14   C15   C16   C17   C18   C19    C2   C20 
#   76   218   767   222   103   503   149   445   521 10908   221    39   931 
#  C21   C22   C23   C24   C25    C3    C4    C5    C6    C7    C8    C9 
# 2715  1646   850   233    88  2668   243    12   482    11  1849    23 

table(t4ArchRSubset$Clusters)
#  C1  C10  C11  C12  C13  C14  C15  C16  C17  C18  C19   C2  C20  C21  C22  C23 
#  50  219  381  142   69  364   99  310  360 8460  174   51  673 2104 1329  712 
# C24  C25   C3   C4   C5   C6   C7   C8   C9 
# 184   77 1915  244   11  376    7 1403   14 

#####

t4_markerList <- getMarkers(t4markers, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25", #returnGR = TRUE)
)

t4_markerList

# 2024-04-08 thresholds: FDR <= 0.01 & abs(Log2FC) >= 1.25
# 2024-04-10 thresholds: FDR <= 0.01 & abs(Log2FC) >= 1.25

t4_markerList$C1

#####

write.csv(t1_markerList, file = "T1_GeneMarker_Master_GrpByCluster_2024-04-10.csv")
# All individual cluster .csv files are empty when using groupBy = "Sample"
write.csv(t1_markerList$C1, file = "T1_C1_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C2, file = "T1_C2_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C3, file = "T1_C3_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C4, file = "T1_C4_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C5, file = "T1_C5_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C6, file = "T1_C6_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C7, file = "T1_C7_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C8, file = "T1_C8_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C9, file = "T1_C9_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C10, file = "T1_C10_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C11, file = "T1_C11_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C12, file = "T1_C12_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C13, file = "T1_C13_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C14, file = "T1_C14_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C15, file = "T1_C15_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C16, file = "T1_C16_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C17, file = "T1_C17_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C18, file = "T1_C18_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C19, file = "T1_C19_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C20, file = "T1_C20_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C21, file = "T1_C21_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C22, file = "T1_C22_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C23, file = "T1_C23_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C24, file = "T1_C24_GeneMarkers_2024-04-10.csv")
write.csv(t1_markerList$C25, file = "T1_C25_GeneMarkers_2024-04-10.csv")

#####

write.csv(t2_markerList, file = "T2_GeneMarker_Master_GrpByCluster_2024-04-10.csv")
# All individual cluster .csv files are empty when using groupBy = "Sample"
write.csv(t2_markerList$C1, file = "T2_C1_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C2, file = "T2_C2_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C3, file = "T2_C3_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C4, file = "T2_C4_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C5, file = "T2_C5_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C6, file = "T2_C6_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C7, file = "T2_C7_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C8, file = "T2_C8_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C9, file = "T2_C9_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C10, file = "T2_C10_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C11, file = "T2_C11_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C12, file = "T2_C12_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C13, file = "T2_C13_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C14, file = "T2_C14_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C15, file = "T2_C15_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C16, file = "T2_C16_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C17, file = "T2_C17_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C18, file = "T2_C18_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C19, file = "T2_C19_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C20, file = "T2_C20_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C21, file = "T2_C21_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C22, file = "T2_C22_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C23, file = "T2_C23_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C24, file = "T2_C24_GeneMarkers_2024-04-10.csv")
write.csv(t2_markerList$C25, file = "T2_C25_GeneMarkers_2024-04-10.csv")

#####

write.csv(t3_markerList, file = "T3_GeneMarker_Master_GrpByCluster_2024-04-10.csv")
# All individual cluster .csv files are empty when using groupBy = "Sample"
write.csv(t3_markerList$C1, file = "T3_C1_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C2, file = "T3_C2_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C3, file = "T3_C3_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C4, file = "T3_C4_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C5, file = "T3_C5_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C6, file = "T3_C6_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C7, file = "T3_C7_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C8, file = "T3_C8_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C9, file = "T3_C9_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C10, file = "T3_C10_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C11, file = "T3_C11_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C12, file = "T3_C12_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C13, file = "T3_C13_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C14, file = "T3_C14_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C15, file = "T3_C15_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C16, file = "T3_C16_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C17, file = "T3_C17_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C18, file = "T3_C18_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C19, file = "T3_C19_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C20, file = "T3_C20_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C21, file = "T3_C21_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C22, file = "T3_C22_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C23, file = "T3_C23_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C24, file = "T3_C24_GeneMarkers_2024-04-10.csv")
write.csv(t3_markerList$C25, file = "T3_C25_GeneMarkers_2024-04-10.csv")

#####

write.csv(t4_markerList, file = "T4_GeneMarker_Master_GrpByCluster_2024-04-10.csv")
# All individual cluster .csv files are empty when using groupBy = "Sample"
write.csv(t4_markerList$C1, file = "T4_C1_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C2, file = "T4_C2_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C3, file = "T4_C3_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C4, file = "T4_C4_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C5, file = "T4_C5_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C6, file = "T4_C6_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C7, file = "T4_C7_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C8, file = "T4_C8_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C9, file = "T4_C9_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C10, file = "T4_C10_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C11, file = "T4_C11_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C12, file = "T4_C12_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C13, file = "T4_C13_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C14, file = "T4_C14_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C15, file = "T4_C15_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C16, file = "T4_C16_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C17, file = "T4_C17_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C18, file = "T4_C18_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C19, file = "T4_C19_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C20, file = "T4_C20_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C21, file = "T4_C21_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C22, file = "T4_C22_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C23, file = "T4_C23_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C24, file = "T4_C24_GeneMarkers_2024-04-10.csv")
write.csv(t4_markerList$C25, file = "T4_C25_GeneMarkers_2024-04-10.csv")

##############################################################################
##############################################################################

markersGenes <- getMarkerFeatures(
    ArchRProj = projMCS5,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersGenes

markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25", #returnGR = TRUE)
)
markerList

# 2024-04-08 thresholds: FDR <= 0.01 & abs(Log2FC) >= 1.25
# 2024-04-26 thresholds: FDR <= 0.01 & abs(Log2FC) >= 1.25

markerList$C1

#####

write.csv(markerList, file = "GeneMarker_List_Master_2024-04-08.csv")
write.csv(markerList$C1, file = "C1_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C2, file = "C2_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C3, file = "C3_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C4, file = "C4_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C5, file = "C5_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C6, file = "C6_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C7, file = "C7_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C8, file = "C8_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C9, file = "C9_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C10, file = "C10_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C11, file = "C11_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C12, file = "C12_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C13, file = "C13_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C14, file = "C14_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C15, file = "C15_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C16, file = "C16_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C17, file = "C17_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C18, file = "C18_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C19, file = "C19_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C20, file = "C20_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C21, file = "C21_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C22, file = "C22_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C23, file = "C23_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C24, file = "C24_GeneMarker_List_2024-04-08.csv")
write.csv(markerList$C25, file = "C25_GeneMarker_List_2024-04-08.csv")

##############################################################################
##############################################################################

#Get lengths of markerlist for each cluster
mlistlens <- lapply(markerList, nrow)

#Remove any clusters with zero markers
markerList <- markerList[mlistlens != 0]


filt_num <- 200
top_genes <- lapply(markerList, function(x){ # for each element of markerList
  if(nrow(x) >=filt_num){ # if there are 5 entries
    x[1:filt_num,]$name   # get the top 5
  }else{ # if there are not 50 entries
    x[1:nrow(x),]$name # get all the entries
  }
})

#Get all the genes unlisted and only the unique ones
unique(unlist(top_genes))
length(unique(unlist(top_genes)))

str(top_genes)
#write.table(top_genes, file = "Top_200_Genes_Master_2024-04-08.csv")
## Error in running command due to difference in row numbers
#read.table(top_genes)
## Error in read.table(top_genes, header = TRUE, sep = "\t") : 


assays(markersGenes)

#Take a look at how this is encoded
markersGenes@assays@data$Mean[1:10,1:3]
markersGenes@assays@data$MeanDiff[1:10,1:3]
markersGenes@assays@data$FDR>0

#Get out mean and relative accessibility
mean_acc <- markersGenes@assays@data$Mean
rel_acc  <- markersGenes@assays@data$MeanDiff
FDR_acc <- markersGenes@assays@data$FDR

#Add rownames
rownames(mean_acc) <- rowData(markersGenes)$name
rownames(rel_acc) <- rowData(markersGenes)$name
rownames(FDR_acc) <- rowData(markersGenes)$name

#Add column names
colnames(mean_acc) <- colnames(markersGenes)
colnames(rel_acc) <- colnames(markersGenes)
colnames(FDR_acc) <- colnames(markersGenes)

#Melt it
melted_mean_acc <- reshape2::melt(as.matrix(mean_acc))
melted_rel_acc <- reshape2::melt(as.matrix(rel_acc))
melted_FDR_acc <- reshape2::melt(as.matrix(FDR_acc))
colnames(melted_mean_acc) <- c("gene", "cluster", "mean")
colnames(melted_rel_acc) <- c("gene", "cluster", "rel")
colnames(melted_FDR_acc) <- c("gene", "cluster", "FDR")
melted_mean_acc$gene <- as.character(melted_mean_acc$gene)
melted_rel_acc$gene <- as.character(melted_rel_acc$gene)
melted_FDR_acc$gene <- as.character(melted_FDR_acc$gene)

#Merge em together
acc_comb1 <- merge(melted_mean_acc, melted_rel_acc)
acc_comb2 <- merge (acc_comb1, melted_FDR_acc)
tail(acc_comb2)

################################################
################################################

#https://davemcg.github.io/post/lets-plot-scrna-dotplots/

markers <- unique(unlist(top_genes))

markers

#Marker_ID <- c(
#  "Olig2", "Neurod6", "Gad2", "Dio2", "C1qa", "Anpep"
#)
#markers_1 <-c (
#  "Vip", "Sst", "Pvalb", "Pdgfra", "Olig2", "Olig1", "Neurod6",
#  "Neurod2", "Neurod1", "Mog", "Gfap", "Gad2", "Gad1", "Drd2",
#  "Drd1", "Dio2", "Cx3cr1", "C1qa", "Anpep"
#)


## This doesn't subset anything
# pdf(file="test_dotplot_smh.pdf", width=10, height=10)  
# ggplot(acc_comb2, aes(x=cluster, y = gene, color = rel, size = mean)) + 
#   geom_point() 
# dev.off()

#Remove dots where there is zero (or near zero) expression
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point()

#Better color, better theme, rotate x axis labels
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  scale_color_viridis_c(name = 'log2 (count + 1)') +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

#Tweak color scaling
acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

acc_comb2

#########################################

# make data square to calculate euclidean distance
mat <- acc_comb2 %>%
  filter(gene %in% markers) %>%
  # filter(mean > 0, abs(rel) > 1.0) %>%
  select(-mean, -FDR) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = rel) %>%
  data.frame() # make df as tibbles -> matrix annoying

row.names(mat) <- mat$gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
mat <- as.matrix(mat)
clust <- hclust(dist(mat)) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

dotplot <- acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

cowplot::plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

dotplot <- acc_comb2 %>% filter(gene %in% markers) %>%
  filter(mean > 0, abs(rel) > 1.0) %>%
  #gene = factor(gene, levels = clust$labels[clust$order]) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = gene, color = rel, size = mean)) +
  geom_point() +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

cowplot::plot_grid(ggtree_plot, NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

#########################################
#########################################
#https://davemcg.github.io/post/lets-plot-scrna-dotplots/

#To group genes (cluster) by similar expression patterns and show the dendrogram

#hclust 
#mutate(gene = factor(markers, levels = YOURNEWGENEORDER))


#################################################
#################################################

pdf(file="GeneMarker_Top50_Heatmap_2024-03-08.pdf", width=15, height=6)
plotMarkerHeatmap(   ### Copied this over from the full function reference - leaving anything not commented at defaul$
  seMarker = markersGenes, ## Set this to the correct object containing the markers
  cutOff = "FDR <= 0.01 & abs(Log2FC >= 2)",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  #Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)
dev.off()

##############################################

hmap_to_subset <- plotMarkerHeatmap(  #Anything not commented left as default
  seMarker = markerGS,   ##Set this to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  ##Label the markers of interest to plot
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  K = 5,
  logFile = createLogFile("plotMarkerHeatmap")
)

## Get the matrix of the heatmap to figure out row indices for the genes of interest

hmap_to_subset_MAT <- plotMarkerHeatmap(   #Anything not commented left as default
  seMarker = markerGS,   ##Set to the correct object containing the markers
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,  #Exclude everything we're not interested in right now
  pal = NULL,
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = unique(unlist(top_genes)),  ##Label the markers of interest we're plotting here
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = TRUE,  #Return the matrix rather than a heatmap
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)


hmap_row_idx <- which(rownames(hmap_to_subset_MAT) %in% unique(unlist(top_genes))) #Get the row numbers of the heatmap matrix that match the known markers

#Plot it out
pdf(file="GeneMarkers_Top50_1-25-2024.pdf", width=7, height=5)
hmap_to_subset[hmap_row_idx, ]
dev.off()
