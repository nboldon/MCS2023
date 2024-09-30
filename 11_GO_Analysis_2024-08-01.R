library(Seurat)
library(Signac)
library(scran)
library(BiocManager)
library(BiocGenerics)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ArchR)
library(pheatmap)
set.seed(1)

#Additional setup
setwd("/Volumes/DataBox/GO_Analysis/")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

projMCS6
getAvailableMatrices(projMCS6)

###################################################
###################################################
###################################################


Cluster.markers <- getMarkerFeatures(
  ArchRProj = projMCS6,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 4000
)


markerList.Cluster.markers <- getMarkers(Cluster.markers,cutOff = "FDR <= 0.1 & abs(Log2FC) >=0.5")


###################################################
###################################################
###################################################


#To obtain top 200 markers in each cluster:

# To get a list of DataFrame objects, one for each cluster, containing the relevant marker features:
markers_clusters <- lapply(markerList.Cluster.markers, function(x) x)

# Assign names of genes from that to an object
geneids_clusters <- lapply(markers_clusters, function(x) x$name)

# Get out the annotations for those genes
test_list_upORdown <- lapply(geneids_clusters, function(x) AnnotationDbi::select(org.Mm.eg.db, keys=x, columns="ENTREZID", keytype="SYMBOL")$ENTREZID)

write.csv(markerList.Cluster.markers, file = "200Cluster_Test")
test_list <- read.table("200Cluster_Test",header = T,sep = ",")
head(test_list)
table(test_list$group)


min(table(test_list$group)) # Here, I will set 150 for the cluster with only 153 genes. If you want to increase it to be 200, just set max number of genes for some clusters with < 200 genes

geneid <- test_list$name

test_list$entrezid <- AnnotationDbi::select(org.Mm.eg.db, keys=geneid, columns="ENTREZID",
                                            keytype="SYMBOL")

groupprofilerlist<-list("1"=subset(test_list$entrezid$ENTREZID, test_list$group ==1)[1:200],
                         "2"=subset(test_list$entrezid$ENTREZID, test_list$group ==2)[1:200],
                          "3"=subset(test_list$entrezid$ENTREZID, test_list$group ==3)[1:200],
                          "4"=subset(test_list$entrezid$ENTREZID, test_list$group ==4)[1:200],
                         "5"=subset(test_list$entrezid$ENTREZID, test_list$group ==5)[1:200],
                         "6"=subset(test_list$entrezid$ENTREZID, test_list$group ==6)[1:200],
                          "7"=subset(test_list$entrezid$ENTREZID, test_list$group ==7)[1:200],
                          "8"=subset(test_list$entrezid$ENTREZID, test_list$group ==8)[1:200],
                         "9"=subset(test_list$entrezid$ENTREZID, test_list$group ==9)[1:200],
                         "10"=subset(test_list$entrezid$ENTREZID, test_list$group ==10)[1:200],
                          "11"=subset(test_list$entrezid$ENTREZID, test_list$group ==11)[1:200],
                          "12"=subset(test_list$entrezid$ENTREZID, test_list$group ==12)[1:200],
                         "13"=subset(test_list$entrezid$ENTREZID, test_list$group ==13)[1:200],
                         "14"=subset(test_list$entrezid$ENTREZID, test_list$group ==14)[1:200],
                          "15"=subset(test_list$entrezid$ENTREZID, test_list$group ==15)[1:200],
                          "16"=subset(test_list$entrezid$ENTREZID, test_list$group ==16)[1:200],
                         "17"=subset(test_list$entrezid$ENTREZID, test_list$group ==17)[1:200],
                         "18"=subset(test_list$entrezid$ENTREZID, test_list$group ==18)[1:200],
                          "19"=subset(test_list$entrezid$ENTREZID, test_list$group ==19)[1:200],
                          "20"=subset(test_list$entrezid$ENTREZID, test_list$group ==20)[1:200],
                         "21"=subset(test_list$entrezid$ENTREZID, test_list$group ==21)[1:200],
                         "22"=subset(test_list$entrezid$ENTREZID, test_list$group ==22)[1:200],
                          "23"=subset(test_list$entrezid$ENTREZID, test_list$group ==23)[1:200],
                          "24"=subset(test_list$entrezid$ENTREZID, test_list$group ==24)[1:200]
)

# Compare gene clusters functional profile
cc.go <- clusterProfiler::compareCluster(geneClusters = groupprofilerlist,
                                         fun = "enrichGO", OrgDb= org.Mm.eg.db, ont= "BP", pvalueCutoff=0.05, pAdjustMethod='BH')
cc.gos <- clusterProfiler::simplify(cc.go, cutoff=0.05, by= "p.adjust")

# you can set showCatergory number to choose how many terms to plot for each cluster
dotplot(cc.go, showCategory = 3)
dotplot(cc.gos, showCategory = 3)


###################################################
################################################### 
###################################################

#To obtain top 200 markers in each cluster by FDR:
#Genes are sorted by FDR by default when they're spit out from getMarkerFeatures & getMarkers

# To get a list of DataFrame objects, one for each cluster, containing the relevant marker features:
markers_clusters <- lapply(markerList.Cluster.markers, function(x) x)

# Assign names of genes from that to an object
geneids_clusters <- lapply(markers_clusters, function(x) x$name)

# Get out the annotations for those genes
test_list_upORdown <- lapply(geneids_clusters, function(x) AnnotationDbi::select(org.Mm.eg.db, keys=x, columns="ENTREZID", keytype="SYMBOL")$ENTREZID)

#write.csv(markerList.Cluster.markers, file = "200Cluster_Test")

markerList <-  markerList.Cluster.markers

mlistlens <- lapply(markerList, nrow)

#Remove any clusters with zero markers
markerList <- markerList[-which(mlistlens == 0)] 

min(unlist(mlistlens[mlistlens>0]))


filt_num <- 200 #How many genes to filter down to
top_genes <- lapply(markerList, function(x){ # for each element of markerList
  if(nrow(x) >=filt_num){ # if there are filt_num entries
    x[1:filt_num,]$name   # get the top filt_num
  }else{ # if there are not filt_num entries
    x[1:nrow(x),]$name # get all the entries
  }
})

test_list_loop <- lapply(top_genes, function(x) AnnotationDbi::select(org.Mm.eg.db, keys=x, columns="ENTREZID", keytype="SYMBOL")$ENTREZID)


sum(!groupprofilerlist$"24"[1:filt_num] %in% test_list_loop$C24)
sum(!test_list_loop$C24 %in% groupprofilerlist$"24"[1:filt_num])
# these make identical lists

# Compare gene clusters functional profile
cc.go <- clusterProfiler::compareCluster(geneClusters = test_list_loop,
                                         fun = "enrichGO", OrgDb= org.Mm.eg.db, ont= "BP", pvalueCutoff=0.05, pAdjustMethod='BH')
cc.gos <- clusterProfiler::simplify(cc.go, cutoff=0.05, by= "p.adjust")

# you can set showCatergory number to choose how many terms to plot for each cluster
dotplot(cc.go, showCategory = 3)
dotplot(cc.gos, showCategory = 3)


###### Get out single dataframe of the top 200 genes per cluster
filt_num <- 200 #How many genes to filter down to
top_genes_all_data <- lapply(markerList, function(x){ # for each element of markerList
  if(nrow(x) >=filt_num){ # if there are filt_num entries
    x[1:filt_num,]   # get the top filt_num
  }else{ # if there are not filt_num entries
    x[1:nrow(x),] # get all the entries
  }
})

top_genes_all_data


for (i in names(top_genes_all_data)){
  top_genes_all_data[[i]]["Cluster"] = i
}
topgenes_data_df <- do.call(rbind, top_genes_all_data)

write.csv(topgenes_data_df, file = "Top-Genes_Master.csv")

### Or could simply write a csv per cluster

for (i in names(top_genes_all_data)){
  write.csv(top_genes_all_data[[i]], file=paste0(i, ".csv"))
}


###################################################
################################################### 
###################################################

#To obtain top 200 markers in each cluster by log2FC:

# To get a list of DataFrame objects, one for each cluster, containing the relevant marker features:
markers_clusters <- lapply(markerList.Cluster.markers, function(x) x)

# Assign names of genes from that to an object
geneids_clusters <- lapply(markers_clusters, function(x) x$name)

# Get out the annotations for those genes
test_list_upORdown <- lapply(geneids_clusters, function(x) AnnotationDbi::select(org.Mm.eg.db, keys=x, columns="ENTREZID", keytype="SYMBOL")$ENTREZID)

#write.csv(markerList.Cluster.markers, file = "200Cluster_Test")

markerList <-  markerList.Cluster.markers

mlistlens <- lapply(markerList, nrow)

#Remove any clusters with zero markers
markerList <- markerList[-which(mlistlens == 0)] 

min(unlist(mlistlens[mlistlens>0]))


filt_num <- 200 #How many genes to filter down to
top_genes <- lapply(markerList, function(x){ # for each element of markerList
  x<-x[order(abs(x$Log2FC), decreasing = TRUE),]
  if(nrow(x) >=filt_num){ # if there are filt_num entries
    x[1:filt_num,]$name   # get the top filt_num
  }else{ # if there are not filt_num entries
    x[1:nrow(x),]$name # get all the entries
  }
})

test_list_loop <- lapply(top_genes, function(x) AnnotationDbi::select(org.Mm.eg.db, keys=x, columns="ENTREZID", keytype="SYMBOL")$ENTREZID)


# Compare gene clusters functional profile
cc.go <- clusterProfiler::compareCluster(geneClusters = test_list_loop,
                                         fun = "enrichGO", OrgDb= org.Mm.eg.db, ont= "BP", pvalueCutoff=0.05, pAdjustMethod='BH')
dotplot(cc.go, showCategory = 3)

write.table(x=cc.go@compareClusterResult, "GO_Cluster-Function.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = T,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

cc.gos <- clusterProfiler::simplify(cc.go, cutoff=0.05, by= "p.adjust")

dotplot(cc.gos, showCategory = 3)

write.table(x=cc.go@compareClusterResult, "GO-S_Cluster-Function.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = T,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


###### Get out single dataframe of the top 200 genes per cluster
filt_num <- 200 #How many genes to filter down to
top_genes_all_data <- lapply(markerList, function(x){ # for each element of markerList
  x<-x[order(abs(x$Log2FC), decreasing = TRUE),]
  if(nrow(x) >=filt_num){ # if there are filt_num entries
    x[1:filt_num,]   # get the top filt_num
  }else{ # if there are not filt_num entries
    x[1:nrow(x),] # get all the entries
  }
})

top_genes_all_data


for (i in names(top_genes_all_data)){
  top_genes_all_data[[i]]["Cluster"] = i
}
topgenes_data_df <- do.call(rbind, top_genes_all_data)

write.csv(topgenes_data_df, file = "Top-Genes_Master.csv")

### Or could simply write a csv per cluster

for (i in names(top_genes_all_data)){
  write.csv(top_genes_all_data[[i]], file=paste0(i, ".csv"))
}




