library(Seurat)
library(Signac)
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
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

projMCS7
getAvailableMatrices(projMCS7)

###################################################
###################################################
###################################################
###################################################
###################################################
###################################################



## Gene markers identified by cell type here: /Volumes/DataBox/MCS2023/Tx_Comp/

# getMarkerFeatures used groupBy = cellTypes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25
# see projMCS7_GeneMarkers_byCellType_2025-01-16.R for additional code (duplicated below)


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
cellType_markerList <- getMarkers(markersGenes, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")

# Get the names of all cell types
celltype_names <- names(cellType_markerList)

# Loop through each cellType and save its marker list to a separate CSV file
for (cellType in celltype_names) {
  # Create a file name based on the cellType
  outputFileName <- paste0(cellType, "_GeneMarker_List_", Sys.Date(), ".csv")
  
  # Save the marker list for the current cellType to a CSV file
  write.csv(markerList[[cellType]], file = outputFileName)
  
  # Print message to track progress
  print(paste0("Saved marker list for cellType: ", cellType))
}



###################################################
###################################################
###################################################


#To obtain top GO terms in each cellType:

# To get a list of DataFrame objects, one for each cluster, containing the relevant marker features:
markers_cellTypes <- lapply(cellType_markerList, function(x) x)

# Assign names of genes from that to an object
geneids_cellTypes <- lapply(cellType_markerList, function(x) x$name)

# Get out the annotations for those genes
geneids_annotated <- lapply(geneids_cellTypes, function(x) {
  AnnotationDbi::select(org.Mm.eg.db, keys = x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
})

groupprofilerlist <- list(
  "glut" = geneids_annotated$glut,
  "glutPrecursor" = geneids_annotated$glutPrecursor,
  "gaba" = geneids_annotated$gaba,
  "microglia" = geneids_annotated$microglia,
  "oligo" = geneids_annotated$oligo,
  "oligoPrecursor" = geneids_annotated$oligoPrecursor,
  "glutOligo" = geneids_annotated$glutOligo,
  "endo" = geneids_annotated$endo,
  "astro" = geneids_annotated$astro,
  "astroPrecursor" = geneids_annotated$astroPrecursor,
  "glutAstro" = geneids_annotated$glutAstro
)

write.csv(cellType_markerList, file = "geneMarkers_byCelltype_2025-01-18.csv")
test_list <- read.table("top200_byCelltype_2025-01-18.csv",header = T,sep = ",")
head(test_list)
table(test_list$group)

cc.go <- clusterProfiler::compareCluster(
  geneClusters = groupprofilerlist,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP", # Biological Process
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)


cc.gos <- clusterProfiler::simplify(
  cc.go, cutoff=0.05, 
  by= "p.adjust")


# you can set showCatergory number to choose how many terms to plot for each cluster
dotplot(cc.go, showCategory = 10) + ggtitle("GO Enrichment Across Cell Types")
dotplot(cc.gos, showCategory = 10) + ggtitle("GO Simplified Enrichment Across Cell Types")



# Create plots and assign them to variables
dotplot_original <- dotplot(cc.go, showCategory = 10) + 
  ggtitle("GO Enrichment Across Cell Types")
dotplot_simplified <- dotplot(cc.gos, showCategory = 10) + 
  ggtitle("GO Simplified Enrichment Across Cell Types")

# Save the plots as image files (e.g., PNG or PDF)
ggsave("GO_Enrichment_Plot_2025-01-18.png", dotplot_original, width = 13, height = 25)
ggsave("GO_Simplified_Enrichment_Plot_2025-01-18.png", dotplot_simplified, width = 13, height =15)

# Optional: Save as PDF for high resolution
ggsave("GO_Enrichment_Plot_2025-01-18.pdf", dotplot_original, width = 13, height =25)
ggsave("GO_Simplified_Enrichment_Plot_2025-01-18.pdf", dotplot_simplified, width = 13, height = 15)



# Convert results to data frames
cc_go_df <- as.data.frame(cc.go)
cc_gos_df <- as.data.frame(cc.gos)

# Save as CSV files
write.csv(cc_go_df, "GO_Enrichment_Results_2025-01-18.csv", row.names = FALSE)
write.csv(cc_gos_df, "GO_Simplified_Enrichment_Results_2025-01-18.csv", row.names = FALSE)

###################################################
################################################### 
###################################################
###################################################
################################################### 
###################################################



## To obtain top 200 markers in each cluster by log2FC:
# Continuing from above code
# Sorts by the magnitude of abs(Log2FC) and does not differentiate between positive or negative


# Filter top 200 genes by abs(Log2FC)
top200_markers <- lapply(cellType_markerList, function(df) {
  df <- df[order(abs(df$Log2FC), decreasing = TRUE), ] # Sort by abs(Log2FC)
  head(df, 200) # Keep the top 200 genes
})

geneids_top200 <- lapply(top200_markers, function(x) x$name)

geneids_annotated_top200 <- lapply(geneids_top200, function(x) {
  AnnotationDbi::select(org.Mm.eg.db, keys = x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID
})

groupprofilerlist_top200 <- list(
  "glut" = geneids_annotated_top200$glut,
  "glutPrecursor" = geneids_annotated_top200$glutPrecursor,
  "gaba" = geneids_annotated_top200$gaba,
  "microglia" = geneids_annotated_top200$microglia,
  "oligo" = geneids_annotated_top200$oligo,
  "oligoPrecursor" = geneids_annotated_top200$oligoPrecursor,
  "glutOligo" = geneids_annotated_top200$glutOligo,
  "endo" = geneids_annotated_top200$endo,
  "astro" = geneids_annotated_top200$astro,
  "astroPrecursor" = geneids_annotated_top200$astroPrecursor,
  "glutAstro" = geneids_annotated_top200$glutAstro
)

cc.go_top200 <- clusterProfiler::compareCluster(
  geneClusters = groupprofilerlist_top200,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP", # Biological Process
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)


cc.gos_top200 <- clusterProfiler::simplify(
  cc.go_top200, cutoff=0.05, 
  by= "p.adjust")


# you can set showCatergory number to choose how many terms to plot for each cluster
dotplot(cc.go_top200, showCategory = 10) + ggtitle("GO Top 200 abs(Log2FC) Enrichment Across Cell Types")
dotplot(cc.gos_top200, showCategory = 10) + ggtitle("GO Top 200 abs(Log2FC) Simplified Enrichment Across Cell Types")



# Create plots and assign them to variables
dotplot_top200_original <- dotplot(cc.go_top200, showCategory = 10) + 
  ggtitle("GO Top 200 abs(Log2FC) Across Cell Types")
dotplot_top200_simplified <- dotplot(cc.gos_top200, showCategory = 10) + 
  ggtitle("GO Top 200 abs(Log2FC) Simplified Enrichment Across Cell Types")

# Save the plots as image files (e.g., PNG or PDF)
ggsave("GO_top200_Plot_2025-01-18.png", dotplot_top200_original, width = 13, height = 25)
ggsave("GO_top200_Simplified_Plot_2025-01-18.png", dotplot_top200_simplified, width = 13, height =15)

# Optional: Save as PDF for high resolution
ggsave("GO_top200_Plot_2025-01-18.pdf", dotplot_top200_original, width = 13, height =25)
ggsave("GO_top200_Simplified_Plot_2025-01-18.pdf", dotplot_top200_simplified, width = 13, height = 15)



# Convert results to data frames
cc_go_top200_df <- as.data.frame(cc.go_top200)
cc_gos_top200_df <- as.data.frame(cc.gos_top200)

# Save as CSV files
write.csv(cc_go_top200_df, "GO_top200_Results_2025-01-18.csv", row.names = FALSE)
write.csv(cc_gos_top200_df, "GO_top200_Simplified_Results_2025-01-18.csv", row.names = FALSE)



