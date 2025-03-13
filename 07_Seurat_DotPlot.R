## 07_Seurat_DotPlot.R
- Code does not run. 




salloc --account=eon -t 0-09:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

conda activate seurat_2024-07

#Load libraries
R
library(Seurat)

setwd("/project/eon/nboldon/MCS2023/Seurat/")

# Load Seurat object
file_path <- "seuratObj_gsProjMCS6.Rds"
if (file.exists(file_path)) {
  seuratObj <- readRDS(file = file_path)
  print("Seurat object loaded successfully!")
  
  # Now you can access and work with the Seurat object using the variable seuratObj
  print(seuratObj)
} else {
  stop("File not found. Please check the file path and try again.")
}

################################################
################################################

# The row names in your Seurat object are currently formatted as "seqnames:start-end:name", 
# which is why the DotPlot function can't find the genes by their names alone. 
# To fix this, we need to modify the Seurat object to use only the gene names as row names. Here's how you can do that:
  
# First, extract the gene names from the current row names:
# Extract gene names from the current row names
gene_names <- sapply(strsplit(rownames(seuratObj), ":"), function(x) x[4])

# Now, set these gene names as the new row names of your Seurat object:
# Set new row names
rownames(seuratObj) <- gene_names
##Error in if (isFALSE(x = v3warn) && any(onames[[1L]] != value[[1L]])) { : 
#missing value where TRUE/FALSE needed

## The error suggests that there might be some inconsistencies in the row names or NA values. 
# Let's try a more robust approach to extract the gene names and set them as row names. 
# We'll handle potential NA values and ensure we're not introducing any issues:

# First, let's check for and handle non-unique gene names:
# Extract gene names
gene_names <- sapply(strsplit(rownames(seuratObj), ":"), function(x) x[4])

# Make gene names unique
unique_gene_names <- make.unique(gene_names)

# Check how many names were changed
sum(unique_gene_names != gene_names)

################################################
################################################

# Verify gene presence:
# Check if these genes are actually present in your dataset. You can use the following command to see all gene names:
rownames(seuratObj)

# Or to check for specific genes:
"Neurod1" %in% rownames(seuratObj)
"Olig1" %in% rownames(seuratObj)
"Olig2" %in% rownames(seuratObj)

# Check for alternative gene names:
# Sometimes genes have multiple names or symbols. Your data might use different ones. Try searching for partial matches:
grep("Neurod", rownames(seuratObj), value = TRUE)
grep("Olig", rownames(seuratObj), value = TRUE)

################################################
################################################

DotPlot(object = seuratObj, features = c("NA:NA-NA:Neurod6", "NA:NA-NA:Foxr2", "NA:NA-NA:Gm2837"))

################################################
################################################

# This example adds color customization, adjusts dot size, splits the plot by a condition (if applicable), 
# and rotates x-axis labels for better readability.

DotPlot(object = seuratObj, 
        features = c("Neurod1", "Olig1", "Olig2"),
        cols = c("lightgrey", "blue"),
        dot.scale = 8,
        split.by = "condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################
################################################

# If you want to plot genes associated with specific pathways or gene sets, 
# you can use the FeaturePlot() function with multiple features:

genes_of_interest <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
DotPlot(object = seurat_object, features = genes_of_interest)

################################################
################################################

