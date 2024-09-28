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

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

projMCS5

############################################
############################################

# Read the CSV file into a dataframe
data <- read.csv("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7/Cell-Abundance_2024-02-08.csv",
	header = TRUE, row.names = 1)

# Extract cluster labels and column names
cluster_labels <- names(data)
sample_names <- rownames(data)

# Define clusters of interest
clusters_of_interest <- c("Cluster1", "Cluster22", "Cluster3", "Cluster11","Cluster14","Cluster17","Cluster8","Cluster21","Cluster18")

# Create a list to store boxplot data
boxplot_data <- list()

# Loop through clusters of interest
for (i in seq_along(clusters_of_interest)) {
  cluster <- clusters_of_interest[i]
  # Extract data for the current cluster
  cluster_data <- subset_data[[cluster]]
  # Convert to a matrix for boxplot
  boxplot_data[[cluster]] <- as.matrix(cluster_data)
}

# Create boxplots for each cluster
for (i in seq_along(clusters_of_interest)) {
  cluster <- clusters_of_interest[i]
  boxplot(boxplot_data[[cluster]], main = cluster, xlab = "Sample", ylab = "Abundance (%)")
}

#############################################
#############################################

## Get the matrix of the heatmap to figure out row indices for the genes of interest

hmap_to_subset_MAT <- plotMarkerHeatmap(   #Anything not commented left as default
  seMarker = markerGS,   ##Set to the correct object containing the markers
  cutOff = "FDR <= 0.1 & Log2FC >= 2",
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

hmap_to_subset_MAT

Quantile(hmap_to_subset_MAT$name

##########################################
##########################################

## Cluster 18 Cell Abundance by Tx Group

groupT1 <- c(0.4344, 0.4153, 0.4155, 0.4451, 0.4317, 0.4513, 0.4591, 0.4201, 0.4613, 0.4133, 0.3920, 0.4504)

groupT2 <- c(0.4193, 0.3690, 0.4167, 0.4080, 0.4227, 0.4404, 0.4145, 0.4440, 0.4150, 0.4632, 0.4485, 0.4273, 0.3957, $

groupT3 <- c(0.3763, 0.4200, 0.4128, 0.4101, 0.4161, 0.4207, 0.3771, 0.4227, 0.4504, 0.4254, 0.4088, 0.3946, 0.4229, $

groupT4 <- c(0.4308, 0.3751, 0.4799, 0.4093, 0.4389, 0.4340, 0.4060, 0.3975, 0.4492, 0.4086, 0.4881, 0.4197)


# Perform ANOVA
result <- aov(c(groupT1,  groupT4) ~ rep(c("groupT1", "groupT4")))

# Summarize ANOVA results
summary(result)


##########################################
##########################################

# Read in the data
data <- read.csv("/Volumes/DataBox/MCS2023/Marker_Analysis_2024-Feb/Cell-Abundance_2024-02-08.csv")

# Create groups of data

groupT1 <- c("C302", "C306", "C309", "C318", "C328", "C332", "C337", "C339", "C346", "C351", "C353", "C360")

groupT2 <- c("C304", "C308", "C312", "C315", "C321", "C324", "C327", "C330", "C333", "C336", "C342", "C348", "C349", $

groupT3 <- c("C305", "C307", "C313", "C316", "C320", "C322", "C323", "C325", "C334", "C340", "C341", "C345", "C350", $

groupT4 <- c("C301", "C303", "C310", "C314", "C319", "C335", "C338", "C344", "C354", "C356", "C361", "C363")


# Perform ANOVA
result <- aov(c(groupT1, groupT2, groupT3, groupT4) ~ rep(c("groupT1", "groupT2", "groupT3", "groupT4"), each = 25))

# Summarize ANOVA results
summary(result)



############################################
###########################################

# Boxplot


