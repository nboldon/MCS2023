## Comparing motifs after subsetting by treatment group

#Load libraries
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
library(RColorBrewer)

#Additional setup
setwd("/Volumes/DataBox/ProjMCS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)


# Load the ArchR project
projMCS7 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS7", force = FALSE, showLogo = FALSE)

# Check the project details and available matrices
projMCS7
getAvailableMatrices(projMCS7)


########################################
########################################
########################################


# Set up treatment groups by renaming sample IDs
treatment <- projMCS7$Sample

# Rename treatment groups based on sample IDs
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

# Assign the modified treatment names back to the project
projMCS7$treatment <- treatment


########################################
########################################

## Subset treatment groups of interest

t1Subset <- projMCS7[
  projMCS7$treatment %in% c("t1"),]
t1Subset

# numberOfCells: 24969
# medianTSS: 14.968
# medianFrags: 7696

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t1Subset),
  outputDirectory = "t1Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t1Subset)


############################################

t2Subset <- projMCS7[
  projMCS7$treatment %in% c("t2"),]
t2Subset

# numberOfCells: 30149
# medianTSS: 14.749
# medianFrags: 6860

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t2Subset),
  outputDirectory = "t2Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t2Subset)

############################################

t3Subset <- projMCS7[
  projMCS7$treatment %in% c("t3"),]
t3Subset

# numberOfCells: 27879
# medianTSS: 14.332
# medianFrags: 6533

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t3Subset),
  outputDirectory = "t3Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t3Subset)

############################################

t4Subset <- projMCS7[
  projMCS7$treatment %in% c("t4"),]
t4Subset

# numberOfCells: 21458
# medianTSS: 13.96
# medianFrags: 5522

subsetArchRProject(
  ArchRProj = projMCS7,
  cells = getCellNames(t4Subset),
  outputDirectory = "t4Subset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

getAvailableMatrices(t4Subset)

############################################
############################################
############################################
############################################


## Repeat the following for each treatment group


t4MarkersPeaks <- getMarkerFeatures(
  ArchRProj = t4Subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

t4MarkersPeaks

# Motif enrichment in marker peaks

t4EnrichMotifs <- peakAnnoEnrichment(
  seMarker = t4MarkersPeaks,
  ArchRProj = t4Subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
)

t4EnrichMotifs

t4HeatmapEM <- plotEnrichHeatmap(t4EnrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(t4HeatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(t4HeatmapEM, name = "t4-Motifs-Enriched-Marker-Heatmap_2025-01-11", width = 8, height = 6, ArchRProj = t4Subset, addDOC = FALSE)

#############
#############
#############

## Unable to save subset projects

# Error in saveArchRProject(ArchRProj = t1Subset, outputDirectory = "/Volumes/DataBox/Save-t1Subset",  : 
# all(file.exists(zfiles)) is not TRUE

saveArchRProject(ArchRProj = t1Subset, outputDirectory = "/Volumes/DataBox/Save-t1Subset", load = FALSE)

t1Subset <- loadArchRProject(path = "/Volumes/DataBox/Save-t1Subset", force = FALSE, showLogo = FALSE)

############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################


## Motif pairwise comparisons by cell type


############################################


# Set up treatment groups by renaming sample IDs
treatment <- projMCS7$Sample

# Rename treatment groups based on sample IDs
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

# Assign the modified treatment names back to the project
projMCS7$treatment <- treatment


############################################


# Define a named vector to map clusters to cell types

cluster_to_cell_type <- c(
  "C1" = "AstroPrecursor",
  "C2" = "Oligodendrocyte",
  "C3" = "Oligodendrocyte",
  "C4" = "GlutOligo",
  "C5" = "OligoPrecursor",
  "C6" = "OligoPrecursor",
  "C7" = "AstroPrecursor",
  "C8" = "Astrocyte",
  "C9" = "AstroPrecursor",
  "C10" = "Microglia",
  "C11" = "Microglia",
  "C12" = "Endo-Vasc",
  "C13" = "Endo-Vasc",
  "C14" = "Endo-Vasc",
  "C15" = "GlutPrecursor",
  "C16" = "GlutPrecursor",
  "C17" = "GlutPrecursor",
  "C18" = "Glutamatergic",
  "C19" = "Glutamatergic",
  "C20" = "GlutPrecursor",
  "C21" = "Glutamatergic",
  "C22" = "Gabaergic",
  "C23" = "Gabaergic",
  "C24" = "GlutAstro",
  "C25" = "GlutPrecursor"
)

# Add cell type information to projMCS7 metadata
projMCS7$cell_type <- cluster_to_cell_type[projMCS7$Clusters]

# Check the distribution of cell types
table(projMCS7$cell_type)

# Optionally, set cell type as a factor for consistent ordering
#projMCS7$cell_type <- factor(
#  projMCS7$cell_type,
#  levels = unique(cluster_to_cell_type)
#)


############################################
############################################
############################################


# Check the unique values in the 'cell_type' and 'treatment' columns
unique(projMCS7$cell_type)
unique(projMCS7$treatment)

# Create combined group labels like "Glutamatergic_t1"
projMCS7$cell_treatment <- paste(projMCS7$cell_type, projMCS7$treatment, sep = "_")

# Check the unique combined labels
unique(projMCS7$cell_treatment)


# Check the New Groups Inspect the new grouping variable to ensure it is properly created:
table(projMCS7$cell_treatment)

# Astrocyte_t1       Astrocyte_t2       Astrocyte_t3       Astrocyte_t4  AstroPrecursor_t1 
# 1982               2530               1985               1523                 87 
# AstroPrecursor_t2  AstroPrecursor_t3  AstroPrecursor_t4       Endo-Vasc_t1       Endo-Vasc_t2 
# 122                114                 81                842               1075 
# Endo-Vasc_t3       Endo-Vasc_t4       Gabaergic_t1       Gabaergic_t2       Gabaergic_t3 
# 879                621               2457               2882               2702 
# Gabaergic_t4   Glutamatergic_t1   Glutamatergic_t2   Glutamatergic_t3   Glutamatergic_t4 
# 2222              13381              16079              14894              11689 
# GlutAstro_t1       GlutAstro_t2       GlutAstro_t3       GlutAstro_t4       GlutOligo_t1 
# 176                248                253                197                239 
# GlutOligo_t2       GlutOligo_t3       GlutOligo_t4   GlutPrecursor_t1   GlutPrecursor_t2 
# 254                258                254               1821               2203 
# GlutPrecursor_t3   GlutPrecursor_t4       Microglia_t1       Microglia_t2       Microglia_t3 
# 2292               1659                984               1171               1069 
# Microglia_t4 Oligodendrocyte_t1 Oligodendrocyte_t2 Oligodendrocyte_t3 Oligodendrocyte_t4 
# 659               2469               2996               2893               2135 
# OligoPrecursor_t1  OligoPrecursor_t2  OligoPrecursor_t3  OligoPrecursor_t4 
# 531                589                540                418 

# Double check assignments:
table(projMCS7$cell_type, projMCS7$treatment)

# t1    t2    t3    t4
# Astrocyte        1982  2530  1985  1523
# AstroPrecursor     87   122   114    81
# Endo-Vasc         842  1075   879   621
# Gabaergic        2457  2882  2702  2222
# Glutamatergic   13381 16079 14894 11689
# GlutAstro         176   248   253   197
# GlutOligo         239   254   258   254
# GlutPrecursor    1821  2203  2292  1659
# Microglia         984  1171  1069   659
# Oligodendrocyte  2469  2996  2893  2135
# OligoPrecursor    531   589   540   418

table(is.na(projMCS7$cell_type), is.na(projMCS7$treatment))
# FALSE
# FALSE 104455

# Use in Analysis Replace Clusters or cell_type with cell_treatment in your grouping arguments during analysis
# For example: When summarizing or visualizing fragment counts
fragmentSummary <- tapply(projMCS7$nFrags, projMCS7$cell_treatment, summary)

# When performing differential analysis
markers <- getMarkerFeatures(
  ArchRProj = projMCS7, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_treatment"
)


############################################
############################################
############################################

## Barplot for cell types and treatment


# Load the viridis package
library(viridis)
# Generate a barplot 
# Save the plot to a PDF with a larger size and adjusted margins
pdf("Cell-Cts_byCellType-Tx_2025-01-13.pdf", width = 12, height = 15)
# Adjust the margins further
par(mar = c(20, 5, 4, 2))  # Increase bottom margin (first value in `mar`)
# Generate the barplot with custom ylim and rotated labels
barplot(table(projMCS7$cell_treatment), 
        las = 2,  # Keep the labels rotated vertically
        col = viridis(length(table(projMCS7$cell_treatment))), 
        main = "Cell Counts by Cell Type and Treatment",
        #xlab = "Cell Type and Treatment",
        ylab = "Cell Count",
        border = NA,  # Optional to remove borders for a cleaner look
        cex.axis = 0.8,  # Smaller axis text size
        ylim = c(0, max(table(projMCS7$cell_treatment)) * 1.3))  # Adjust y-axis to fit the tallest bars
dev.off()


############################################
############################################
############################################
############################################
############################################
############################################


## Boxplot for cell type and treatment


# Load required libraries
library(dplyr)
library(tidyr)

# Extract the cellColData as a data frame
cell_data <- as.data.frame(projMCS7@cellColData)

# Summarize data: count the number of cells per cluster for each sample
summary_data <- cell_data %>%
  group_by(Sample, Clusters) %>%
  summarise(CellCount = n(), .groups = "drop")

# Convert to wide format: each sample as a column
wide_data <- summary_data %>%
  pivot_wider(names_from = Sample, values_from = CellCount, values_fill = 0)

# View the resulting dataframe
print("Wide Data Before Adding Treatment Groups:")
print(wide_data)



# Create a new row to map sample IDs to treatments

treatment_mapping <- c(
  "C302_" = "t1", "C306_" = "t1", "C309_" = "t1", "C318_" = "t1",
  "C323_" = "t1", "C328_" = "t1", "C332_" = "t1", "C337_" = "t1",
  "C339_" = "t1", "C346_" = "t1", "C351_" = "t1", "C353_" = "t1",
  "C360_" = "t1", "C304_" = "t2", "C308_" = "t2", "C312_" = "t2",
  "C349_" = "t2", "C315_" = "t2", "C321_" = "t2", "C324_" = "t2",
  "C355_" = "t2", "C327_" = "t2", "C330_" = "t2", "C333_" = "t2",
  "C358_" = "t2", "C336_" = "t2", "C342_" = "t2", "C348_" = "t2",
  "C362_" = "t2", "C305_" = "t3", "C307_" = "t3", "C313_" = "t3",
  "C350_" = "t3", "C316_" = "t3", "C320_" = "t3", "C322_" = "t3",
  "C352_" = "t3", "C325_" = "t3", "C334_" = "t3", "C359_" = "t3",
  "C340_" = "t3", "C341_" = "t3", "C345_" = "t3", "C364_" = "t3",
  "C301_" = "t4", "C303_" = "t4", "C310_" = "t4", "C314_" = "t4",
  "C319_" = "t4", "C335_" = "t4", "C338_" = "t4", "C344_" = "t4",
  "C354_" = "t4", "C356_" = "t4", "C361_" = "t4", "C363_" = "t4"
)

# Extract the sample names from the wide_data columns (ignoring the first column, which contains clusters)
sample_names <- colnames(wide_data)[-1]

# Map sample names to treatments
treatment_row <- sapply(sample_names, function(sample) treatment_mapping[[sample]])

# Check if treatment_row is populated correctly
print(treatment_row)

# Convert all columns (except Clusters) to character type
wide_data <- wide_data %>%
  mutate(across(-Clusters, as.character))

# Add the treatment row
wide_data <- wide_data %>%
  add_row(Clusters = "TreatmentGroup", !!!as.list(treatment_row))

# View the updated wide_data
print("Wide Data After Adding Treatment Groups:")
print(wide_data)

treatment_row_check <- wide_data %>% filter(Clusters == "TreatmentGroup")
print(treatment_row_check)

tail(wide_data)



# Create a new column for cell types
cluster_to_cell_type <- c(
  "C1" = "AstroPrecursor",
  "C2" = "Oligodendrocyte",
  "C3" = "Oligodendrocyte",
  "C4" = "GlutOligo",
  "C5" = "OligoPrecursor",
  "C6" = "OligoPrecursor",
  "C7" = "AstroPrecursor",
  "C8" = "Astrocyte",
  "C9" = "AstroPrecursor",
  "C10" = "Microglia",
  "C11" = "Microglia",
  "C12" = "Endo-Vasc",
  "C13" = "Endo-Vasc",
  "C14" = "Endo-Vasc",
  "C15" = "GlutPrecursor",
  "C16" = "GlutPrecursor",
  "C17" = "GlutPrecursor",
  "C18" = "Glutamatergic",
  "C19" = "Glutamatergic",
  "C20" = "GlutPrecursor",
  "C21" = "Glutamatergic",
  "C22" = "Gabaergic",
  "C23" = "Gabaergic",
  "C24" = "GlutAstro",
  "C25" = "GlutPrecursor"
)

# Check if cluster names in wide_data match the keys in the mapping
print(unique(wide_data$Clusters))

# Assign cell types based on cluster
wide_data$cell_type <- cluster_to_cell_type[as.character(wide_data$Clusters)]

# View the updated table
head(wide_data)



# Save the wide_data as a .csv file
write.csv(wide_data, "cellCount_byTx-CellType_2025-01-14.csv", row.names = FALSE)


#############################

## Create boxplots for individual .csv file


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Read in the data
data <- read.csv("microglia.csv", header = TRUE)

# Check the column names to see the issue clearly
colnames(data)

# Reshape the data to long format, treating all columns as 'Treatment'
data_long <- data %>%
  gather(key = "Treatment", value = "CellCount") %>%
  # Extract the main treatment name (before the .1, .2, etc.) using sub() function
  mutate(Treatment = sub("\\..*", "", Treatment)) %>%
  mutate(Treatment = factor(Treatment, levels = c("t1", "t2", "t3", "t4")))

# Check the reshaped data
head(data_long)

# Create the boxplot with individual points, mean line, and the viridis color palette
ggplot(data_long, aes(x = Treatment, y = CellCount, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(color = "black", width = 0.2, size = 2) +  # Add individual points as dots
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "red") +  # Add mean as red points
  scale_fill_viridis(discrete = TRUE) +  # Apply viridis color palette
  labs(x = "Microglia", y = "Total Cell Counts") +  # Axis labels
  theme_minimal() +  # Clean theme
  theme(legend.position = "none")  # Remove the legend


#############################

## Create boxplots for all .csv files side by side


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)

# Define a function to create the boxplot for each dataset
create_boxplot <- function(file_name, title) {
  # Read in the data
  data <- read.csv(file_name, header = TRUE)
  
  # Reshape the data to long format and clean up treatment names
  data_long <- data %>%
    gather(key = "Treatment", value = "CellCount") %>%
    # Remove suffixes (e.g., .1, .2) from treatment names
    mutate(Treatment = sub("\\..*", "", Treatment)) %>%
    mutate(Treatment = factor(Treatment, levels = c("t1", "t2", "t3", "t4")))
  
  # Create the boxplot
  plot <- ggplot(data_long, aes(x = Treatment, y = CellCount, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(color = "black", width = 0.2, size = 2) +  # Add individual points as dots
    stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "red") +  # Add mean as red points
    scale_fill_viridis(discrete = TRUE) +  # Apply viridis color palette
    labs(x = title) +  # Set the title for the x-axis (dataset name)
    theme_minimal() +  # Clean theme
    theme(legend.position = "none", axis.title.y = element_blank())  # Remove y-axis label for all except first plot
  
  return(plot)
}

# Create boxplots for all the datasets
files <- c("Glut.csv", "Glut_Precursor.csv", "Glut_Astro.csv", "Glut_Oligo.csv", "GABA.csv", 
           "Oligo.csv", "Oligo_Precursor.csv", "Microglia.csv", "Astrocytes.csv", "Astro_Precursor.csv", 
           "Endo-Vasc.csv")

# Create plots for each dataset and store them in a list
plots <- lapply(seq_along(files), function(i) {
  title <- gsub(".csv", "", files[i])  # Extract the title from the filename (remove .csv)
  plot <- create_boxplot(files[i], title)
  
  # Add the y-axis label only for the first plot
  if (i == 1) {
    plot <- plot + labs(y = "Total Cell Counts")  # Add y-axis label for the first plot
  }
  
  return(plot)
})

# Save the plots to a PDF file
pdf("CellCount_byTx-CellType_2025-01-15.pdf", width = 16, height = 12)  # Adjust size if needed
grid.arrange(grobs = plots, ncol = 3)  # Arrange the plots in a grid
dev.off()  # Close the PDF device


################################


## To display boxplots all in a line: 

# Define a function to create the boxplot for each dataset
create_boxplot <- function(file_name, title, add_y_label = FALSE) {
  # Read in the data
  data <- read.csv(file_name, header = TRUE)
  
  # Reshape the data to long format and clean up treatment names
  data_long <- data %>%
    gather(key = "Treatment", value = "CellCount") %>%
    # Remove suffixes (e.g., .1, .2) from treatment names
    mutate(Treatment = sub("\\..*", "", Treatment)) %>%
    mutate(Treatment = factor(Treatment, levels = c("t1", "t2", "t3", "t4")))
  
  # Create the boxplot
  plot <- ggplot(data_long, aes(x = Treatment, y = CellCount, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(color = "black", width = 0.2, size = 2) +  # Add individual points as dots
    stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "red") +  # Add mean as red points
    scale_fill_viridis(discrete = TRUE) +  # Apply viridis color palette
    labs(x = title) +  # Set the title for the x-axis (dataset name)
    theme_minimal() +  # Clean theme
    theme(legend.position = "none", axis.title.y = element_blank())  # Remove y-axis label for all except first plot
  
  # Add y-axis label only if specified
  if (add_y_label) {
    plot <- plot + labs(y = "Total Cell Counts")  # Add y-axis label for the first plot
  }
  
  return(plot)
}

# Create boxplots for all the datasets
files <- c("Glut.csv", "Glut_Precursor.csv", "Glut_Astro.csv", "Glut_Oligo.csv", "GABA.csv", 
           "Oligo.csv", "Oligo_Precursor.csv", "Microglia.csv", "Astrocytes.csv", "Astro_Precursor.csv", 
           "Endo-Vasc.csv")

# Create plots for each dataset and store them in a list
plots <- lapply(seq_along(files), function(i) {
  title <- gsub(".csv", "", files[i])  # Extract the title from the filename (remove .csv)
  
  # Add y-axis label only to the first plot
  add_y_label <- ifelse(i == 1, TRUE, FALSE)
  
  plot <- create_boxplot(files[i], title, add_y_label)
  return(plot)
})

# Save the plots to a PDF file
pdf("CellCount_byTx-CellType_2025-01-15.pdf", width = 16, height = 4)  # Adjust height for a continuous line
grid.arrange(grobs = plots, ncol = length(files))  # Arrange the plots in a single row
dev.off()  # Close the PDF device



#############################


## To include Anova tests to test for significant differences between treatment groups


create_boxplot <- function(file_name, title, add_y_label = FALSE) {
  # Read in the data
  data <- read.csv(file_name, header = TRUE)
  
  # Reshape the data to long format and clean up treatment names
  data_long <- data %>%
    gather(key = "Treatment", value = "CellCount") %>%
    mutate(Treatment = sub("\\..*", "", Treatment)) %>%
    mutate(Treatment = factor(Treatment, levels = c("t1", "t2", "t3", "t4")))
  
  # Perform ANOVA to test for differences between treatment groups
  anova_result <- aov(CellCount ~ Treatment, data = data_long)
  
  # Extract the p-value from ANOVA (proper indexing to avoid NA)
  anova_summary <- summary(anova_result)
  anova_p_value <- anova_summary[[1]]$`Pr(>F)`[1]  # Correctly extract the p-value for Treatment
  
  # Create the boxplot
  plot <- ggplot(data_long, aes(x = Treatment, y = CellCount, fill = Treatment)) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(color = "black", width = 0.2, size = 2) +  # Add individual points as dots
    stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "red") +  # Add mean as red points
    scale_fill_viridis(discrete = TRUE) +  # Apply viridis color palette
    labs(x = title) +  # Set the title for the x-axis (dataset name)
    theme_minimal() +  # Clean theme
    theme(legend.position = "none", axis.title.y = element_blank())  # Remove y-axis label for all except first plot
  
  # Add y-axis label only if specified
  if (add_y_label) {
    plot <- plot + labs(y = "Total Cell Counts")  # Add y-axis label for the first plot
  }
  
  # Add ANOVA p-value to the plot (without additional information)
  y_max <- max(data_long$CellCount) + 10  # Adjust y position for text
  plot <- plot +
    geom_text(aes(x = 2, y = y_max, label = paste("p =", round(anova_p_value, 3))),
              color = "black", size = 3, vjust = 0)  # Display only p-value
  
  return(plot)
}

# Create boxplots for all the datasets
files <- c("Glut.csv", "Glut_Precursor.csv", "Glut_Astro.csv", "Glut_Oligo.csv", "GABA.csv", 
           "Oligo.csv", "Oligo_Precursor.csv", "Microglia.csv", "Astrocytes.csv", "Astro_Precursor.csv", 
           "Endo-Vasc.csv")

# Create plots for each dataset and store them in a list
plots <- lapply(seq_along(files), function(i) {
  title <- gsub(".csv", "", files[i])  # Extract the title from the filename (remove .csv)
  
  # Add y-axis label only to the first plot
  add_y_label <- ifelse(i == 1, TRUE, FALSE)
  
  plot <- create_boxplot(files[i], title, add_y_label)
  return(plot)
})

# Save the plots to a PDF file
pdf("CellCount_byTx-CellType_2025-01-15_with_ANOVA_p_values_only.pdf", width = 16, height = 4)  # Adjust height for a continuous line
grid.arrange(grobs = plots, ncol = length(files))  # Arrange the plots in a single row
dev.off()  # Close the PDF device





############################################
############################################
############################################
############################################
############################################
############################################


## Comparing t1 to t3 Glutamatergic cells


t1_t3_glut_markers <- getMarkerFeatures(
  ArchRProj = projMCS7,
  useMatrix = "PeakMatrix",
  groupBy = "cell_treatment",  # Use the combined cell type and treatment groups
  useGroups = "Glutamatergic_t1",  # Specify the group of interest
  bgdGroups = "Glutamatergic_t3",  # Specify the background group
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


t1_t3_glut_EnrichMotifs <- peakAnnoEnrichment(
  seMarker = t1_t3_glut_markers,
  ArchRProj = projMCS7,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
)

t1_t3_glut_EnrichMotifs

t1_t3_glut_HeatmapEM <- plotEnrichHeatmap(t1_t3_glut_EnrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(t1_t3_glut_HeatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(t1_t3_glut_HeatmapEM, name = "t1_t3_glut_-Motifs-Enriched-Marker-Heatmap_2025-01-13", width = 8, height = 6, ArchRProj = projMCS7, addDOC = FALSE)



############################################
############################################
############################################


## Loop to run through all cell types and treatment group comparisons


# Extract unique cell types and treatments
cell_types <- unique(projMCS7$cell_type)
treatments <- unique(projMCS7$treatment)

# Create an empty list to store results
results <- list()

# Loop through all combinations of cell types and treatments
for (cell_type in cell_types) {
  for (treatment_1 in treatments) {
    for (treatment_2 in treatments) {
      
      # Ensure you're comparing different treatments
      if (treatment_1 != treatment_2) {
        
        # Create the combined label for each comparison
        useGroup <- paste(cell_type, treatment_1, sep = "_")
        bgdGroup <- paste(cell_type, treatment_2, sep = "_")
        
        # Run the marker analysis within a tryCatch block to handle errors
        tryCatch({
          # Run marker features analysis
          markers <- getMarkerFeatures(
            ArchRProj = projMCS7,
            useMatrix = "PeakMatrix",
            groupBy = "cell_treatment",
            useGroups = useGroup,
            bgdGroups = bgdGroup,
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "wilcoxon"
          )
          
          # Perform enrichment analysis
          enrichMotifs <- peakAnnoEnrichment(
            seMarker = markers,
            ArchRProj = projMCS7,
            peakAnnotation = "Motif",
            cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5"
          )
          
          # Check if enrichments exist
          if (nrow(enrichMotifs) == 0) {
            message(paste("No enrichments found for", useGroup, "vs", bgdGroup))
          } else {
            # Store the enrichment results
            results[[paste(useGroup, bgdGroup, sep = "_")]] <- enrichMotifs
            # Plot the heatmap
            enrichHeatmap <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
            ComplexHeatmap::draw(enrichHeatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
            plotPDF(enrichHeatmap, 
                    name = paste(useGroup, bgdGroup, "Motifs-Enriched-Marker-Heatmap", Sys.Date(), sep = "_"),
                    width = 8, height = 6, ArchRProj = projMCS7, addDOC = FALSE)
          }
          
        }, error = function(e) {
          message(paste("Error with", useGroup, "vs", bgdGroup, ":", e))
        })
      }
    }
  }
}


## No enrichments found for the following comparisons:
# GlutPrecursor_t4 vs GlutPrecursor_t1
# GlutPrecursor_t4 vs GlutPrecursor_t2
# GlutPrecursor_t4 vs GlutPrecursor_t3
# GlutPrecursor_t1 vs GlutPrecursor_t4
# GlutPrecursor_t1 vs GlutPrecursor_t2
# GlutPrecursor_t1 vs GlutPrecursor_t3
# GlutPrecursor_t2 vs GlutPrecursor_t4
# GlutPrecursor_t2 vs GlutPrecursor_t1
# GlutPrecursor_t2 vs GlutPrecursor_t3
# GlutPrecursor_t3 vs GlutPrecursor_t4
# GlutPrecursor_t3 vs GlutPrecursor_t1
# GlutPrecursor_t3 vs GlutPrecursor_t2
# Glutamatergic_t4 vs Glutamatergic_t1
# Glutamatergic_t4 vs Glutamatergic_t2
# Glutamatergic_t4 vs Glutamatergic_t3
# Glutamatergic_t1 vs Glutamatergic_t4
# Glutamatergic_t1 vs Glutamatergic_t2
# Glutamatergic_t1 vs Glutamatergic_t3
# Glutamatergic_t2 vs Glutamatergic_t4
# Glutamatergic_t2 vs Glutamatergic_t1
# Glutamatergic_t2 vs Glutamatergic_t3
# Glutamatergic_t3 vs Glutamatergic_t4
# Glutamatergic_t3 vs Glutamatergic_t1
# Glutamatergic_t3 vs Glutamatergic_t2
# GlutOligo_t4 vs GlutOligo_t1
# GlutOligo_t4 vs GlutOligo_t2 
# GlutOligo_t4 vs GlutOligo_t3
# GlutOligo_t1 vs GlutOligo_t4
# GlutOligo_t1 vs GlutOligo_t2
# GlutOligo_t1 vs GlutOligo_t3
# GlutOligo_t2 vs GlutOligo_t4
# GlutOligo_t2 vs GlutOligo_t1
# GlutOligo_t2 vs GlutOligo_t3
# GlutOligo_t3 vs GlutOligo_t4
# GlutOligo_t3 vs GlutOligo_t1
# GlutOligo_t3 vs GlutOligo_t2
# Endo-Vasc_t4 vs Endo-Vasc_t1
# Endo-Vasc_t4 vs Endo-Vasc_t2
# Endo-Vasc_t4 vs Endo-Vasc_t3
# Endo-Vasc_t1 vs Endo-Vasc_t4
# Endo-Vasc_t1 vs Endo-Vasc_t2
# Endo-Vasc_t1 vs Endo-Vasc_t3
# Endo-Vasc_t2 vs Endo-Vasc_t4
# Endo-Vasc_t2 vs Endo-Vasc_t1
# Endo-Vasc_t2 vs Endo-Vasc_t3
# Endo-Vasc_t3 vs Endo-Vasc_t4 
# Endo-Vasc_t3 vs Endo-Vasc_t1
# Endo-Vasc_t3 vs Endo-Vasc_t2
# Gabaergic_t4 vs Gabaergic_t1
# Gabaergic_t4 vs Gabaergic_t2
# Gabaergic_t4 vs Gabaergic_t3
# Gabaergic_t1 vs Gabaergic_t4
# Gabaergic_t1 vs Gabaergic_t2
# Gabaergic_t1 vs Gabaergic_t3
# Gabaergic_t2 vs Gabaergic_t4
# Gabaergic_t2 vs Gabaergic_t1
# Gabaergic_t2 vs Gabaergic_t3
# Gabaergic_t3 vs Gabaergic_t4
# Gabaergic_t3 vs Gabaergic_t1
# Gabaergic_t3 vs Gabaergic_t2
# Oligodendrocyte_t4 vs Oligodendrocyte_t1
# Oligodendrocyte_t4 vs Oligodendrocyte_t2
# Oligodendrocyte_t4 vs Oligodendrocyte_t3
# Oligodendrocyte_t1 vs Oligodendrocyte_t4
# Oligodendrocyte_t1 vs Oligodendrocyte_t2
# Oligodendrocyte_t1 vs Oligodendrocyte_t3
# Oligodendrocyte_t2 vs Oligodendrocyte_t4
# Oligodendrocyte_t2 vs Oligodendrocyte_t1 
# Oligodendrocyte_t2 vs Oligodendrocyte_t3
# Oligodendrocyte_t3 vs Oligodendrocyte_t4
# Oligodendrocyte_t3 vs Oligodendrocyte_t1
# Oligodendrocyte_t3 vs Oligodendrocyte_t2
# Astrocyte_t4 vs Astrocyte_t1
# Astrocyte_t4 vs Astrocyte_t2
# Astrocyte_t4 vs Astrocyte_t3
# Astrocyte_t1 vs Astrocyte_t4 
# Astrocyte_t1 vs Astrocyte_t2
# Astrocyte_t1 vs Astrocyte_t3
# Astrocyte_t2 vs Astrocyte_t4
# Astrocyte_t2 vs Astrocyte_t1
# Astrocyte_t2 vs Astrocyte_t3
# Astrocyte_t3 vs Astrocyte_t4
# 
