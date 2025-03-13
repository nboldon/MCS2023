library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(ggplot2)

#Additional setup
setwd("/Volumes/DataBox/ProjMCS6/Violin_Plots")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)

projMCS6

######################################################
######################################################
######################################################

# Subset by Cluster
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

############################################
############################################
############################################

## Cluster 1 

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

p1 <- plotGroups(
  ArchRProj = C25ArchRSubset, 
  groupBy = "treatment", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p1

C25_nFrags-byTx_2024-08-13

####################################
####################################
####################################

p2 <- plotGroups(
  ArchRProj = C25ArchRSubset,
  groupBy = "treatment",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

C25_TSS-byTx_2024-08-13

####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################
####################################


## Loop for all clusters, prints in 1 pdf


# Load necessary libraries
library(ArchR)
library(ggplot2)
library(cowplot)
library(viridis)

# Set working directory
setwd("/Volumes/DataBox/ProjMCS6/Violin_Plots")

# Load ArchR Project
projMCS6 <- loadArchRProject(path = "/Volumes/DataBox/Save-ProjMCS6", force = FALSE, showLogo = FALSE)
          
            
            # Define the color palette using viridis
            # You can choose from the default viridis palette, or other palettes such as "magma", "inferno", "plasma"
            treatment_colors <- viridis(4, option = "D")  # Generates a viridis palette with 4 colors
            
            # Open a PDF to save the plots
            pdf("QC-TSS-UniqFrags_byTx-Cluster_ViolinPlots_2024-12-16.pdf", width = 16, height = 12)
            
            # Initialize lists to store plots
            tss_plots <- list()
            frag_plots <- list()
            
            # Loop through all clusters
            # Loop through all clusters
            for (cluster in sort(unique(as.numeric(gsub("C", "", projMCS6$Clusters))))) {
              
              # Subset the project for the current cluster
              subsetProj <- projMCS6[projMCS6$Clusters %in% paste0("C", cluster)]
              
              # Assign treatments to the subset (this part is unchanged)
              treatment <- subsetProj$Sample
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
              subsetProj$treatment <- treatment
              
              # Create the TSS Enrichment plot with viridis colors for both fill and outline
              tss_plot <- plotGroups(
                ArchRProj = subsetProj,
                groupBy = "treatment",
                colorBy = "cellColData",
                name = "TSSEnrichment",
                plotAs = "violin",
                alpha = 0.4,
                addBoxPlot = TRUE
              ) + scale_color_viridis(discrete = TRUE) + 
                scale_fill_viridis(discrete = TRUE) +   # Set fill color scale
                ggtitle(paste("Cluster", cluster, "- TSSEnrichment by Treatment")) +
                theme(legend.position = "none")  # Remove legend
              
              # Add the plot to the TSS plot list
              tss_plots[[as.character(cluster)]] <- tss_plot
              
              # Create the log10(nFrags) plot with viridis colors for both fill and outline
              frag_plot <- plotGroups(
                ArchRProj = subsetProj,
                groupBy = "treatment",
                colorBy = "cellColData",
                name = "log10(nFrags)",
                plotAs = "violin",
                alpha = 0.4,
                addBoxPlot = TRUE
              ) + scale_color_viridis(discrete = TRUE) + 
                scale_fill_viridis(discrete = TRUE) +   # Set fill color scale
                ggtitle(paste("Cluster", cluster, "- log10(nFrags) by Treatment")) +
                theme(legend.position = "none")  # Remove legend
              
              # Add the plot to the Frags plot list
              frag_plots[[as.character(cluster)]] <- frag_plot
            }
            
            # Combine TSS plots into a single page with lexicographical order
            tss_combined <- plot_grid(
              plotlist = tss_plots[order(as.numeric(gsub("C", "", names(tss_plots))))], 
              ncol = 5, 
              align = "v"
            )
            print(tss_combined)
            
            # Combine Frags plots into a single page with lexicographical order
            frag_combined <- plot_grid(
              plotlist = frag_plots[order(as.numeric(gsub("C", "", names(frag_plots))))], 
              ncol = 5, 
              align = "v"
            )
            print(frag_combined)
            
            # Close the PDF
            dev.off()
            
