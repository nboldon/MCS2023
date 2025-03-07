

## Abundance boxplots normalized by total cell count in sample


setwd("/Volumes/DataBox/MCS2023/Cell_Abundance")

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

library(ggplot2)
library(viridis)

# Ensure 'cluster' is a factor with correct levels
boxData$cluster <- factor(boxData$cluster, levels = unique(boxData$cluster))  # Ensure 'cluster' is a factor

########################################

# Ensure clusters are ordered numerically
boxData$cluster <- factor(boxData$cluster, levels = unique(boxData$cluster[order(as.numeric(gsub("C", "", boxData$cluster)))]))

# Plot without mean points and apply the Viridis color palette
ggplot(boxData, aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16, 17, 18)) + 
  scale_color_viridis_d() +  # Apply Viridis color palette for discrete 'treatment' values
  facet_wrap(~cluster, scales = "free") +  # Facet by cluster (correct order)
  labs(
    color = "Treatment",  # Update the color legend title
    y = "Normalized Cell Abundance",  # Update the y-axis label
    x = "Cluster (C)",  # Add x-axis label
    title = "Cell Abundance by Cluster"  # Add plot title
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(color = "black", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10)  # Adjust legend text size
  ) +
  guides(shape = "none")  # Remove duplicate legend for 'shape'

# Save the plot
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_2024-12-09.jpg", device = "jpg", width = 14, height = 9)

########################################

# Ensure clusters are ordered numerically
boxData$cluster <- factor(boxData$cluster, levels = unique(boxData$cluster[order(as.numeric(gsub("C", "", boxData$cluster)))]))

# Plot with mean points and apply the Viridis color palette
ggplot(boxData, aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16, 17, 18)) + 
  scale_color_viridis_d() +  # Apply Viridis color palette for discrete 'treatment' values
  facet_wrap(~cluster, scales = "free") +  # Facet by cluster (correct order)
  geom_boxplot(alpha = 0.05) +  # Add transparency to boxplot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +  # Add mean points
  labs(
    color = "Treatment",  # Update the color legend title
    y = "Normalized Cell Abundance",  # Update the y-axis label
    x = "Cluster (C)",  # Add x-axis label
    title = "Cell Abundance by Cluster with Mean"  # Add plot title
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(color = "black", size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10)  # Adjust legend text size
  ) +
  guides(shape = "none")  # Remove duplicate legend for 'shape'

# Save the plot with mean points
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_wMean_2024-12-09.jpg", device = "jpg", width = 14, height = 9)




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



## Cell abundance boxplots by cell type, normalized by total cells in sample


# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Multiply the cellAbund values by 100 to convert percentages to representative values
boxData$cellAbund <- boxData$cellAbund * 100

# Create a new 'cell_type' column based on the 'cluster' column
boxData$cell_type <- NA

# Manually map clusters to cell types
boxData$cell_type[boxData$cluster %in% c("C18", "C19", "C21")] <- "Glutaminergic"
boxData$cell_type[boxData$cluster %in% c("C15", "C16", "C17", "C20", "C25")] <- "GlutPrecursor"
boxData$cell_type[boxData$cluster %in% c("C22", "C23")] <- "GABAergic"
boxData$cell_type[boxData$cluster %in% c("C10", "C11")] <- "Microglia"
boxData$cell_type[boxData$cluster %in% c("C2", "C3")] <- "Oligodendrocyte"
boxData$cell_type[boxData$cluster %in% c("C5", "C6")] <- "OligoPrecursor"
boxData$cell_type[boxData$cluster %in% c("C4")] <- "GlutOligo"
boxData$cell_type[boxData$cluster %in% c("C12", "C13", "C14")] <- "EndoVasc"
boxData$cell_type[boxData$cluster %in% c("C8")] <- "Astrocyte"
boxData$cell_type[boxData$cluster %in% c("C1", "C7", "C9")] <- "AstrocytePrecursor"
boxData$cell_type[boxData$cluster %in% c("C24")] <- "GlutAstro"



# Manually specify the order of cell types
cell_type_order <- c(
  "Glutaminergic", "GlutPrecursor", "GABAergic", "Microglia", 
  "Oligodendrocyte", "OligoPrecursor", "GlutOligo", "EndoVasc",
  "Astrocyte", "AstrocytePrecursor", "GlutAstro"
)

# Ensure that 'cell_type' is a factor with the desired order
boxData$cell_type <- factor(boxData$cell_type, levels = cell_type_order)

# Ensure 'treatment' is a factor
boxData$treatment <- factor(boxData$treatment)

# Now, create the plot with the specified order
abundance_boxplot <- ggplot(boxData, aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # Set shape for different treatments
  scale_color_manual(values = c("2N" = "#440154", "2N+" = "#31688EFF", "Ts" = "#35B779FF", "Ts+" = "#FDE725FF")) +  # Adjusted color choices
  facet_wrap(~cell_type, scales = "free") +  # Split into separate panels by cell_type
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  labs(
    color = "Treatment",  # Update the color legend title
    y = "Normalized Cell Abundance",  # Update the y-axis label
    x = "Cell Types",  # Add x-axis label
    title = "Cell Abundance by Cell Type with Mean"  # Add plot title
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adjust panel borders
    panel.background = element_rect(fill = "white", color = NA),  # Adjust panel background
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Adjust facet strip background
    strip.text = element_text(color = "black", size = 10),  # Adjust facet strip text size and color
    axis.text = element_text(size = 10),  # Adjust axis text size
    axis.title = element_text(size = 12)  # Adjust axis title text size
  ) +
  guides(shape = "none")  # Remove duplicate legend for 'shape'

# Save the plot
ggsave(
  filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_ByCellType_2024-12-09.jpg",
  plot = abundance_boxplot,
  device = "jpg",  
  width = 14,
  height = 9
)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


## Loop to create boxplots for all treatment group combinations by cell type 


# Load necessary libraries
library(ggpubr)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Multiply the cellAbund values by 100 to convert percentages to representative values
boxData$cellAbund <- boxData$cellAbund * 100

# Create a new 'cell_type' column based on the 'cluster' column
boxData$cell_type <- NA

# Manually map clusters to cell types
boxData$cell_type[boxData$cluster %in% c("C18", "C19", "C21")] <- "Glutaminergic"
boxData$cell_type[boxData$cluster %in% c("C15", "C16", "C17", "C20", "C25")] <- "GlutPrecursor"
boxData$cell_type[boxData$cluster %in% c("C22", "C23")] <- "GABAergic"
boxData$cell_type[boxData$cluster %in% c("C10", "C11")] <- "Microglia"
boxData$cell_type[boxData$cluster %in% c("C2", "C3")] <- "Oligodendrocyte"
boxData$cell_type[boxData$cluster %in% c("C5", "C6")] <- "OligoPrecursor"
boxData$cell_type[boxData$cluster %in% c("C4")] <- "GlutOligo"
boxData$cell_type[boxData$cluster %in% c("C12", "C13", "C14")] <- "EndoVasc"
boxData$cell_type[boxData$cluster %in% c("C8")] <- "Astrocyte"
boxData$cell_type[boxData$cluster %in% c("C1", "C7", "C9")] <- "AstrocytePrecursor"
boxData$cell_type[boxData$cluster %in% c("C24")] <- "GlutAstro"

# Manually specify the order of cell types
cell_type_order <- c(
  "Glutaminergic", "GlutPrecursor", "GABAergic", "Microglia", 
  "Oligodendrocyte", "OligoPrecursor", "GlutOligo", "EndoVasc",
  "Astrocyte", "AstrocytePrecursor", "GlutAstro"
)

# Ensure that 'cell_type' is a factor with the desired order
boxData$cell_type <- factor(boxData$cell_type, levels = cell_type_order)

# Ensure 'treatment' is a factor
boxData$treatment <- factor(boxData$treatment)

# Create all pairwise combinations of treatments
treatment_combinations <- combn(levels(boxData$treatment), 2, simplify = FALSE)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Iterate through each combination and save plots with dynamic file names
for (combo in treatment_combinations) {
  # Filter data for the current pairwise treatment combination
  filtered_data <- subset(boxData, treatment %in% combo)
  
  # Create the file name dynamically
  file_name <- paste0(
    "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_",
    paste(combo, collapse = "_vs_"),
    "_2024-12-09.pdf"
  )
  
  # Adjust the plot height based on the number of facets
  num_facets <- length(unique(filtered_data$cell_type))
  plot_height <- max(4, num_facets * 0.8)  # Adjust height dynamically
  
  # Create the boxplot
  abundance_boxplot <- ggplot(filtered_data, aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), cex = 1.2) + 
    scale_shape_manual(values = c(15, 16)) +  # Set shape for different treatments
    scale_color_manual(values = colors) +  # Use the fixed colors for both treatments
    facet_wrap(~cell_type, scales = "free", ncol = 2) +  # Split into separate panels by cell_type in specified order
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
    labs(
      color = "Treatment",  # Update the color legend title
      y = "Normalized Cell Abundance",  # Update the y-axis label
      x = "Cell Types",  # Add x-axis label
      title = paste("Cell Abundance by Cell Type with Mean: ", paste(combo, collapse = " vs. "))  # Add plot title
    ) + 
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adjust panel borders
      panel.background = element_rect(fill = "white", color = NA),  # Adjust panel background
      strip.background = element_rect(fill = "lightgray", color = "black"),  # Adjust facet strip background
      strip.text = element_text(color = "black", size = 10),  # Adjust facet strip text size and color
      axis.text = element_text(size = 10),  # Adjust axis text size
      axis.title = element_text(size = 12)  # Adjust axis title text size
    ) + 
    guides(shape = "none")  # Remove duplicate legend for 'shape'
  
  # Save the plot to a PDF file with adjusted height
  ggsave(
    filename = file_name,
    plot = abundance_boxplot,
    device = "pdf",
    width = 14,
    height = 17
  )
}



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



# Loop to create boxplots for genotype comparisons

library(ggpubr)
library(ggplot2)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Multiply the cellAbund values by 100 to convert percentages to representative values
boxData$cellAbund <- boxData$cellAbund * 100

# Create a new 'cell_type' column based on the 'cluster' column
boxData$cell_type <- NA

# Manually map clusters to cell types
boxData$cell_type[boxData$cluster %in% c("C18", "C19", "C21")] <- "Glutaminergic"
boxData$cell_type[boxData$cluster %in% c("C15", "C16", "C17", "C20", "C25")] <- "GlutPrecursor"
boxData$cell_type[boxData$cluster %in% c("C22", "C23")] <- "GABAergic"
boxData$cell_type[boxData$cluster %in% c("C10", "C11")] <- "Microglia"
boxData$cell_type[boxData$cluster %in% c("C2", "C3")] <- "Oligodendrocyte"
boxData$cell_type[boxData$cluster %in% c("C5", "C6")] <- "OligoPrecursor"
boxData$cell_type[boxData$cluster %in% c("C4")] <- "GlutOligo"
boxData$cell_type[boxData$cluster %in% c("C12", "C13", "C14")] <- "EndoVasc"
boxData$cell_type[boxData$cluster %in% c("C8")] <- "Astrocyte"
boxData$cell_type[boxData$cluster %in% c("C1", "C7", "C9")] <- "AstrocytePrecursor"
boxData$cell_type[boxData$cluster %in% c("C24")] <- "GlutAstro"

# Manually specify the order of cell types
cell_type_order <- c(
  "Glutaminergic", "GlutPrecursor", "GABAergic", "Microglia", 
  "Oligodendrocyte", "OligoPrecursor", "GlutOligo", "EndoVasc",
  "Astrocyte", "AstrocytePrecursor", "GlutAstro"
)

# Ensure that 'cell_type' is a factor with the desired order
boxData$cell_type <- factor(boxData$cell_type, levels = cell_type_order)

# Group the treatment groups
boxData$treatment_group <- as.factor(
  ifelse(boxData$treatment %in% c("Ts", "Ts+"), "TsCombined", 
         ifelse(boxData$treatment %in% c("2N", "2N+"), "2NCombined", boxData$treatment))
)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Plot with mean points and apply the Viridis color palette
ggplot(boxData, aes(x = treatment_group, y = cellAbund, shape = treatment_group, color = treatment_group)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16)) +  # Set shape for different treatments
  scale_color_manual(values = colors) +  # Use the fixed colors for both treatments
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type (correct order)
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +  # Add mean points
  labs(
    color = "Treatment Group",  # Update the color legend title
    y = "Normalized Cell Abundance",  # Update the y-axis label
    x = "Genotype",  # Add x-axis label
    title = "Cell Abundance by Genotype and Cell Type"  # Add plot title
  ) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adjust panel borders
    panel.background = element_rect(fill = "white", color = NA),  # Adjust panel background
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Adjust facet strip background
    strip.text = element_text(color = "black", size = 10),  # Adjust facet strip text size and color
    axis.text = element_text(size = 10),  # Adjust axis text size
    axis.title = element_text(size = 12),  # Adjust axis title text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10)  # Adjust legend text size
  ) + 
  guides(shape = "none")  # Remove duplicate legend for 'shape'

# Save the plot to a file
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_by_Genotype-CellType_2024-12-09.jpg", 
       device = "jpg", width = 14, height = 9)

cat("Plot saved as CellAbund_Boxplots_by_CellType_2024-12-09.jpg\n")


######################################
######################################



## Plots one treatment group with all cell types


library(ggpubr)
library(ggplot2)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Multiply the cellAbund values by 100 to convert percentages to representative values
boxData$cellAbund <- boxData$cellAbund * 100

# Create a new 'cell_type' column based on the 'cluster' column
boxData$cell_type <- NA

# Manually map clusters to cell types
boxData$cell_type[boxData$cluster %in% c("C18", "C19", "C21")] <- "Glutaminergic"
boxData$cell_type[boxData$cluster %in% c("C15", "C16", "C17", "C20", "C25")] <- "GlutPrecursor"
boxData$cell_type[boxData$cluster %in% c("C22", "C23")] <- "GABAergic"
boxData$cell_type[boxData$cluster %in% c("C10", "C11")] <- "Microglia"
boxData$cell_type[boxData$cluster %in% c("C2", "C3")] <- "Oligodendrocyte"
boxData$cell_type[boxData$cluster %in% c("C5", "C6")] <- "OligoPrecursor"
boxData$cell_type[boxData$cluster %in% c("C4")] <- "GlutOligo"
boxData$cell_type[boxData$cluster %in% c("C12", "C13", "C14")] <- "EndoVasc"
boxData$cell_type[boxData$cluster %in% c("C8")] <- "Astrocyte"
boxData$cell_type[boxData$cluster %in% c("C1", "C7", "C9")] <- "AstrocytePrecursor"
boxData$cell_type[boxData$cluster %in% c("C24")] <- "GlutAstro"

# Manually specify the order of cell types
cell_type_order <- c(
  "Glutaminergic", "GlutPrecursor", "GABAergic", "Microglia", 
  "Oligodendrocyte", "OligoPrecursor", "GlutOligo", "EndoVasc",
  "Astrocyte", "AstrocytePrecursor", "GlutAstro"
)

# Ensure that 'cell_type' is a factor with the desired order
boxData$cell_type <- factor(boxData$cell_type, levels = cell_type_order)

# Group the treatment groups
boxData$treatment_group <- as.factor(
  ifelse(boxData$treatment %in% c("Ts", "Ts+"), "TsCombined", 
         ifelse(boxData$treatment %in% c("2N", "2N+"), "2NCombined", boxData$treatment))
)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Loop to create boxplots for each treatment group comparison
unique_comparisons <- unique(boxData$treatment_group)
plot_height <- 8  # Define plot height

for (comparison in unique_comparisons) {
  
  # Filter data for the current comparison
  filtered_data <- subset(boxData, treatment_group == comparison)
  
  # Create the boxplot for all cell types in the current comparison
  abundance_boxplot <- ggplot(filtered_data, aes(
    x = cell_type, 
    y = cellAbund, 
    color = treatment_group
  )) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), cex = 1.2) + 
    stat_summary(fun = mean, geom = "point", size = 3, shape = 4, color = "black") + 
    scale_color_manual(values = colors) +
    labs(
      color = "Treatment Group",
      y = "Normalized Cell Abundance",
      x = "Cell Type",
      title = paste("Cell Abundance by Treatment Group: ", comparison)
    ) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save the plot for the current comparison
  file_name <- paste0(
    "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplot_", 
    comparison, "_2024-12-09.pdf"
  )
  
  ggsave(
    filename = file_name,
    plot = abundance_boxplot,
    device = "pdf",
    width = 14,
    height = plot_height
  )
}




######################################
######################################



# Loop to create boxplots for diet comparisons

library(ggpubr)
library(ggplot2)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Multiply the cellAbund values by 100 to convert percentages to representative values
boxData$cellAbund <- boxData$cellAbund * 100

# Create a new 'cell_type' column based on the 'cluster' column
boxData$cell_type <- NA

# Manually map clusters to cell types
boxData$cell_type[boxData$cluster %in% c("C18", "C19", "C21")] <- "Glutaminergic"
boxData$cell_type[boxData$cluster %in% c("C15", "C16", "C17", "C20", "C25")] <- "GlutPrecursor"
boxData$cell_type[boxData$cluster %in% c("C22", "C23")] <- "GABAergic"
boxData$cell_type[boxData$cluster %in% c("C10", "C11")] <- "Microglia"
boxData$cell_type[boxData$cluster %in% c("C2", "C3")] <- "Oligodendrocyte"
boxData$cell_type[boxData$cluster %in% c("C5", "C6")] <- "OligoPrecursor"
boxData$cell_type[boxData$cluster %in% c("C4")] <- "GlutOligo"
boxData$cell_type[boxData$cluster %in% c("C12", "C13", "C14")] <- "EndoVasc"
boxData$cell_type[boxData$cluster %in% c("C8")] <- "Astrocyte"
boxData$cell_type[boxData$cluster %in% c("C1", "C7", "C9")] <- "AstrocytePrecursor"
boxData$cell_type[boxData$cluster %in% c("C24")] <- "GlutAstro"

# Manually specify the order of cell types
cell_type_order <- c(
  "Glutaminergic", "GlutPrecursor", "GABAergic", "Microglia", 
  "Oligodendrocyte", "OligoPrecursor", "GlutOligo", "EndoVasc",
  "Astrocyte", "AstrocytePrecursor", "GlutAstro"
)

# Ensure that 'cell_type' is a factor with the desired order
boxData$cell_type <- factor(boxData$cell_type, levels = cell_type_order)

# Group the treatment groups
boxData$treatment_group <- as.factor(
  ifelse(boxData$treatment %in% c("2N+", "Ts+"), "MCS", 
         ifelse(boxData$treatment %in% c("2N", "Ts"), "CON", boxData$treatment))
)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Plot with mean points and apply the Viridis color palette
ggplot(boxData, aes(x = treatment_group, y = cellAbund, shape = treatment_group, color = treatment_group)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16)) +  # Set shape for different treatments
  scale_color_manual(values = colors) +  # Use the fixed colors for both treatments
  facet_wrap(~cell_type, scales = "free") +  # Facet by cell type (correct order)
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +  # Add mean points
  labs(
    color = "Treatment Group",  # Update the color legend title
    y = "Normalized Cell Abundance",  # Update the y-axis label
    x = "Maternal Diet",  # Add x-axis label
    title = "Cell Abundance by Diet and Cell Type"  # Add plot title
  ) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adjust panel borders
    panel.background = element_rect(fill = "white", color = NA),  # Adjust panel background
    strip.background = element_rect(fill = "lightgray", color = "black"),  # Adjust facet strip background
    strip.text = element_text(color = "black", size = 10),  # Adjust facet strip text size and color
    axis.text = element_text(size = 10),  # Adjust axis text size
    axis.title = element_text(size = 12),  # Adjust axis title text size
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10)  # Adjust legend text size
  ) + 
  guides(shape = "none")  # Remove duplicate legend for 'shape'

# Save the plot to a file
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_by_Diet-CellType_2024-12-09.jpg", 
       device = "jpg", width = 14, height = 9)

cat("Plot saved as CellAbund_Boxplots_by_CellType_2024-12-09.jpg\n")



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



### Additional viridis colors 


## "viridis": This is the default color palette, ranging from dark purple to yellow.
viridis_colors <- viridis(10, option = "viridis")
# Print the colors for magma in hex format
print(viridis_colors)
"#440154FF" "#482878FF" "#3E4A89FF" "#31688EFF" "#26828EFF" "#1F9E89FF" "#35B779FF" "#6DCD59FF" "#B4DE2CFF" "#FDE725FF"

#"2N" (#440154) → R: 68, G: 1, B: 84
#"2N+" (#31688EFF) → R: 49, G: 104, B: 142
#"Ts" (#35B779FF) → R: 53, G: 183, B: 121
#"Ts+" (#FDE725FF) → R: 253, G: 231, B: 37
        

## "magma": A color palette with a range from dark purple to bright yellow, but with a more subdued yellow.
# Generate a 10-color palette for 'magma'
magma_colors <- viridis(10, option = "magma")
print(magma_colors)
 "#000004FF" "#180F3EFF" "#451077FF" "#721F81FF" "#9F2F7FFF" "#CD4071FF" "#F1605DFF" "#FD9567FF" "#FEC98DFF" "#FCFDBFFF"

## "plasma": A bright and high-contrast palette, ranging from deep purple to yellow-orange.
# Generate a 10-color palette for 'plasma'
plasma_colors <- viridis(10, option = "plasma")
print(plasma_colors)
 "#0D0887FF" "#47039FFF" "#7301A8FF" "#9C179EFF" "#BD3786FF" "#D8576BFF" "#ED7953FF" "#FA9E3BFF" "#FDC926FF" "#F0F921FF"

## "inferno": A color palette that starts from dark purple and ends in bright yellow, but with more intense and warmer colors than viridis.
 # Generate a 10-color palette for 'inferno'
 inferno_colors <- viridis(10, option = "inferno")
 print(inferno_colors)
 "#000004FF" "#1B0C42FF" "#4B0C6BFF" "#781C6DFF" "#A52C60FF" "#CF4446FF" "#ED6925FF" "#FB9A06FF" "#F7D03CFF" "#FCFFA4FF"

## "cividis": A color-blind-friendly palette, designed for those with color vision deficiencies, ranging from dark blue to bright yellow.
 cividis_colors <- viridis(10, option = "cividis")
 print(cividis_colors)
 "#00204DFF" "#00336FFF" "#39486BFF" "#575C6DFF" "#707173FF" "#8A8779FF" "#A69D75FF" "#C4B56CFF" "#E4CF5BFF" "#FFEA46FF"






