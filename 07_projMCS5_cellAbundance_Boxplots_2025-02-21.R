

## Cell abundance by cluster


## Abundance boxplots normalized by total cell count in sample


setwd("/Volumes/DataBox/MCS2023/Cell_Abundance")

# Load required libraries
library(ggplot2)
library(viridis)
library(rstatix)
library(ggpubr)
library(dplyr)
library(tibble)
library(purrr)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell_Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

# Ensure clusters are ordered numerically
boxData$cluster <- factor(boxData$cluster,
                          levels = unique(boxData$cluster[order(as.numeric(gsub("C", "", boxData$cluster)))]))

# Ensure 'treatment' is a factor
boxData$treatment <- factor(boxData$treatment)

# Perform statistical analysis
# Run ANOVA for each cluster
stats_results <- boxData %>%
  group_by(cluster) %>%
  anova_test(cellAbund ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# Perform pairwise comparisons for each cluster
pairwise_tests <- boxData %>%
  group_by(cluster) %>%
  pairwise_t_test(cellAbund ~ treatment, p.adjust.method = "bonferroni") %>%
  add_significance()

# Save statistical results to CSV files
write.csv(stats_results, "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_ANOVA_Stats_byCluster_2025-02-21.csv", row.names = FALSE)
write.csv(pairwise_tests, "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Pairwise_Stats_byCluster_2025-02-21.csv", row.names = FALSE)

# Filter for significant comparisons only and create comparison list by cluster
sig_comparisons <- pairwise_tests %>%
  filter(p.adj.signif != "ns") %>%
  group_by(cluster) %>%
  summarise(comparisons = list(map2(group1, group2, c))) %>%
  {setNames(.$comparisons, .$cluster)}

# Calculate y-position for significance brackets with more space
y_positions <- boxData %>%
  group_by(cluster) %>%
  summarise(max_y = max(cellAbund) * 1.4)

# Function to create plot with improved significance markers and comparison brackets
create_improved_plot <- function(include_mean = FALSE) {
  # Create base plot
  p <- ggplot(boxData, aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), cex = 1.2) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_color_viridis_d() +
    facet_wrap(~cluster, scales = "free")
  
  # Add comparison brackets for each cluster only if significant comparisons exist
  for(cluster_name in names(sig_comparisons)) {
    if(length(sig_comparisons[[cluster_name]]) > 0) {
      cluster_data <- boxData %>% filter(cluster == cluster_name)
      p <- p + stat_compare_means(
        data = cluster_data,
        comparisons = sig_comparisons[[cluster_name]],
        method = "t.test",
        label = "p.signif",
        label.y = y_positions$max_y[y_positions$cluster == cluster_name],
        bracket.size = 0.5,
        step.increase = 0.1,
        tip.length = 0.01
      )
    }
  }
  
  # Add remaining plot elements
  p <- p + labs(
    color = "Treatment",
    y = "Normalized Cell Abundance",
    x = "Cluster (C)",
    title = if(include_mean) "Cell Abundance by Cluster with Mean and Statistics" else "Cell Abundance by Cluster with Statistics"
  ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      strip.text = element_text(color = "black", size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    guides(shape = "none")
  
  if(include_mean) {
    p <- p + stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black")
  }
  
  return(p)
}

# Create and save both plots
p1 <- create_improved_plot(include_mean = FALSE)
p2 <- create_improved_plot(include_mean = TRUE)

# Save plots with higher resolution
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_Stats_byCluster_2025-02-21.jpg",
       plot = p1, device = "jpg", width = 14, height = 9, dpi = 300)

ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_Stats_wMean_byCluster_2025-02-21.jpg",
       plot = p2, device = "jpg", width = 14, height = 9, dpi = 300)


###################################################
###################################################
###################################################
###################################################
###################################################



## Cell abundance by cell type


## Load necessary libraries
library(ggpubr)
library(rstatix)
library(dplyr)

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

# Function to perform statistical tests and return results
perform_stats <- function(data) {
  # Perform one-way ANOVA for each cell type
  stats_results <- data %>%
    group_by(cell_type) %>%
    anova_test(cellAbund ~ treatment) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  
  # If ANOVA is significant, perform pairwise t-tests
  pairwise_tests <- data %>%
    group_by(cell_type) %>%
    pairwise_t_test(cellAbund ~ treatment, p.adjust.method = "bonferroni") %>%
    add_significance()
  
  return(list(anova = stats_results, pairwise = pairwise_tests))
}

# Iterate through each combination and save plots with dynamic file names
for (combo in treatment_combinations) {
  # Filter data for the current pairwise treatment combination
  filtered_data <- subset(boxData, treatment %in% combo)
  
  # Perform statistical tests
  stats_results <- perform_stats(filtered_data)
  
  # Calculate maximum value for each cell type for better y-axis scaling
  cell_type_max <- filtered_data %>%
    group_by(cell_type) %>%
    summarise(max_val = max(cellAbund))
  
  # Filter for significant comparisons only
  sig_comparisons <- stats_results$pairwise %>%
    filter(p.adj.signif != "ns") %>%
    group_by(cell_type) %>%
    summarise(comparisons = list(map2(group1, group2, c))) %>%
    {setNames(.$comparisons, .$cell_type)}
  
  # Calculate better y-positions for significance brackets
  y_positions <- filtered_data %>%
    group_by(cell_type) %>%
    summarise(max_y = max(cellAbund) * 1.05) # Reduced from 1.2 to 1.05
  
  # Calculate number of columns based on number of unique cell types
  n_cell_types <- length(unique(filtered_data$cell_type))
  n_cols <- min(4, ceiling(sqrt(n_cell_types))) # Adjust the max number of columns as needed
  n_rows <- ceiling(n_cell_types / n_cols)
  
  # Create the file name dynamically
  file_name <- paste0(
    "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_",
    paste(combo, collapse = "_vs_"),
    "_cellType_2025-02-21.pdf"
  )
  
  # Create the boxplot with statistical annotations
  abundance_boxplot <- ggplot(filtered_data, 
                              aes(x = treatment, y = cellAbund, shape = treatment, color = treatment)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), cex = 1.2) +
    scale_shape_manual(values = c(15, 16)) +
    scale_color_manual(values = colors) +
    # Use dynamic number of columns for facet_wrap
    facet_wrap(~cell_type, scales = "free_y", ncol = n_cols) +
    # Add explicit y-axis control
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)) # Controls padding at bottom (0) and top (0.1)
    )
  
  # Add significant comparison brackets for each cell type
  for(cell_type_name in names(sig_comparisons)) {
    if(length(sig_comparisons[[cell_type_name]]) > 0) {
      cell_data <- filtered_data %>% filter(cell_type == cell_type_name)
      abundance_boxplot <- abundance_boxplot + 
        stat_compare_means(
          data = cell_data,
          comparisons = sig_comparisons[[cell_type_name]],
          method = "t.test",
          label = "p.signif",
          label.y = y_positions$max_y[y_positions$cell_type == cell_type_name],
          bracket.size = 0.5,
          step.increase = 0.1,
          tip.length = 0.01
        )
    }
  }
  
  # Add remaining plot elements
  abundance_boxplot <- abundance_boxplot +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
    labs(
      color = "Treatment",
      y = "Normalized Cell Abundance",
      x = "Cell Types",
      title = paste("Cell Abundance by Cell Type with Mean: ", paste(combo, collapse = " vs. "))
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      strip.text = element_text(color = "black", size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Center the title
      plot.margin = margin(10, 10, 10, 10) # Add some margin around the plot
    ) +
    guides(shape = "none")
  
  # Save the plot to a PDF file with dynamic dimensions
  ggsave(
    filename = file_name,
    plot = abundance_boxplot,
    device = "pdf",
    width = min(14, 3.5 * n_cols), # 3.5 inches per column, max 14
    height = min(17, 3.5 * n_rows) # 3.5 inches per row, max 17
  )
  
  # Save statistical results to CSV files
  stats_anova_file <- gsub("\\.pdf$", "_ANOVA_cellType_stats_2025-02-21.csv", file_name)
  stats_pairwise_file <- gsub("\\.pdf$", "_pairwise_cellType_stats_2025-02-21.csv", file_name)
  
  # Save ANOVA results
  write.csv(stats_results$anova, stats_anova_file, row.names = FALSE)
  
  # Save pairwise test results
  write.csv(stats_results$pairwise, stats_pairwise_file, row.names = FALSE)
}
###################################################
###################################################
###################################################
###################################################
###################################################


## Cell abundance by genotype


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(rstatix)
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

# Group the treatment groups
boxData$treatment_group <- as.factor(
  ifelse(boxData$treatment %in% c("Ts", "Ts+"), "TsCombined",
         ifelse(boxData$treatment %in% c("2N", "2N+"), "2NCombined", boxData$treatment))
)

# Perform statistical analysis
# Run t-tests for each cell type and create a more detailed results dataframe
stats_results <- boxData %>%
  group_by(cell_type) %>%
  t_test(cellAbund ~ treatment_group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  select(cell_type,
         group1,
         group2,
         n1,
         n2,
         statistic,
         df,
         p,
         p.adj,
         p.adj.signif) %>%
  mutate(
    mean_group1 = sapply(seq_len(n()), function(i) {
      mean(boxData$cellAbund[boxData$treatment_group == group1[i] &
                               boxData$cell_type == cell_type[i]])
    }),
    mean_group2 = sapply(seq_len(n()), function(i) {
      mean(boxData$cellAbund[boxData$treatment_group == group2[i] &
                               boxData$cell_type == cell_type[i]])
    }),
    sd_group1 = sapply(seq_len(n()), function(i) {
      sd(boxData$cellAbund[boxData$treatment_group == group1[i] &
                             boxData$cell_type == cell_type[i]])
    }),
    sd_group2 = sapply(seq_len(n()), function(i) {
      sd(boxData$cellAbund[boxData$treatment_group == group2[i] &
                             boxData$cell_type == cell_type[i]])
    })
  )

# Save statistical results to a CSV file
write.csv(stats_results, 
          file = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Genotype_Stats_2025-02-21.csv",
          row.names = FALSE)

# Calculate number of columns based on number of unique cell types
n_cell_types <- length(unique(boxData$cell_type))
n_cols <- min(4, ceiling(sqrt(n_cell_types)))
n_rows <- ceiling(n_cell_types / n_cols)

# Create a function to manually determine appropriate y-position for each cell type
get_optimal_y_position <- function(data) {
  cell_type_stats <- data %>%
    group_by(cell_type) %>%
    summarise(
      max_val = max(cellAbund),
      q3 = quantile(cellAbund, 0.75),
      iqr = IQR(cellAbund)
    )
  
  # Set y-position just slightly above the maximum value
  # Use 1.02 multiplier (very close to max) for better visualization
  cell_type_stats$y_pos <- cell_type_stats$max_val * 1.02
  
  return(cell_type_stats)
}

# Get optimal y-positions
y_position_data <- get_optimal_y_position(boxData)

# Filter for significant comparisons
sig_tests <- stats_results %>% 
  filter(p.adj.signif != "ns") %>%
  select(cell_type, group1, group2, p.adj.signif)

# Create custom comparisons list for each cell type
comparison_list <- lapply(unique(boxData$cell_type), function(ct) {
  filtered <- sig_tests %>% filter(cell_type == ct)
  if(nrow(filtered) > 0) {
    return(list(c(filtered$group1[1], filtered$group2[1])))
  } else {
    return(NULL)
  }
})
names(comparison_list) <- unique(boxData$cell_type)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Create the base plot
p <- ggplot(boxData, aes(x = treatment_group, y = cellAbund, shape = treatment_group, color = treatment_group)) + 
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, color = "black") +  # Box outline only, no fill
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),  
              alpha = 0.6, size = 1.5, aes(color = treatment_group)) +  # Use color only for jitter
  scale_color_manual(values = colors) +  # Apply color scheme to color only
  facet_wrap(~cell_type, scales = "free_y", ncol = n_cols) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "black") + 
  labs(
    color = "Treatment Group",  # Keep only the color legend
    y = "Normalized Cell Abundance",
    x = "Genotype",
    title = "Cell Abundance by Genotype and Cell Type with Statistics"
  ) +
  guides(
    fill = "none",  # Remove the fill legend
    shape = "none"  # Remove the shape legend if you don't want it either
  )


# Add the statistical annotations manually for each cell type
for(ct in unique(boxData$cell_type)) {
  # Skip if no significant comparison
  if(is.null(comparison_list[[ct]])) next
  
  # Get the y position for this cell type
  y_pos <- y_position_data$y_pos[y_position_data$cell_type == ct]
  
  # Add stat annotation for this cell type only
  p <- p + stat_compare_means(
    data = boxData %>% filter(cell_type == ct),
    comparisons = comparison_list[[ct]],
    method = "t.test",
    label = "p.signif",
    label.y = y_pos,
    bracket.size = 0.4,
    step.increase = 0.05,
    tip.length = 0.01
  )
}

# Apply more aggressive y-axis limits for each facet
p <- p + scale_y_continuous(
  # This dramatically reduces white space by setting expansion to nearly zero
  expand = expansion(mult = c(0.05, 0.1))
)

# Apply improved theme
p <- p + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  panel.background = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "lightgray", color = "black"),
  strip.text = element_text(color = "black", size = 10, face = "bold"),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12, face = "bold"),
  legend.position = "right",
  legend.title = element_text(size = 12, face = "bold"),
  legend.text = element_text(size = 10),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  plot.margin = margin(10, 10, 10, 10),
  # Add tight spacing between facets
  panel.spacing = unit(0.5, "lines")
)

# Save the plot to a file with dynamic dimensions
ggsave(
  filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_by_Genotype-CellType_Stats_2025-02-21.pdf",
  plot = p, 
  device = "pdf", 
  width = min(14, 3.5 * n_cols),
  height = min(9, 3.0 * n_rows),  # Slightly smaller height per row
  limitsize = FALSE
)

# Also save as jpg if needed
ggsave(
  filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_by_Genotype-CellType_Stats_2025-02-21.jpg",
  plot = p, 
  device = "jpg", 
  width = min(14, 3.5 * n_cols),
  height = min(9, 3.0 * n_rows),
  dpi = 300
)



###################################################
###################################################
###################################################
###################################################
###################################################


## Cell abundance by diet


library(ggpubr)
library(ggplot2)
library(rstatix)

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

# Perform statistical analysis
# Run t-tests for each cell type and create a more detailed results dataframe
stats_results <- boxData %>%
  group_by(cell_type) %>%
  t_test(cellAbund ~ treatment_group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  select(cell_type, 
         group1, 
         group2, 
         n1, 
         n2, 
         statistic, 
         df, 
         p, 
         p.adj, 
         p.adj.signif) %>%
  mutate(
    mean_group1 = sapply(seq_len(n()), function(i) {
      mean(boxData$cellAbund[boxData$treatment_group == group1[i] & 
                               boxData$cell_type == cell_type[i]])
    }),
    mean_group2 = sapply(seq_len(n()), function(i) {
      mean(boxData$cellAbund[boxData$treatment_group == group2[i] & 
                               boxData$cell_type == cell_type[i]])
    }),
    sd_group1 = sapply(seq_len(n()), function(i) {
      sd(boxData$cellAbund[boxData$treatment_group == group1[i] & 
                             boxData$cell_type == cell_type[i]])
    }),
    sd_group2 = sapply(seq_len(n()), function(i) {
      sd(boxData$cellAbund[boxData$treatment_group == group2[i] & 
                             boxData$cell_type == cell_type[i]])
    })
  )

# Save statistical results to a CSV file
write.csv(stats_results, 
          file = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Diet_Stats_2025-02-21.csv", 
          row.names = FALSE)

# Calculate y-position for significance brackets
y_positions <- boxData %>%
  group_by(cell_type) %>%
  summarise(max_y = max(cellAbund) * 1.1)

# Define the colors you want to use
colors <- c("#781C6DFF", "#1F9E89FF")

# Create the plot with statistical annotations
p <- ggplot(boxData, aes(x = treatment_group, y = cellAbund, shape = treatment_group, color = treatment_group)) + 
  geom_boxplot(position = position_dodge(0.8)) +
  geom_jitter(position = position_dodge(0.8), cex = 1.2) +
  scale_shape_manual(values = c(15, 16)) +
  scale_color_manual(values = colors) +
  facet_wrap(~cell_type, scales = "free") +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.y = y_positions$max_y
  ) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
  labs(
    color = "Treatment Group",
    y = "Normalized Cell Abundance",
    x = "Maternal Diet",
    title = "Cell Abundance by Diet and Cell Type with Statistics"
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
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) + 
  guides(shape = "none")

# Save the plot to a file
ggsave(filename = "/Volumes/DataBox/MCS2023/Cell_Abundance/CellAbund_Boxplots_by_Diet-CellType_Stats_2025-02-21.jpg", 
       plot = p, device = "jpg", width = 14, height = 9)

# Print the completion message
cat("Plot saved with statistical analysis\n")


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






