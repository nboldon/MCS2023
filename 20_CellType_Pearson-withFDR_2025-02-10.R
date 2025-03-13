


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_Cluster_zscore-Stats")


library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(patchwork)


## Complete for each cell type: astro, gaba, glut, oligo, endo, microglia

## Complete for each behavioral test: TaskA-45, TaskB-59, TaskC-48, TaskD-149



# Read the combined data
df <- read.csv("astro_TaskA-45_StatsCombined_2025-02-03.csv")

# Define metrics
behavioral_metrics <- c("Correct_Response", "Premature_Response",
                        "Missed_Response_Window", "Wrong_Choice")

# Get gene columns
gene_columns <- names(df)[!names(df) %in% 
                            c("Sample", "Treatment", "Correct_Response",
                              "Premature_Response", "Missed_Response_Window",
                              "Wrong_Choice")]

# Create scatter plot function (same as before)
create_gene_scatter <- function(data, gene, behavioral_metric) {
  plot_data <- data %>%
    select(all_of(c("Treatment", behavioral_metric, gene))) %>%
    drop_na()
  
  correlation <- tryCatch({
    cor.test(plot_data[[gene]], plot_data[[behavioral_metric]])
  }, error = function(e) list(estimate = NA, p.value = NA))
  
  p <- ggplot(plot_data, aes_string(x = gene, y = behavioral_metric, color = "Treatment", shape = "Treatment")) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis_d() +
    scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +
    theme_minimal() +
    labs(
      title = paste("Astro", gene, "vs", "TaskA-45", behavioral_metric),
      x = paste(gene, "Z-Score"),
      y = behavioral_metric,
      subtitle = sprintf("Pearson's R: %.3f (p = %.3f)",
                         correlation$estimate,
                         correlation$p.value)
    ) +
    theme(plot.subtitle = element_text(size = 10))
  
  return(p)
}

# Calculate correlations and p-values
correlation_results <- data.frame()
for (gene in gene_columns) {
  for (metric in behavioral_metrics) {
    test <- cor.test(df[[gene]], df[[metric]])
    correlation_results <- rbind(correlation_results, 
                                 data.frame(Gene = gene,
                                            Metric = metric,
                                            Correlation = test$estimate,
                                            P_value = test$p.value))
  }
}

# Calculate FDR-adjusted p-values
correlation_results$FDR_adjusted_p <- p.adjust(correlation_results$P_value, method = "BH")

# Create a more readable output
summary_df <- correlation_results %>%
  arrange(FDR_adjusted_p) %>%
  mutate(
    Significant = FDR_adjusted_p < 0.05,
    Correlation = round(Correlation, 3),
    P_value = format(P_value, scientific = TRUE, digits = 12),
    FDR_adjusted_p = format(FDR_adjusted_p, scientific = TRUE, digits = 12)
  )

# Save results
write.csv(summary_df, 
          file = "astro_TaskA-45_CorrelationStats_withFDR_2025-02-12.csv", 
          row.names = FALSE)

# Create plots (modified to include FDR-adjusted p-values)
plot_list <- list()
for (gene in gene_columns) {
  for (metric in behavioral_metrics) {
    # Get correlation stats
    stats <- correlation_results %>%
      filter(Gene == gene & Metric == metric)
    
    fdr_value <- as.numeric(stats$FDR_adjusted_p)
    
    # Create plot data
    plot_data <- df %>%
      select(all_of(c("Treatment", metric, gene))) %>%
      drop_na()
    
    # Generate base plot
    p <- ggplot(plot_data, aes_string(x = gene, y = metric, color = "Treatment", shape = "Treatment")) +
      geom_point(size = 3, alpha = 0.7) +
      scale_color_viridis_d() +
      scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +
      theme_minimal() +
      labs(
        title = paste("Astro", gene, "vs", "TaskA-45", metric),
        x = paste(gene, "Z-Score"),
        y = metric,
        subtitle = sprintf("R: %.3f (p = %.3e, FDR = %.3e)", 
                           stats$Correlation, 
                           as.numeric(stats$P_value), 
                           fdr_value)
      ) +
      theme(plot.subtitle = element_text(size = 10))
    
    # Highlight the p-value if significant
    if (fdr_value < 0.05) {
      p <- p + annotate(
        "text", x = min(plot_data[[gene]], na.rm = TRUE), 
        y = max(plot_data[[metric]], na.rm = TRUE), 
        label = sprintf("Significant: FDR = %.3e", fdr_value),
        color = "purple", fontface = "bold", hjust = 0, vjust = 1, size = 5
      )
    }
    
    plot_list[[paste(gene, metric)]] <- p
  }
}


# Arrange and save plots
final_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(title = "Astro Gene Z-Scores vs Task A-45 Behavioral Metrics")

ggsave("astro-zscore_TaskA-45_scatterplots_withFDR_2025-02-12.pdf", 
       final_plot,
       width = 20, 
       height = 5 * length(gene_columns),
       limitsize = FALSE)



################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################
################################################




### loop does not run 

## Loop for allfiles with the same file name ending
# Highlights sig FDR p-values on graphs

library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(patchwork)

# Get list of all combined stats files
files <- list.files(pattern = "*StatsCombined.*2025-02-03.csv")

# Define behavioral metrics
behavioral_metrics <- c("Correct_Response", "Premature_Response",
                        "Missed_Response_Window", "Wrong_Choice")

# Function to create scatter plot with highlighted significant points
create_gene_scatter <- function(data, gene, behavioral_metric, cell_type, task, stats) {
  plot_data <- data %>%
    select(all_of(c("Treatment", behavioral_metric, gene))) %>%
    drop_na()
  
  # Add significance information to plot data
  is_significant <- as.numeric(stats$FDR_adjusted_p) < 0.05
  
  p <- ggplot(plot_data, aes_string(x = gene, y = behavioral_metric)) +
    # Add points with different aesthetics based on significance
    geom_point(aes(color = Treatment, shape = Treatment),
               size = 3, 
               alpha = 0.7,
               # Add black border to significant points
               stroke = if(is_significant) 1.5 else 0.5,
               color = if(is_significant) "black" else NA) +
    # Add colored points on top
    geom_point(aes(color = Treatment, shape = Treatment),
               size = 2.5,
               alpha = 0.7) +
    scale_color_viridis_d() +
    scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +
    theme_minimal() +
    labs(
      title = paste(cell_type, gene, "vs", task, behavioral_metric),
      x = paste(gene, "Z-Score"),
      y = behavioral_metric,
      subtitle = sprintf("R: %.3f (p = %.3e, FDR = %.3e)%s",
                         stats$Correlation,
                         as.numeric(stats$P_value),
                         as.numeric(stats$FDR_adjusted_p),
                         if(is_significant) " *" else "")  # Add asterisk if significant
    ) +
    theme(plot.subtitle = element_text(size = 10))
  
  # Add significance annotation if significant
  if(is_significant) {
    p <- p + annotate("text", 
                      x = min(plot_data[[gene]], na.rm = TRUE),
                      y = max(plot_data[[behavioral_metric]], na.rm = TRUE),
                      label = "FDR < 0.05",
                      hjust = 0,
                      vjust = 1,
                      color = "red",
                      fontface = "bold")
  }
  
  return(p)
}

# Loop through each file
for (file in files) {
  # Extract cell type and task from filename
  cell_type <- gsub("(.*)_Task.*", "\\1", file)
  task <- gsub(".*_(Task[A-D]-[0-9]+)_.*", "\\1", file)
  
  # Read the data
  df <- read.csv(file)
  
  # Get gene columns for this dataset
  gene_columns <- names(df)[!names(df) %in% 
                              c("Sample", "Treatment", "Correct_Response",
                                "Premature_Response", "Missed_Response_Window",
                                "Wrong_Choice")]
  
  # Calculate correlations and p-values
  correlation_results <- data.frame()
  for (gene in gene_columns) {
    for (metric in behavioral_metrics) {
      test <- cor.test(df[[gene]], df[[metric]])
      correlation_results <- rbind(correlation_results, 
                                   data.frame(Gene = gene,
                                              Metric = metric,
                                              Correlation = test$estimate,
                                              P_value = test$p.value))
    }
  }
  
  # Calculate FDR-adjusted p-values
  correlation_results$FDR_adjusted_p <- p.adjust(correlation_results$P_value, method = "BH")
  
  # Create summary dataframe
  summary_df <- correlation_results %>%
    arrange(FDR_adjusted_p) %>%
    mutate(
      Significant = FDR_adjusted_p < 0.05,
      Correlation = round(Correlation, 3),
      P_value = format(P_value, scientific = TRUE, digits = 3),
      FDR_adjusted_p = format(FDR_adjusted_p, scientific = TRUE, digits = 3)
    )
  
  # Create output filename for stats
  stats_output <- gsub("StatsCombined", "CorrelationStats_withFDR", file)
  write.csv(summary_df, file = stats_output, row.names = FALSE)
  
  # Create plots
  plot_list <- list()
  for (gene in gene_columns) {
    for (metric in behavioral_metrics) {
      # Get correlation stats
      stats <- correlation_results %>%
        filter(Gene == gene & Metric == metric)
      
      # Create plot
      p <- create_gene_scatter(df, gene, metric, cell_type, task, stats)
      plot_list[[paste(gene, metric)]] <- p
    }
  }
  
  # Arrange and save plots
  final_plot <- wrap_plots(plot_list, ncol = 4) +
    plot_annotation(
      title = paste(cell_type, "Gene Z-Scores vs", task, "Behavioral Metrics"),
      subtitle = "Points with black borders and * indicate FDR < 0.05"
    )
  
  # Create output filename for plots
  plot_output <- gsub("StatsCombined", "scatterplots_withFDR", file)
  plot_output <- gsub(".csv", ".pdf", plot_output)
  
  ggsave(plot_output, 
         final_plot,
         width = 20, 
         height = 5 * length(gene_columns),
         limitsize = FALSE)
  
  # Print progress
  cat(sprintf("Completed analysis for %s\n", file))
}
