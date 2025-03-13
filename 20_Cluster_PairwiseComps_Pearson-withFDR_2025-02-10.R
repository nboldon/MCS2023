


# Set working directory
setwd("/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps")


## Combine spreadsheets for downstream analysis


## Run for each cluster and "AllStats" behavioral study 
# Behavior stats: TaskA-45_AllStats.csv, B-59, C-48, D-149
# Zscores: C3_byTx_zscores_2024-10-18.csv, C8, C11, C14, C18, C22



# Read spreadsheets into data frames
df1 = subset(read.csv("/Volumes/DataBox/MCS2023/Stats/Behavior_StatsFiles/TaskD-149_AllStats.csv"), 
             select = -c(Response_Trials))
df2 = read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/Tx_Comps/C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

# Rename columns in first dataframe
names(df1) = c("Sample", "Treatment", "Correct_Response", "Premature_Response", 
               "Missed_Response_Window", "Wrong_Choice")

# Reshape the second dataframe to long format
df2_long <- df2 %>%
  pivot_longer(cols = c(t1, t2, t3, t4), 
               names_to = "Treatment", 
               values_to = "Gene_Score") %>%
  rename(Gene = X)

# Spread the data to create columns for each gene
df2_wide <- df2_long %>%
  pivot_wider(names_from = Gene, 
              values_from = Gene_Score)

# Merge the dataframes
combined_df = merge(df1, df2, all.x = TRUE)

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C22_TaskD-149_StatsCombined_2025-02-03.csv", row.names = FALSE)



##########################################################




library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(patchwork)


## Complete for clusters: C3, C8, C11, C14, C18, C22

## Complete for each behavioral test: TaskA-45, TaskB-59, TaskC-48, TaskD-149



# Read the combined data
df <- read.csv("C3_TaskD-149_StatsCombined_2025-02-03.csv")

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
      title = paste("C3", gene, "vs", "TaskD-149", behavioral_metric),
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
          file = "C3_TaskD-149_CorrelationStats_withFDR_2025-02-12.csv", 
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
        title = paste("C3", gene, "vs", "TaskD-149", metric),
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
  plot_annotation(title = "C3 Gene Z-Scores vs Task D-149 Behavioral Metrics")

ggsave("C3-zscore_TaskD-149_scatterplots_withFDR_2025-02-12.pdf", 
       final_plot,
       width = 20, 
       height = 5 * length(gene_columns),
       limitsize = FALSE)




sessionInfo()


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

