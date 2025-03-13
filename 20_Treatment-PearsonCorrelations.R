

# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_Cluster_zscore-Stats")


# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)  # for plot arrangement


# Read the combined data
df <- read.csv("C3_TaskD-149_StatsCombined_2025-02-03.csv")

# Gather columns for behavioral metrics
behavioral_metrics <- c("Correct_Response", "Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")

# Function to create treatment-specific correlation analysis
analyze_treatment_correlations <- function(data, gene, behavioral_metric) {
  # Create a list to store results
  treatment_correlations <- list()
  
  # Get unique treatment groups
  treatments <- unique(data$Treatment)
  
  # Calculate correlations for each treatment group
  for (treatment in treatments) {
    # Filter data for specific treatment
    treatment_data <- data[data$Treatment == treatment, ]
    
    # Calculate correlation
    correlation_result <- tryCatch({
      cor.test(treatment_data[[gene]], treatment_data[[behavioral_metric]])
    }, error = function(e) {
      list(estimate = NA, p.value = NA, conf.int = c(NA, NA))
    })
    
    # Store results
    treatment_correlations[[treatment]] <- list(
      correlation = correlation_result$estimate,
      p_value = correlation_result$p.value,
      ci_lower = correlation_result$conf.int[1],
      ci_upper = correlation_result$conf.int[2]
    )
  }
  
  return(treatment_correlations)
}

# Get all gene columns (excluding non-gene columns)
gene_columns <- names(df)[!names(df) %in% 
                            c("Sample", "Treatment", "Correct_Response", 
                              "Premature_Response", "Missed_Response_Window", 
                              "Wrong_Choice")]

# Create comprehensive correlation analysis
comprehensive_correlation_results <- list()

# Perform analysis for each gene and behavioral metric
for (gene in gene_columns) {
  gene_results <- list()
  for (metric in behavioral_metrics) {
    # Analyze correlations
    treatment_specific_correlations <- analyze_treatment_correlations(df, gene, metric)
    
    # Prepare results for each treatment
    treatment_summary <- lapply(names(treatment_specific_correlations), function(treatment) {
      result <- treatment_specific_correlations[[treatment]]
      data.frame(
        Gene = gene,
        Behavioral_Metric = metric,
        Treatment = treatment,
        Correlation = result$correlation,
        P_Value = result$p_value,
        CI_Lower = result$ci_lower,
        CI_Upper = result$ci_upper
      )
    })
    
    # Combine results
    gene_results[[metric]] <- do.call(rbind, treatment_summary)
  }
  
  # Combine all metrics for this gene
  comprehensive_correlation_results[[gene]] <- do.call(rbind, gene_results)
}

# Combine all results into a single dataframe
final_correlation_results <- do.call(rbind, comprehensive_correlation_results)

# Save results to CSV
write.csv(final_correlation_results, 
          file = "C3_TaskD-149_TreatmentSpecific_CorrelationStats_2025-02-03.csv", 
          row.names = FALSE)

# Create treatment-specific scatter plots
create_treatment_scatter <- function(data, gene, behavioral_metric) {
  # Prepare the data
  plot_data <- data %>%
    select(all_of(c("Treatment", behavioral_metric, gene))) %>%
    drop_na()
  
  # Create the plot
  p <- ggplot(plot_data, aes_string(x = gene, y = behavioral_metric, color = "Treatment", shape = "Treatment")) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = Treatment)) +
    scale_color_viridis_d() +
    scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +
    theme_minimal() +
    labs(
      title = paste("Treatment-Specific:", gene, "vs", behavioral_metric),
      x = paste(gene, "Z-Score"),
      y = behavioral_metric
    )
  
  return(p)
}

# Create and save treatment-specific plots
plot_list <- list()
for (gene in gene_columns) {
  for (metric in behavioral_metrics) {
    plot_list[[paste(gene, metric)]] <- create_treatment_scatter(df, gene, metric)
  }
}

# Arrange plots
library(patchwork)
final_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(title = "Treatment-Specific: C3 Gene Z-Scores vs Task D-149 Behavioral Metrics")

# Save the plot
ggsave("C3-zscore_TaskD-149_TreatmentSpecific_scatterplots_2025-02-03.pdf", 
       final_plot, 
       width = 20, 
       height = 5 * length(gene_columns), 
       limitsize = FALSE)

# Print summary for quick reference
print(final_correlation_results)
