

# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_Stats")

library(tidyr)
library(dplyr)


## Run for each cluster and "AllStats" behavioral study 
# Behavior stats: TaskA-45_AllStats.csv, B-59, C-48, D-149
# Zscores: gaba_byTx_zscores_2025-01-16.csv, glut, oligo, endo, astro, microglia



# Read spreadsheets into data frames
df1 = subset(read.csv("/Volumes/DataBox/MCS2023/Stats/Behavior_StatsFiles/TaskD-149_AllStats.csv"), 
             select = -c(Response_Trials))
df2 = read.csv("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_Stats/astro_byTx_zscores_2025-01-16.csv")

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
combined_df = merge(df1, df2_wide, by = c("Treatment"), all.x = TRUE)

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "astro_TaskD-149_StatsCombined_2025-02-03.csv", row.names = FALSE)



##########################################################



library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

# Read the combined data
df <- read.csv("astro_TaskD-149_StatsCombined_2025-02-03.csv")

# Gather columns for behavioral metrics
behavioral_metrics <- c("Correct_Response", "Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")

# Function to create scatter plot for a specific gene and behavioral metric
create_gene_scatter <- function(data, gene, behavioral_metric) {
  # Prepare the data
  plot_data <- data %>%
    # Select only columns needed
    select(all_of(c("Treatment", behavioral_metric, gene))) %>%
    # Remove rows with NA
    drop_na()
  
  # Calculate correlation
  correlation <- tryCatch({
    cor.test(plot_data[[gene]], plot_data[[behavioral_metric]])
  }, error = function(e) list(estimate = NA, p.value = NA))
  
  # Create the plot
  p <- ggplot(plot_data, aes_string(x = gene, y = behavioral_metric, color = "Treatment", shape = "Treatment")) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis_d() +  # Viridis color palette for discrete data
    scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +  # Set shapes manually
    theme_minimal() +
    labs(
      title = paste("Astrocyte", gene, "vs", "TaskD-149", behavioral_metric),
      x = paste(gene, "Z-Score"),
      y = behavioral_metric,
      subtitle = sprintf("Pearson's R: %.3f (p = %.3f)", 
                         correlation$estimate, 
                         correlation$p.value)
    ) +
    theme(plot.subtitle = element_text(size = 10))
  
  return(p)
}

# Get all gene columns (excluding non-gene columns)
gene_columns <- names(df)[!names(df) %in% 
                            c("Sample", "Treatment", "Correct_Response", 
                              "Premature_Response", "Missed_Response_Window", 
                              "Wrong_Choice")]

# Create plots for all genes and behavioral metrics
plot_list <- list()
for (gene in gene_columns) {
  for (metric in behavioral_metrics) {
    plot_list[[paste(gene, metric)]] <- create_gene_scatter(df, gene, metric)
  }
}

# Arrange plots
library(patchwork)
final_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(title = "Astrocyte Gene Z-Scores vs Task D-149 Behavioral Metrics")

# Save the plot
ggsave("astro-zscore_TaskD-149_scatterplots_2025-02-03.pdf", final_plot, 
       width = 20, height = 5 * length(gene_columns), 
       limitsize = FALSE)

# Optionally, print correlation summary
correlation_summary <- lapply(gene_columns, function(gene) {
  sapply(behavioral_metrics, function(metric) {
    cor.test(df[[gene]], df[[metric]])$estimate
  })
})
names(correlation_summary) <- gene_columns
print(correlation_summary)

write.csv(correlation_summary, file = "astro_TaskD-149_CorrelationStats_byTreatment.csv")


######################


# Optionally, print correlation summary with p-values
correlation_summary <- lapply(gene_columns, function(gene) {
  sapply(behavioral_metrics, function(metric) {
    test <- cor.test(df[[gene]], df[[metric]])
    c(Correlation = test$estimate, P_value = test$p.value)
  })
})

names(correlation_summary) <- gene_columns

# Convert to data frame for easier visualization
correlation_summary_df <- do.call(rbind, lapply(correlation_summary, as.data.frame))
correlation_summary_df <- tibble::rownames_to_column(correlation_summary_df, "Gene")

# Print summary
print(correlation_summary_df)

# Save as CSV
write.csv(correlation_summary_df, file = "glut_TaskA-45_CorrelationStats-pValues_byTreatment.csv", row.names = FALSE)


######################


# General Interpretation of Pearson R:
Positive correlation (R > 0): Higher gene expression is associated with an increase in the behavioral metric.
Negative correlation (R < 0): Higher gene expression is associated with a decrease in the behavioral metric.
Strength of correlation:
  0.1 to 0.3 → Weak correlation
0.3 to 0.5 → Moderate correlation
> 0.5 → Strong correlation
< -0.1 follows the same pattern but in the opposite direction.



##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################




# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_Stats")

library(tidyr)
library(dplyr)


## Run for each cluster and "WrongStats" behavioral study 
# Behavior stats: TaskA-45_WrongStats.csv, B-59, C-48, D-149
# Zscores: gaba_byTx_zscores_2025-01-16.csv, glut, oligo, endo, astro, microglia



# Read spreadsheets into data frames
df1 = subset(read.csv("/Volumes/DataBox/MCS2023/Stats/Behavior_StatsFiles/TaskC-48_WrongStats.csv"))
df2 = read.csv("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_Stats/astro_byTx_zscores_2025-01-16.csv")

# Rename columns in first dataframe
names(df1) = c("Sample", "Treatment", "Premature_Response", 
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
combined_df = merge(df1, df2_wide, by = c("Treatment"), all.x = TRUE)

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "astro_TaskC-48_WrongStatsCombined_2025-02-03.csv", row.names = FALSE)



##########################################################



library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

# Read the combined data
df <- read.csv("astro_TaskC-48_WrongStatsCombined_2025-02-03.csv")

# Gather columns for behavioral metrics
behavioral_metrics <- c("Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")

# Function to create scatter plot for a specific gene and behavioral metric
create_gene_scatter <- function(data, gene, behavioral_metric) {
  # Prepare the data
  plot_data <- data %>%
    # Select only columns needed
    select(all_of(c("Treatment", behavioral_metric, gene))) %>%
    # Remove rows with NA
    drop_na()
  
  # Calculate correlation
  correlation <- tryCatch({
    cor.test(plot_data[[gene]], plot_data[[behavioral_metric]])
  }, error = function(e) list(estimate = NA, p.value = NA))
  
  # Create the plot
  p <- ggplot(plot_data, aes_string(x = gene, y = behavioral_metric, color = "Treatment", shape = "Treatment")) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis_d() +  # Viridis color palette for discrete data
    scale_shape_manual(values = c("t1" = 16, "t2" = 17, "t3" = 15, "t4" = 18)) +  # Set shapes manually
    theme_minimal() +
    labs(
      title = paste("Astrocyte", gene, "vs", "TaskC-48 Wrong", behavioral_metric),
      x = paste(gene, "Z-Score"),
      y = behavioral_metric,
      subtitle = sprintf("Pearson's R: %.3f (p = %.3f)", 
                         correlation$estimate, 
                         correlation$p.value)
    ) +
    theme(plot.subtitle = element_text(size = 10))
  
  return(p)
}

# Get all gene columns (excluding non-gene columns)
gene_columns <- names(df)[!names(df) %in% 
                            c("Sample", "Treatment",
                              "Premature_Response", "Missed_Response_Window", 
                              "Wrong_Choice")]

# Create plots for all genes and behavioral metrics
plot_list <- list()
for (gene in gene_columns) {
  for (metric in behavioral_metrics) {
    plot_list[[paste(gene, metric)]] <- create_gene_scatter(df, gene, metric)
  }
}

# Arrange plots
library(patchwork)
final_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(title = "Astrocyte Gene Z-Scores vs Task C-48 Wrong Stats Behavioral Metrics")

# Save the plot
ggsave("astro-zscore_TaskC-48-Wrong_scatterplots_2025-02-03.pdf", final_plot, 
       width = 20, height = 5 * length(gene_columns), 
       limitsize = FALSE)

# Optionally, print correlation summary
correlation_summary <- lapply(gene_columns, function(gene) {
  sapply(behavioral_metrics, function(metric) {
    cor.test(df[[gene]], df[[metric]])$estimate
  })
})
names(correlation_summary) <- gene_columns
print(correlation_summary)

write.csv(correlation_summary, file = "astro_TaskC-48-Wrong_CorrelationStats_byTreatment.csv")




##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################



## Search through Pearson R correlation results for values abs(.0.4)


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_zscore-Stats")

# Read the file 
data <- read.csv("endo_TaskC-48-Wrong_CorrelationStats_byTreatment.csv", stringsAsFactors = FALSE)

# Identify columns with gene names (from 2nd column onward)
gene_cols <- names(data)[2:ncol(data)]

# Create a logical vector for rows with any gene column having abs value > 0.4
high_value_mask <- apply(data[, gene_cols], 1, function(x) any(abs(x) > 0.4))

# Select rows based on the mask
high_value_rows <- data[high_value_mask, c("X", gene_cols)]

# Write the results to a CSV file
write.csv(high_value_rows, "high_value_genes.csv", row.names = FALSE)




##########################################################
##########################################################
##########################################################
##########################################################
##########################################################



## Loop through all files in the folder:


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_zscore-Stats")


# Get list of files ending with "CorrelationStats_byTreatment.csv"
files <- list.files(pattern = "CorrelationStats_byTreatment.csv$")

# Loop through each file
for (file in files) {
  # Read the file 
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract file identifier (everything before "CorrelationStats_byTreatment.csv")
  file_identifier <- sub("_CorrelationStats_byTreatment\\.csv$", "", file)
  
  # Identify columns with gene names (from 2nd column onward)
  gene_cols <- names(data)[2:ncol(data)]
  
  # Create a logical vector for rows with any gene column having abs value > 0.4
  high_value_mask <- apply(data[, gene_cols], 1, function(x) any(abs(x) > 0.4))
  
  # Select rows based on the mask
  high_value_rows <- data[high_value_mask, c("X", gene_cols)]
  
  # Create output filename
  output_filename <- paste0(file_identifier, "_high_value_genes.csv")
  
  # Write the results to a CSV file
  write.csv(high_value_rows, output_filename, row.names = FALSE)
  
  # Optional: Print a message to confirm processing
  cat("Processed", file, "- Output saved as", output_filename, "\n")
}





##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################




## Tried to save all results abs(>0.4) into 1 file, but it didn't save properly:



# Load required library
library(dplyr)
library(readr)

# Get list of files ending with "StatsCombined_2025-02-03.csv"
files <- list.files(pattern = "StatsCombined_2025-02-03.csv$")

# Initialize summary dataframe
summary_df <- tibble(
  File = character(),
  CellType = character(),
  Task = character(),
  HighValueRowsCount = numeric()
)

# Process files and combine results
all_results <- lapply(files, function(file) {
  # Read the file 
  data <- read.csv(file, stringsAsFactors = FALSE)
  
  # Extract file identifier parts
  file_parts <- strsplit(sub("_StatsCombined_2025-02-03\\.csv$", "", file), "_")[[1]]
  
  # Determine cell type and task
  cell_type <- file_parts[1]  # First part of filename
  task <- file_parts[length(file_parts)]  # Last part of filename
  
  # Identify columns with gene names (from 6th column onward)
  gene_cols <- names(data)[6:ncol(data)]
  
  # Create a logical vector for rows with any gene column having abs value > 0.4
  high_value_mask <- apply(data[, gene_cols], 1, function(x) any(abs(x) > 0.4))
  
  # Select rows based on the mask
  high_value_rows <- data[high_value_mask, c("Treatment", "Sample", gene_cols)]
  
  # Add cell type and task columns
  high_value_rows$CellType <- cell_type
  high_value_rows$Task <- task
  
  # Create summary entry
  summary_entry <- tibble(
    File = file,
    CellType = cell_type,
    Task = task,
    HighValueRowsCount = nrow(high_value_rows)
  )
  
  # Attach summary to global summary_df
  summary_df <<- bind_rows(summary_df, summary_entry)
  
  return(high_value_rows)
}) %>% 
  bind_rows()

# Write consolidated results
write_csv(all_results, "consolidated_high_value_genes.csv")

# Write summary
write_csv(summary_df, "processed_files_summary.csv")

# Optional: Print summary to console
print(summary_df)
























##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
