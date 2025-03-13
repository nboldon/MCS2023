# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats")

# Load required libraries
library(ggplot2)
library(dplyr)

# Load your data
data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")


################################################################
################################################################
################################################################


## Pearson correlations
## To plot each comparison seperately, includes all 4 treatment groups

# Convert 'Treatment' to a factor if not already
data$Treatment <- as.factor(data$Treatment)

# Define the behavioral and gene variables
behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")
gene_vars <- c("Aard", "Htr4", "Htr7", "B430306N03Rik", "4930546C10Rik", "Kcns1", "Wfdc5", "Scgb2b1", "Kcnj11")

# Create a directory to save the plots
if (!dir.exists("C19_pearson-correlation_plots")) {
  dir.create("C19_pearson-correlation_plots")
}

# Loop through each behavioral variable and gene variable
for (behavior_var in behavior_vars) {
  for (gene_var in gene_vars) {
    
    # Check if gene column exists in the data
    if (!(gene_var %in% colnames(data))) {
      message(paste("Gene variable", gene_var, "not found in the data. Skipping..."))
      next
    }
    
    # Initialize an empty list to store the correlation values
    correlation_results <- data.frame(Treatment = character(), Correlation = numeric())
    
    # Loop through each treatment group and calculate the correlation
    for (treatment in unique(data$Treatment)) {
      subset_data <- data %>% filter(Treatment == treatment)
      
      # Skip if the subset is empty
      if (nrow(subset_data) == 0) next
      
      # Ensure both behavior and gene variables are numeric using mutate
      subset_data <- subset_data %>%
        mutate(
          !!behavior_var := as.numeric(as.character(.data[[behavior_var]])),
          !!gene_var := as.numeric(as.character(.data[[gene_var]]))
        )
      
      # Check for NA values and remove them if necessary
      subset_data <- subset_data %>%
        filter(!is.na(!!sym(behavior_var)), !is.na(!!sym(gene_var)))
      
      # Skip if there is not enough data to compute correlation (e.g., not enough rows after filtering)
      if (nrow(subset_data) < 2) next  # Need at least 2 data points to compute correlation
      
      # Calculate correlation
      corr <- cor(subset_data[[behavior_var]], subset_data[[gene_var]], method = "pearson")
      
      # Store the correlation result
      correlation_results <- rbind(correlation_results, data.frame(Treatment = treatment, Correlation = corr))
    }
    
    # Plot the correlations for the current behavior and gene pair across treatment groups
    plot_file <- paste0("C19_pearson-correlation_plots/", behavior_var, "_vs_", gene_var, "_2024-11-14.png")
    png(plot_file)
    
    # Create scatter plot with regression lines for each treatment
    p <- ggplot(data, aes(x = !!sym(behavior_var), y = !!sym(gene_var), color = Treatment)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = paste("Correlation of", behavior_var, "vs", gene_var),
           x = behavior_var,
           y = gene_var) +
      theme_minimal()
    
    # Add Pearson's correlation coefficient to the plot (only for the first treatment group)
    corr_text <- paste("Pearson Correlation: ", round(correlation_results$Correlation[1], 2), sep = "")
    p <- p + annotate("text", x = Inf, y = -Inf, label = corr_text, hjust = 1.1, vjust = -0.1, size = 4)
    
    # Print the plot
    print(p)
    
    dev.off()  # Close the file device
  }
}


# Gene variable 4930546C10Rik not found in the data.
# Gene variable 4930546C10Rik not found in the data. 
# Gene variable 4930546C10Rik not found in the data
# Gene variable 4930546C10Rik not found in the data.



################################################################
################################################################
################################################################


## Pearson correlations
## Plots individual comparisons on the same page, includes all 4 treatment groups

# Load the required libraries
library(tidyr)
library(ggplot2)

# Check if all gene names exist in the data
missing_genes <- setdiff(gene_vars, colnames(data))
if (length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

# Filter gene_vars to only include genes that exist in the dataset
valid_gene_vars <- intersect(gene_vars, colnames(data))

# Loop through each behavioral variable
for (behavior_var in behavior_vars) {
  # Reshape the data from wide format to long format, using only valid genes
  long_data <- data %>%
    select(Treatment, all_of(behavior_var), all_of(valid_gene_vars)) %>%
    pivot_longer(cols = all_of(valid_gene_vars), 
                 names_to = "Gene", 
                 values_to = "GeneExpression")
  
  # Ensure that the behavior variable is numeric
  long_data[[behavior_var]] <- as.numeric(long_data[[behavior_var]])
  
  # Create a directory to save the plots
  if (!dir.exists("C19_all_genes_correlation_plots")) {
    dir.create("C19_all_genes_correlation_plots")
  }
  
  # Create the plot for the current behavior variable across all genes
  plot_file <- paste0("C19_all_genes_correlation_plots/", behavior_var, "_all_genes_vs_", "2024-11-14.png")
  png(plot_file)
  
  # Create scatter plot for all genes across treatment groups
  p <- ggplot(long_data, aes(x = !!sym(behavior_var), y = GeneExpression, color = Treatment)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(group = Gene)) +  # Add a regression line for each gene
    facet_wrap(~Gene, scales = "free_y") +  # Create separate plots for each gene
    labs(title = paste("Correlation of", behavior_var, "vs Genes by Treatment"),
         x = behavior_var,
         y = "Gene Expression") +
    theme_minimal()
  
  # Print the plot
  print(p)
  
  dev.off()  # Close the file device
}


# Missing genes: 4930546C10Rik 




################################################################
################################################################
################################################################





## Step-by-Step Debugging

# Check Column Names in filtered_data
print(names(filtered_data))

# Check for NA Values in Behavior Columns
print(sapply(filtered_data[behavior_vars], function(x) sum(is.na(x))))

# Check for Any Non-Numeric Values
print(sapply(filtered_data[behavior_vars], function(x) sum(!is.numeric(x))))

# Check for Special Characters or Whitespaces in Column Names
print(names(filtered_data))
# Clean up column names by trimming whitespace
names(filtered_data) <- trimws(names(filtered_data))


# Check long_data Creation
print(head(long_data))

# Check column type after reshaping
long_data <- long_data %>%
  mutate(!!behavior_var := as.numeric(!!sym(behavior_var)))

# Double check dataframe creation
print(head(long_data))


###########################################################################



## Pearson correlations
## T1 vs T3 gene correlations - Task on x-axis, gene accessibility on y-axis

# Step 1: Set working directory and load required libraries
setwd("/Volumes/DataBox/MCS2023/Stats")  
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 2: Load your data
data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")
# Warning message: Incomplete final line is not critical, it's just indicating the file doesn't end with a new line.

# Step 3: Define the variables for behavior and gene expression
behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")  # Adjust as needed
gene_vars <- c("Aard", "Htr4", "Htr7", "B430306N03Rik", "4930546C10Rik", "Kcns1", "Wfdc5", "Scgb2b1", "Kcnj11")  # Adjust as needed

# Step 4: Replace any mismatched gene names (e.g., "4930546C10Rik" should match "X4930546C10Rik")
gene_vars <- gsub("^4930546C10Rik$", "X4930546C10Rik", gene_vars)

# Step 5: Subset the data to include only treatment groups t1 and t3
filtered_data <- data %>% filter(Treatment %in% c("t1", "t3"))  # Select only t1 and t3 groups

# Step 6: Loop through each behavioral variable
for (behavior_var in behavior_vars) {
  
  # Step 7: Reshape the data from wide format to long format
  long_data <- filtered_data %>%
    select(Treatment, all_of(behavior_var), all_of(gene_vars)) %>%
    pivot_longer(cols = all_of(gene_vars), 
                 names_to = "Gene", 
                 values_to = "GeneExpression")
  
  # Step 8: Ensure the behavior variable is numeric
  long_data[[behavior_var]] <- as.numeric(long_data[[behavior_var]])
  
  # Step 9: Remove rows where behavior or gene expression is NA
  long_data <- long_data %>%
    filter(!is.na(!!sym(behavior_var)) & !is.na(GeneExpression))
  
  # Step 10: Calculate Pearson correlation for each gene and treatment group (t1, t3)
  correlation_results <- long_data %>%
    group_by(Treatment, Gene) %>%
    summarise(Correlation = cor(get(behavior_var), GeneExpression, method = "pearson", use = "complete.obs"),
              .groups = "drop") %>%
    filter(Treatment %in% c("t1", "t3"))  # Keep only correlations for t1 and t3
  
  # Step 11: Check for any correlations that are still NA
  correlation_results <- correlation_results %>%
    filter(!is.na(Correlation))
  
  # Step 12: Create a directory to save the plots if it doesn't exist
  if (!dir.exists("C19_t1_t3_gene_correlation_plots")) {
    dir.create("C19_t1_t3_gene_correlation_plots")
  }
  
  # Step 13: Create the plot for the current behavior variable across all genes for t1 and t3
  plot_file <- paste0("C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_t3_gene_correlation.png")
  png(plot_file)
  
  # Step 14: Scatter plot for all genes across treatment groups t1 and t3
  p <- ggplot(long_data, aes(x = !!sym(behavior_var), y = GeneExpression, color = Gene, shape = Treatment)) +
    geom_point(alpha = 0.6) +  # Add points for gene expression
    geom_smooth(method = "lm", se = FALSE, aes(group = interaction(Gene, Treatment)), color = "black") +  # Regression lines
    labs(title = paste("Correlation of", behavior_var, "vs Genes for t1 and t3 Treatment Groups"),
         x = behavior_var,
         y = "Gene Expression") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Step 15: Annotate Pearson's correlation on the plot
  for (i in 1:nrow(correlation_results)) {
    corr_text <- paste0("Pearson Correlation for ", correlation_results$Gene[i], ": ", 
                        round(correlation_results$Correlation[i], 2))
    
    p <- p + 
      annotate("text", x = 0.8, y = 0.2 + (i * 0.1), label = corr_text, size = 3, hjust = 0) # Adjust position as needed
  }
  
  # Step 16: Print the plot
  print(p)
  
  dev.off()  # Close the file device
  
  # Step 17: Optionally save correlation results to a CSV for further analysis
  write.csv(correlation_results, paste0("C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_t3_correlations.csv"))
}


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################


## Scatterplot for t1 on x-axis and t3 on y-axis, comparing gene accessibility scores


# Set working directory and load required libraries
setwd("/Volumes/DataBox/MCS2023/Stats")  # Set your working directory
library(dplyr)
library(tidyr)
library(ggplot2)

# Load your data
data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")

# Convert 'Treatment' to a factor if not already
data$Treatment <- as.factor(data$Treatment)

# Define the variables for behavior and gene expression
behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")  # Adjust as needed
gene_vars <- c("Aard", "Htr4", "Htr7", "B430306N03Rik", "4930546C10Rik", "Kcns1", "Wfdc5", "Scgb2b1", "Kcnj11")  # Adjust as needed

# Replace any mismatched gene names (e.g., "4930546C10Rik" should match "X4930546C10Rik")
gene_vars <- gsub("^4930546C10Rik$", "X4930546C10Rik", gene_vars)

# Subset the data to include only treatment groups t1 and t3
filtered_data <- data %>% filter(Treatment %in% c("t1", "t3"))  # Select only t1 and t3 groups

library(tidyr)

# Assuming 'filtered_data' is your original dataset
long_data <- filtered_data %>%
  pivot_longer(cols = -Treatment,  # Keeping 'Treatment' as an identifier
               names_to = "Gene",  # Column for gene names
               values_to = "Expression")  # Column for expression values

# Create a wide format version of the data (using pivot_wider)
long_data_wide <- long_data %>%
  pivot_wider(names_from = Treatment, values_from = Expression)

# Check the transformed wide format data
head(long_data_wide)

# Continue with the rest of your analysis
# Assuming the long_data_wide is already created
library(ggplot2)

# Loop through each behavioral variable
for (behavior_var in behavior_vars) {
  
  # Create a file name for saving the plot
  plot_file <- paste0("C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_vs_t3_scatter_2024-11-14.png")
  
  # Specify plot size using ggsave for better control over dimensions
  png(plot_file, width = 20, height = 28, units = "in", res = 300)  # You can adjust width, height
  
  # Scatter plot for t1 (x-axis) vs t3 (y-axis) for each gene
  p <- ggplot(long_data_wide, aes(x = t1, y = t3, color = Gene)) +
    geom_point(alpha = 0.6, size = 3) +  # Adjust the point size (increased from default)
    geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add regression line
    labs(title = paste("Scatter plot of Gene Accessibility: 2N vs Ts for", behavior_var),
         x = "2N Gene Accessibility",
         y = "Ts Gene Accessibility") +
    theme_minimal() +
    theme(
      legend.position = "bottom",  # Place legend at the bottom
      legend.title = element_text(size = 12),  # Increase legend title size
      legend.text = element_text(size = 10),  # Increase legend item size
      plot.margin = margin(1, 1, 1, 1, "cm")  # Add margin to prevent clipping
    )
  
  # Print the plot
  print(p)
  
  dev.off()  # Close the file device
  
  # Optionally save correlation results to a CSV
  write.csv(long_data_wide, paste0("C19_t1_t3_gene_correlation_plots/", behavior_var, "_t1_vs_t3_gene_accessibility.csv"))
}



#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################


## Scatterplot for t1 on x-axis and t3 on y-axis, comparing behavioral test scores (1 plot for each comparison)


# Set working directory and load required libraries
setwd("/Volumes/DataBox/MCS2023/Stats")   
library(dplyr)
library(tidyr)
library(ggplot2)

# Load your data
long_data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")

# Convert 'Treatment' to a factor if not already
long_data$Treatment <- as.factor(long_data$Treatment)

# Define the variables for behavior (Test Scores)
behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")  # Adjust as needed

# Subset the data to include only treatment groups t1 and t3
filtered_data <- long_data %>% filter(Treatment %in% c("t1", "t3"))  # Select only t1 and t3 groups

# Select only relevant columns for analysis: Treatment and the behavior columns
long_data_filtered <- long_data %>%
  select(Treatment, all_of(behavior_vars))  # Select Treatment and all behavior variables

# Create a wide format version of the data (using pivot_wider)
long_data_wide <- long_data_filtered %>%
  pivot_wider(names_from = Treatment, values_from = all_of(behavior_vars))  # Pivot to wide format based on treatments

# Check the transformed wide format data
head(long_data_wide)

# Ensure 'Behavior' column is included in wide format data - Update this part
# You don't need 'Behavior' or 'Score' columns, as we are already selecting behavior columns
# directly in the wide format pivoting

# Create the directory if it doesn't exist
if (!dir.exists("C19_t1_t3_behavior_correlation_plots")) {
  dir.create("C19_t1_t3_behavior_correlation_plots")
}

# Loop through each behavioral variable for plotting
for (behavior_var in behavior_vars) {
  
  # Create a file name for saving the plot
  plot_file <- paste0("C19_t1_t3_behavior_correlation_plots/", behavior_var, "_t1_vs_t3_scatter.png")
  
  # Specify plot size using ggsave for better control over dimensions
  png(plot_file, width = 20, height = 28, units = "in", res = 300)  # Adjust width, height
  
  # Remove rows with NA values for the current behavior variable
  long_data_wide_clean <- long_data_wide %>% drop_na()
  
  # Scatter plot for t1 (x-axis) vs t3 (y-axis) for each behavior score
  p <- ggplot(long_data_wide_clean, aes(x = !!sym(paste0(behavior_var, "_t1")), 
                                        y = !!sym(paste0(behavior_var, "_t3")))) +
    geom_point(alpha = 0.6, size = 3) +  # Adjust the point size (increased from default)
    geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add regression line
    labs(title = paste("Scatter plot of Behavioral Test Scores: 2N vs Ts for", behavior_var),
         x = "Behavior Score (2N)",
         y = "Behavior Score (Ts)") + 
    theme_minimal() + 
    theme(
      legend.position = "bottom",  # Place legend at the bottom
      legend.title = element_text(size = 12),  # Increase legend title size
      legend.text = element_text(size = 10),  # Increase legend item size
      plot.margin = margin(1, 1, 1, 1, "cm")  # Add margin to prevent clipping
    )
  
  # Print the plot
  print(p)
  
  # Close the file device
  dev.off()
  
  # Optionally save correlation results to a CSV
  write.csv(long_data_wide_clean, paste0("C19_t1_t3_behavior_correlation_plots/", behavior_var, "_t1_vs_t3_behavior_scores.csv"))
}


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################


## Scatterplot for t1 on x-axis and t3 on y-axis, comparing behavioral test scores (1 plot for all comparisons)


setwd("/Volumes/DataBox/MCS2023/Stats")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)  # Load viridis package for color palettes

# Load your data
long_data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")

# Convert 'Treatment' to a factor if not already
long_data$Treatment <- as.factor(long_data$Treatment)

# Define the variables for behavior (Test Scores)
behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")

# Reshape the data to long format: one row per behavior and treatment group
long_data_long <- long_data %>%
  select(Treatment, all_of(behavior_vars)) %>%
  pivot_longer(cols = all_of(behavior_vars), 
               names_to = "Behavior", 
               values_to = "Score")

# Check the reshaped data
head(long_data_long)

# Separate the data for t1 and t3 treatments
long_data_t1 <- long_data_long %>%
  filter(Treatment == "t1") %>%
  rename(Score_t1 = Score) %>%
  select(Behavior, Score_t1)

long_data_t3 <- long_data_long %>%
  filter(Treatment == "t3") %>%
  rename(Score_t3 = Score) %>%
  select(Behavior, Score_t3)

# Merge t1 and t3 data by Behavior to plot them against each other
merged_data <- merge(long_data_t1, long_data_t3, by = "Behavior")

# Now reshape the data to make a long format where each row represents a behavior with both t1 and t3
long_data_merged <- merged_data %>%
  pivot_longer(cols = c(Score_t1, Score_t3), 
               names_to = "Treatment", 
               values_to = "Score")

# Check the reshaped data again
head(long_data_merged)

# Create the directory if it doesn't exist
if (!dir.exists("C19_t1_t3_behavior_correlation_plots")) {
  dir.create("C19_t1_t3_behavior_correlation_plots")
}

# Create a plot for all behavior comparisons in one plot
plot_file <- "C19_t1_t3_behavior_comparison_2024-11-14.png"

# Specify plot size using ggsave for better control over dimensions
png(plot_file, width = 12, height = 8, units = "in", res = 300)  # Adjust width, height

# Create the plot with Viridis color palette
p <- ggplot(long_data_merged, aes(x = Score, y = Score, color = Behavior, shape = Treatment)) +
  geom_point(size = 4) +  # Points for each behavioral test
  labs(title = "Behavioral Test Scores: t1 vs t3",
       x = "t1 Score (Treatment 1)",
       y = "t3 Score (Treatment 3)") + 
  scale_color_viridis(discrete = TRUE) +  # Use viridis color palette
  scale_shape_manual(values = c(16, 17)) +  # Different shapes for t1 and t3
  theme_minimal() + 
  theme(
    legend.position = "right",  # Position legend to the right
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 10),  # Increase legend item size
    plot.margin = margin(1, 1, 1, 1, "cm")  # Add margin to prevent clipping
  )

# Print the plot
print(p)

# Close the file device
dev.off()

# Optionally save the reshaped data to a CSV (if you need it)
write.csv(long_data_merged, "C19_t1_t3_merged_behavior_data_2024-11-14.csv")



#######################################################
#######################################################
#######################################################
#######################################################
#######################################################



## Scatterplot for t2-t4 on x-axis and t1-t3 on y-axis, comparing gene zscores scores (1 plot for all comparisons)


setwd("/Volumes/DataBox/MCS2023/Stats")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(viridis)  # Load viridis package for color palettes

# Load your data
# long_data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")

# Using the same .csv file as Ken's MiniTab scatterplots
long_data <- read.csv("C19_significant_t2-t4_t1-t3_zscores_2024-10-18.csv")

# Convert 'Treatment' to a factor if not already
long_data$Treatment <- as.factor(long_data$Treatment)

# Define the variables for genes
#behavior_vars <- c("Task_B_Correct", "Task_B_Premature", "Task_B_Missed", "Task_B_Wrong")
#gene_vars <- c("Aard", "Htr4", "Htr7", "B430306N03Rik", "4930546C10Rik", "Kcns1", "Wfdc5", "Scgb2b1", "Kcnj11")

# Using same genes as Ken's MiniTab scatterplots
gene_vars <- c("Slc12a7", "Smim6", "Hist1h1e", "Cyp8b1", "Ush1c", "Lipo4", "Lrrc15", "B430306N03Rik", "Svs5",
  "Ces1f", "St8sia6", "Pigr", "Teddm1a", "Trem3", "Treml2", "Ch25h", "Gal3st2", "Neil2", "Lipo2", "Mir1264",
  "Olfr874", "Gm5166", "Gm9112", "Tnp1", "Scgb2b1", "Clec7a", "Mir3970")
# Did not include: 4921517D22Rik, 1700028P14Rik, Cldn34-ps

# Replace any mismatched gene names (e.g., "4930546C10Rik" should match "X4930546C10Rik")
#gene_vars <- gsub("^4930546C10Rik$", "X4930546C10Rik", gene_vars)

# Reshape the data to long format: one row per behavior and treatment group
long_data_long <- long_data %>%
  select(Treatment, all_of(gene_vars)) %>%
  pivot_longer(cols = all_of(gene_vars), 
               names_to = "Genes", 
               values_to = "Score")

# Check the reshaped data
head(long_data_long)

# Separate the data for t1 and t3 treatments
long_data_t1 <- long_data_long %>%
  filter(Treatment == "t1") %>%
  rename(Score_t1 = Score) %>%
  select(Genes, Score_t1)

long_data_t2 <- long_data_long %>%
  filter(Treatment == "t2") %>%
  rename(Score_t2 = Score) %>%
  select(Genes, Score_t2)

long_data_t3 <- long_data_long %>%
  filter(Treatment == "t3") %>%
  rename(Score_t3 = Score) %>%
  select(Genes, Score_t3)

long_data_t4 <- long_data_long %>%
  filter(Treatment == "t4") %>%
  rename(Score_t4 = Score) %>%
  select(Genes, Score_t4)



# Merge t1 and t3 data by Behavior to plot them against each other
# The merge() function expects two data frames at a time and does not support merging multiple data frames directly with a single call.
#merged_data <- merge(long_data_t1, long_data_t2, long_data_t3, long_data_t4, by = "Genes")

# List of data frames to merge
data_list <- list(long_data_t1, long_data_t2, long_data_t3, long_data_t4)

# Merge them by "Genes"
merged_data <- reduce(data_list, ~ merge(.x, .y, by = "Genes"))

# Now reshape the data to make a long format where each row represents a behavior with both t1 and t3
long_data_merged <- merged_data %>%
  pivot_longer(cols = c(Score_t1, Score_t2, Score_t3, Score_t4), 
               names_to = "Treatment", 
               values_to = "Score")

# Check the reshaped data again
head(long_data_merged)



# Create a plot for all behavior comparisons in one plot
plot_file <- "C19_t2-t4_t1-t3_gene_scatterplot_2024-11-19.png"

# Specify plot size using ggsave for better control over dimensions
png(plot_file, width = 12, height = 8, units = "in", res = 300)  # Adjust width, height

# Create the plot with Viridis color palette
p <- ggplot(long_data_merged, aes(x = Score, y = Score, color = Genes, shape = Treatment)) +
  geom_point(size = 4) +  # Points for each gene test
  labs(title = "Gene Z-Scores: t2-t4 vs t1-t3",
       x = "t2-t4 Score (2N+ & Ts+)",
       y = "t1-t3 Score (2N & Ts)") + 
  scale_color_viridis(discrete = TRUE) +  # Use viridis color palette
  scale_shape_manual(values = c(16, 17, 18, 19)) +  # Different shapes for t1, t2, t3, and t4
  theme_minimal() + 
  theme(
    legend.position = "right",  # Position legend to the right
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 10),  # Increase legend item size
    plot.margin = margin(1, 1, 1, 1, "cm")  # Add margin to prevent clipping
  )

# Print the plot
print(p)

# Close the file device
dev.off()

# Optionally save the reshaped data to a CSV (if you need it)
write.csv(long_data_merged, "C19_t2-t4_t1-t3_scatterplot_gene_data_2024-11-19.csv")



#######################################################


## The following loop breaks up and plots tx group comps, instead of using all tx groups (as above)


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# Load the data
long_data <- read.csv("C19_significant_t2-t4_t1-t3_zscores_2024-10-18.csv", stringsAsFactors = FALSE)

# Convert 'Treatment' to a factor
long_data$Treatment <- as.factor(long_data$Treatment)

# Define the variables for genes
gene_vars <- c("Slc12a7", "Smim6", "Hist1h1e", "Cyp8b1", "Ush1c", "Lipo4", "Lrrc15", "B430306N03Rik", "Svs5",
               "Ces1f", "St8sia6", "Pigr", "Teddm1a", "Trem3", "Treml2", "Ch25h", "Gal3st2", "Neil2", "Lipo2",
               "Mir1264", "Olfr874", "Gm5166", "Gm9112", "Tnp1", "Scgb2b1", "Clec7a", "Mir3970")

# Reshape the data to long format
long_data_long <- long_data %>%
  select(Treatment, all_of(gene_vars)) %>%
  pivot_longer(cols = all_of(gene_vars),
               names_to = "Genes",
               values_to = "Score")

# Create all pairwise treatment comparisons
treatment_pairs <- combn(unique(long_data_long$Treatment), 2, simplify = FALSE)

# Loop through each treatment pair to create scatterplots
for (pair in treatment_pairs) {
  x_group <- pair[1]
  y_group <- pair[2]
  
  # Extract data for the two groups
  data_x <- long_data_long %>% 
    filter(Treatment == x_group) %>%
    rename(Score_x = Score) %>%
    mutate(Treatment_Group = x_group)  # Assign Treatment_Group
  
  data_y <- long_data_long %>% 
    filter(Treatment == y_group) %>%
    rename(Score_y = Score) %>%
    mutate(Treatment_Group = y_group)  # Assign Treatment_Group
  
  # Check if the genes are the same in both datasets before merging
  common_genes <- intersect(data_x$Genes, data_y$Genes)
  
  # Filter both data_x and data_y to include only common genes
  data_x <- data_x %>% filter(Genes %in% common_genes)
  data_y <- data_y %>% filter(Genes %in% common_genes)
  
  # Merge data by 'Genes' and ensure Treatment_Group is included in the merged data
  merged_data <- merge(data_x, data_y, by = "Genes", suffixes = c("_x", "_y"))
  
  # Add Treatment_Group based on the Treatment of x and y
  merged_data$Treatment_Group <- factor(rep(c(x_group, y_group), length.out = nrow(merged_data)))
  
  # Create the scatter plot with Score_x vs Score_y, differentiated by shape for Treatment_Group
  scatter_plot <- ggplot(merged_data, aes(x = Score_x, y = Score_y)) +
    geom_point(aes(color = Genes, shape = Treatment_Group), alpha = 0.7, size = 3) +  # Added shape by Treatment_Group
    labs(title = paste("Scatterplot of", y_group, "vs", x_group),
         x = paste("Score for", x_group),
         y = paste("Score for", y_group)) +
    scale_color_viridis(discrete = TRUE) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save plot to file
  plot_file <- paste0("scatterplot_", x_group, "_vs_", y_group, ".png")
  png(plot_file, width = 12, height = 8, units = "in", res = 300)
  print(scatter_plot)  # Explicitly print the plot
  dev.off()
}



#######################################################


## Streamline the above code to read in only the columns of interest and auto generate the gene_vars list,
## then prints all plots in one pdf file.


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# Load the data
long_data <- read.csv("C19_significant_t2-t4_t1-t3_zscores_2024-10-18.csv", stringsAsFactors = FALSE)

# Extract gene names dynamically by excluding the 'Treatment' column
# We assume the first column is the Treatment and we capture all the remaining columns for genes
gene_vars <- colnames(long_data)[2:ncol(long_data)]  # Skip the first column (Treatment)

# Reshape the data to long format
long_data_long <- long_data %>%
  select(Treatment, all_of(gene_vars)) %>%
  pivot_longer(cols = all_of(gene_vars),
               names_to = "Genes",
               values_to = "Score")

# Create all pairwise treatment comparisons
treatment_pairs <- combn(unique(long_data_long$Treatment), 2, simplify = FALSE)

# Define the PDF file to save all plots
pdf("C19_sig_t2-t4_t1-t3_zscores_scatterplots_2024-11-19.pdf", width = 12, height = 8)

# Loop through each treatment pair to create scatterplots
for (pair in treatment_pairs) {
  x_group <- pair[1]
  y_group <- pair[2]
  
  # Extract data for the two groups
  data_x <- long_data_long %>% 
    filter(Treatment == x_group) %>%
    rename(Score_x = Score) %>%
    mutate(Treatment_Group = x_group)  # Assign Treatment_Group
  
  data_y <- long_data_long %>% 
    filter(Treatment == y_group) %>%
    rename(Score_y = Score) %>%
    mutate(Treatment_Group = y_group)  # Assign Treatment_Group
  
  # Check if the genes are the same in both datasets before merging
  common_genes <- intersect(data_x$Genes, data_y$Genes)
  
  # Filter both data_x and data_y to include only common genes
  data_x <- data_x %>% filter(Genes %in% common_genes)
  data_y <- data_y %>% filter(Genes %in% common_genes)
  
  # Merge data by 'Genes' and ensure Treatment_Group is included in the merged data
  merged_data <- merge(data_x, data_y, by = "Genes", suffixes = c("_x", "_y"))
  
  # Add Treatment_Group based on the Treatment of x and y
  merged_data$Treatment_Group <- factor(rep(c(x_group, y_group), length.out = nrow(merged_data)))
  
  # Create the scatter plot with Score_x vs Score_y, differentiated by shape for Treatment_Group
  scatter_plot <- ggplot(merged_data, aes(x = Score_x, y = Score_y)) +
    geom_point(aes(color = Genes, shape = Treatment_Group), alpha = 0.7, size = 3) +  # Added shape by Treatment_Group
    labs(title = paste("Scatterplot of", y_group, "vs", x_group),
         x = paste("Score for", x_group),
         y = paste("Score for", y_group)) +
    scale_color_viridis(discrete = TRUE) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Print the plot to the PDF
  print(scatter_plot)
}

# Close the PDF device
dev.off()


#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################
#######################################################


## Pearson correlation values between "Task_B_Correct" and gene expression values across treatments (t1 and t3), 
## with t1 on the x-axis and t3 on the y-axis. 

###### THIS DOES NOT RUN

# Load libraries
library(dplyr)
library(ggplot2)
library(viridis)

# Load your data
data <- read.csv("C19_byTx_zscores-behavior_2024-11-14.csv")

# Filter for relevant treatments (t1 and t3)
filtered_data <- data %>%
  filter(Treatment %in% c("t1", "t3"))

# Separate data for each treatment
data_t1 <- filtered_data %>% filter(Treatment == "t1")
data_t3 <- filtered_data %>% filter(Treatment == "t3")

# Calculate Pearson correlation for each gene with Task_B_Correct for both t1 and t3
correlations <- data.frame(Gene = colnames(data_t1)[5:ncol(data_t1)])  # assuming genes start from the 5th column
correlations$t1_correlation <- sapply(correlations$Gene, function(gene) {
  cor(data_t1[[gene]], data_t1$Task_B_Correct, use = "complete.obs")
})
correlations$t3_correlation <- sapply(correlations$Gene, function(gene) {
  cor(data_t3[[gene]], data_t3$Task_B_Correct, use = "complete.obs")
})

# Plot t1 vs t3 Pearson correlations
plot <- ggplot(correlations, aes(x = t1_correlation, y = t3_correlation, color = Gene)) +
  geom_point(size = 3) +
  labs(title = "Pearson Correlation of Genes with Task_B_Correct (t1 vs t3)",
       x = "Correlation with Task_B_Correct (t1)",
       y = "Correlation with Task_B_Correct (t3)") +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Display the plot
print(plot)

# Save the plot
ggsave("Gene_Correlation_Task_B_Correct_t1_vs_t3.png", plot = plot, width = 28, height = 49, dpi = 300)

