
setwd("/Volumes/DataBox/Heatmap_Comps")


#################################################################

## Rainbow color pallete 

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the Excel file
file_path <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_TxComp_GeneUnions_2024-07-24.xlsx"
data <- read_excel(file_path, sheet = 1)

# Select only comparison columns, excluding any columns that start with "Log2FC" or "FDR"
data_filtered <- data %>%
  select(-starts_with("Log2FC"), -starts_with("FDR"))

# Count unique genes in each comparison, filtering out NULL or empty values
gene_counts <- data_filtered %>%
  summarise_all(~ sum(!is.na(.) & . != "NULL" & . != "")) %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Gene_Count")

# Create the plot
gene_plot <- ggplot(gene_counts, aes(x = Comparison, y = Gene_Count, fill = Comparison)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Counts by Comparison",
       x = "Comparison Group",
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 10, 10, 10))

# Display the plot
print(gene_plot)

# Save the plot as a JPG file
ggsave("gene_counts_by_comparison.jpg", plot = gene_plot, width = 10, height = 6, dpi = 300)




#####################################################

## Color blind-friendly palette


# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

# Load the Excel file
file_path <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_TxComp_GeneUnions_2024-07-24.xlsx"
data <- read_excel(file_path, sheet = 1)

# Select only comparison columns, excluding any columns that start with "Log2FC" or "FDR"
data_filtered <- data %>%
  select(-starts_with("Log2FC"), -starts_with("FDR"))

# Count unique genes in each comparison, filtering out NULL or empty values
gene_counts <- data_filtered %>%
  summarise_all(~ sum(!is.na(.) & . != "NULL" & . != "")) %>%
  pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Gene_Count")

# Create the plot using viridis color palette
gene_plot <- ggplot(gene_counts, aes(x = Comparison, y = Gene_Count, fill = Comparison)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "D") +  # Colorblind-friendly palette
  theme_minimal() +
  labs(title = "Gene Counts by Comparison",
       x = "Comparison Group",
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(10, 10, 10, 10))

# Display the plot
print(gene_plot)

# Save the plot as a JPG file
ggsave("gene_counts_by_comparison.jpg", plot = gene_plot, width = 10, height = 6, dpi = 300)



##########################################################
##########################################################
##########################################################


## Loop for all clusters 


setwd("/Volumes/DataBox/Heatmap_Comps")

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

# Define the path and file naming convention
base_path <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/"
file_template <- "C%d_TxComp_GeneUnions_2024-07-24.xlsx"

# Loop through clusters 1 to 25
for (cluster in 1:25) {
  
  # Construct the file path for the current cluster
  file_path <- sprintf(paste0(base_path, file_template), cluster)
  
  # Check if the file exists
  if (file.exists(file_path)) {
    
    # Load the Excel file
    data <- read_excel(file_path, sheet = 1)
    
    # Select only comparison columns, excluding any columns that start with "Log2FC" or "FDR"
    data_filtered <- data %>%
      select(-starts_with("Log2FC"), -starts_with("FDR"))
    
    # Count unique genes in each comparison, filtering out NULL or empty values
    gene_counts <- data_filtered %>%
      summarise_all(~ sum(!is.na(.) & . != "NULL" & . != "")) %>%
      pivot_longer(cols = everything(), names_to = "Comparison", values_to = "Gene_Count")
    
    # Create the plot using the viridis color palette
    gene_plot <- ggplot(gene_counts, aes(x = Comparison, y = Gene_Count, fill = Comparison)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      theme_minimal() +
      labs(title = paste("Gene Counts by Pairwise Comparison - Cluster", cluster),
           x = "Comparison Group",
           y = "Number of Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.margin = margin(10, 10, 10, 10))
    
    # Display the plot
    print(gene_plot)
    
    # Save the plot as a JPG file with a unique name for each cluster
    output_file <- sprintf("C%d_gene_counts_by_comps_2024-11-08.jpg", cluster)
    ggsave(output_file, plot = gene_plot, width = 10, height = 6, dpi = 300)
    
    cat("Processed and saved plot for Cluster", cluster, "\n")
    
  } else {
    # If the file does not exist, print a message and continue to the next cluster
    cat("File for Cluster", cluster, "not found. Skipping...\n")
  }
}


## Files not found for Cluster 1, 2, 5, 7, 9, 19



#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
