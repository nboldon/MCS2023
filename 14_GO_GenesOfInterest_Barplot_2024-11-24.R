
## Creates barplots showing the number of GO genes of interest per cluster, by pairwise comparison

setwd("/Volumes/DataBox/GO_Analysis")


# Load necessary libraries
library(tidyverse)
library(gtools)  # For mixedsort
library(viridis)

# File path
file_path <- "GO_Genes-Of-Interest_t1-t2_t3-t4_2024-11-08.csv"
file_path <- "GO_Genes-Of-Interest_t2-t4_t1-t3_2024-11-08.csv"
file_path <- "GO_Genes-Of-Interest_t3_t4_2024-11-08.csv"

# Step 1: Load the data correctly
# Read the data and handle missing values
data <- read.delim(file_path, sep = "\t", header = TRUE, na.strings = c("NULL", "", " "))

# Step 2: Check if the data is in one concatenated column, then split
if (ncol(data) == 1) {
  # Split the single column into individual cluster columns, assuming each cluster is separated by commas
  data <- data %>%
    separate(col = 1, into = paste0("C", 1:25), sep = ",", fill = "right", convert = TRUE)
}

# Step 3: Replace "NULL" values with NA so they are treated as missing data
data <- data %>%
  mutate(across(everything(), ~ na_if(., "NULL")))

# Step 4: Clean and count the number of valid genes per cluster column
# Remove rows with NA, empty, or "NULL" values, then count genes per cluster
gene_counts <- data %>%
  # Remove rows where all values are NA, empty, or "NULL"
  filter(rowSums(!is.na(.) & . != "" & . != " ") > 0) %>%
  # Count the number of non-empty genes for each cluster
  summarise(across(everything(), ~ sum(!is.na(.) & . != "" & . != " "))) %>%
  pivot_longer(cols = everything(), names_to = "Cluster", values_to = "Count")

# Step 5: Ensure clusters are in the correct order (C1, C2, C3, ...)
gene_counts <- gene_counts %>%
  mutate(Cluster = factor(Cluster, levels = mixedsort(unique(Cluster))))

# Step 6: Create the bar plot
plot <- ggplot(gene_counts, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(
    title = "Number of GO Genes of Interest in Each Cluster - Ts vs. Ts+ Pairwise Comps",
    x = "Cluster",
    y = "Number of GO Genes of Interest"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Step 7: Save the plot as a file
output_file <- "GO_GenesOfInterest-byCluster_t3_t4_BarPlot_2024-11-24.png"
ggsave(output_file, plot, width = 10, height = 6, dpi = 300)

# Print the plot to view it in the console
print(plot)
