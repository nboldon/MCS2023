


# Set working directory
setwd("/Volumes/DataBox/Heatmap_Comps")

# Load required libraries
library(viridis)
library(ggplot2)

# Define the file path and pattern
file_path <- "/Volumes/DataBox/ProjMCS6/ResearchQuestions/"
cluster_files <- list.files(file_path, pattern = "GenotypeComp_FDR-0-1_Log2FC-0-5", full.names = TRUE)

# Initialize a data frame to store the counts
cluster_counts <- data.frame(cluster = integer(0), g1_count = integer(0), g2_count = integer(0))

# Loop through each file and count "g1" and "g2" rows
for (file in cluster_files) {
  tryCatch({
    # Read the data from the file
    data <- read.csv(file)
    
    # Count the number of "g1" and "g2" rows
    g1_count <- sum(data$group_name == "g1")
    g2_count <- sum(data$group_name == "g2")
    
    # Extract the cluster number from the file name using regex (assuming cluster is a number after 'C')
    cluster_num <- as.numeric(sub("C(\\d+)_.*", "\\1", basename(file)))
    
    # Append the counts to the data frame
    cluster_counts <- rbind(cluster_counts, data.frame(cluster = cluster_num, g1_count = g1_count, g2_count = g2_count))
  }, error = function(e) {
    # Print an error message and skip to the next file if there's an error
    cat("Error with file", file, ":", e$message, "\n")
  })
}

# Reshape the data for plotting
cluster_counts_long <- reshape(cluster_counts, 
                               varying = c("g1_count", "g2_count"), 
                               v.names = "count", 
                               times = c("g1", "g2"), 
                               timevar = "group", 
                               direction = "long")

# Plot the barplot
p <- ggplot(cluster_counts_long, aes(x = factor(cluster), y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = TRUE, labels = c("2N/2N+", "Ts/Ts+")) +  # Customize the legend labels
  labs(x = "Cluster", 
       y = "Total Gene Count",  # Update y-axis label
       fill = "Group", 
       title = "Total Gene Count by Cluster for Genotype Pairwise Comparisons") +  # Add plot title
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a .jpg file
ggsave("/Volumes/DataBox/Heatmap_Comps/genotyp_gene-counts_byCluster_Barplot_2024-11-08.jpg", plot = p, width = 10, height = 6)


############################################################

