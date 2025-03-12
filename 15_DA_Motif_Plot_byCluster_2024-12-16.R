
# FDR.<=.0.1


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)  # Ensure viridis package is loaded

# Load your data
motifs_data <- read.csv("motifs.csv")

# Rename columns for easier reference
colnames(motifs_data) <- gsub(" ", ".", colnames(motifs_data)) # Replaces spaces with periods

# Reshape data from wide to long format for ggplot
motifs_long <- melt(motifs_data, 
                    id.vars = c("Cluster", "Cell.Type.ID", "Cell.Count"), 
                    measure.vars = c("T1.vs.T2.DA.Motif.Count", 
                                     "T1.vs.T3.DA.Motif..Count", 
                                     "T1.vs.T4.DA.Motif..Count", 
                                     "T2.vs.T1.DA.Motif..Count", 
                                     "T2.vs.T3.DA.Motif..Count"), 
                    variable.name = "Comparison", 
                    value.name = "DA_Motif_Count")

# Ensure the Cluster column is treated as a factor for proper ordering
motifs_long$Cluster <- factor(motifs_long$Cluster, levels = sort(unique(motifs_long$Cluster)))

# Create the plot with viridis color palette
motif_plot <- ggplot(motifs_long, aes(x = Cluster, y = DA_Motif_Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "DA Motif Count Comparison Across Clusters",
       x = "Cluster",
       y = "DA Motif Count",
       fill = "Comparison") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels for better readability
        axis.title.x = element_text(margin = margin(t = 10)),           # Add space above x-axis title
        axis.title.y = element_text(margin = margin(r = 10))) +         # Add space next to y-axis title
  scale_fill_viridis(discrete = TRUE, option = "D")  # Use viridis discrete color palette with a specific option

# Display the plot
print(motif_plot)

# Save the plot to a file (e.g., PNG format)
ggsave("DA_Motif_Count_Comparison_2024-12-16.png", plot = motif_plot, width = 10, height = 6, dpi = 300)

