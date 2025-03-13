
# Load necessary libraries
library(tidyverse)

# Total cell counts per cluster
total_cells <- c(88, 161, 396, 588, 6348, 1758, 3478, 2929, 919, 4351, 1026, 
                 7411, 4907, 1064, 1611, 2925, 1029, 3392)

# Create the data frame for DA gene comparisons
data <- data.frame(
  Gene = as.character(c("Total", "Total", "Total")),
  Comparison = as.character(c("2N/2N+ vs. Ts/Ts+", "2N+/Ts+ vs. 2N/Ts", "Ts vs. Ts+")),
  C1 = c(51, 34, 75),
  C2 = c(57, 51, 75),
  C3 = c(49, 61, 90),
  C4 = c(66, 51, 101),
  C5 = c(68, 41, 93),
  C6 = c(65, 48, 104),
  C7 = c(19, 24, 40),
  C8 = c(73, 41, 110),
  C9 = c(57, 47, 101),
  C10 = c(74, 42, 100),
  C11 = c(14, 22, 28),
  C12 = c(41, 40, 60),
  C13 = c(15, 8, 20),
  C14 = c(41, 44, 64),
  C15 = c(18, 14, 21),
  C16 = c(64, 48, 71),
  C17 = c(43, 58, 104),
  C18 = c(66, 47, 97),
  stringsAsFactors = FALSE
)

# Reshape the comparison data to long format
data_long <- data %>%
  pivot_longer(cols = C1:C18, names_to = "Cluster", values_to = "Count")

# Create a data frame for the total cell counts, aligning clusters
total_cells_data <- data.frame(
  Gene = "Total",
  Comparison = "Total Cells",
  Cluster = paste0("C", 1:18),
  Count = total_cells
)

# Combine the comparison and total cell count data
data_combined <- bind_rows(data_long, total_cells_data)

# Ensure proper cluster ordering
data_combined$Cluster <- factor(data_combined$Cluster,
                                levels = paste0("C", 1:18))

# Convert Comparison to a factor to control legend order
data_combined$Comparison <- factor(data_combined$Comparison,
                                   levels = c("2N/2N+ vs. Ts/Ts+", "2N+/Ts+ vs. 2N/Ts",
                                              "Ts vs. Ts+", "Total Cells"))

# Define custom colors for the comparisons and total cells
colors <- c("2N/2N+ vs. Ts/Ts+" = "#440154",
            "2N+/Ts+ vs. 2N/Ts" = "#31688EFF",
            "Ts vs. Ts+" = "#35B779FF",
            "Total Cells" = "#FDE725FF")

# Create new positions for the bars
data_combined <- data_combined %>%
  mutate(x_position = as.numeric(Cluster) +
           case_when(
             Comparison == "2N/2N+ vs. Ts/Ts+" ~ -0.25,  # Adjusted position
             Comparison == "2N+/Ts+ vs. 2N/Ts" ~ 0,
             Comparison == "Ts vs. Ts+" ~ 0.25,  # Adjusted position
             Comparison == "Total Cells" ~ 0.5   # Adjusted position
           ))

# Create the plot
ggplot(data_combined, aes(fill = Comparison)) +
  geom_bar(data = subset(data_combined, Comparison != "Total Cells"),
           aes(x = x_position, y = Count),
           stat = "identity",
           width = 0.28) +
  geom_bar(data = subset(data_combined, Comparison == "Total Cells"),
           aes(x = x_position,
               y = Count / max(total_cells) * max(data_long$Count)),
           stat = "identity",
           width = 0.15) +
  scale_x_continuous(
    breaks = 1:25,
    labels = paste0("C", 1:25)
  ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(
    name = "Number of DA Genes",
    sec.axis = sec_axis(~ . * max(total_cells) / max(data_long$Count),
                        name = "Total Cell Count")
  ) +
  labs(
    title = "C18 Subset DA Genes and Total Cell Counts by Cluster",
    x = "Cluster",
    fill = "Comparison"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

