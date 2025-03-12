
# Color Universal Design Pallete: 
# https://jfly.uni-koeln.de/color/

Blue: Hex "#0072B8"
Purple: Hex "#9B59B6"
Yellow: Hex "#F0E442"
Green: Hex "#009E73"

library(tidyverse)
library(colorBlindness)
library(here)

# from https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


colors <- c("#56B4E9", "#F0E442", "#009E73", "#DDA0DD")
#colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
#colors <- c("#FF6F61", "#4CAF50", "#008C9E", "#9B59B6")

######################################
######################################
######################################


# Load necessary library
library(tidyverse)

# Create the data frame again with explicit types
data <- data.frame(
  Gene = as.character(c("Total", "Total", "Total")),
  Comparison = as.character(c("2N/2N+ vs. Ts/Ts+", "2N+/Ts+ vs. 2N/Ts", "Ts vs. Ts+")),
  C1 = c(71, 41, 92),
  C2 = c(56, 47, 92),
  C3 = c(72, 37, 87),
  C4 = c(18, 19, 44),
  C5 = c(41, 23, 55),
  C6 = c(67, 47, 102),
  C7 = c(32, 32, 57),
  C8 = c(93, 47, 78),
  C9 = c(66, 55, 106),
  C10 = c(67,51, 86),
  C11 = c(78,50,101),
  C12 = c(62, 55, 109),
  C13 = c(48, 42, 100),
  C14 = c(52, 64, 115),
  C15 = c(54, 60, 109),
  C16 = c(58, 50, 103),
  C17 = c(6, 2, 11),
  C18 = c(57, 19, 54),
  C19 = c(37, 30, 65),
  C20 = c(70, 47, 98),
  C21 = c(76, 42, 99),
  C22 = c(72, 51, 84),
  C23 = c(61, 42, 97),
  C24 = c(25, 20, 40),
  C25 = c(28, 38, 47),
  stringsAsFactors = FALSE # Prevent automatic factor conversion
)

# Check the structure to ensure types are correct
str(data)

# Remove the 'Comparison' column temporarily for reshaping
data_temp <- data %>% select(-Comparison)

# Reshape the data to long format
data_long <- data_temp %>%
  pivot_longer(cols = starts_with("C"), names_to = "Cluster", values_to = "Total")

# Re-add the 'Comparison' column
data_long <- data_long %>%
  mutate(Comparison = rep(data$Comparison, each = ncol(data_temp) - 1))

# Check the structure of the long format data
str(data_long)

# Define your custom color palette
colors <- c("#56B4E9", "#F0E442", "#009E73", "#DDA0DD")
#colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
#colors <- c("#FF6F61", "#4CAF50", "#008C9E", "#9B59B6")

# Create the bar graph with the custom color palette
ggplot(data_long, aes(x = Cluster, y = Total, fill = Comparison)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total Number of Significant DA Genes by Cluster",
       x = "Cluster",
       y = "Total DA Genes") +
  scale_fill_manual(values = colors) +  # Apply the custom color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability



########################################

## Add y-axis for total cell counts per cluster

# Load necessary libraries
library(tidyverse)

# Total cell counts per cluster
total_cells <- c(257, 180, 10313, 1005, 51, 2027, 61, 8020, 86, 975, 2908, 872,
                 409, 2136, 560, 1677, 1916, 44381, 976, 3426, 10686, 6750, 3513, 874, 396)

# Create the data frame for DA gene comparisons
data <- data.frame(
  Gene = as.character(c("Total", "Total", "Total")),
  Comparison = as.character(c("2N/2N+ vs. Ts/Ts+", "2N+/Ts+ vs. 2N/Ts", "Ts vs. Ts+")),
  C1 = c(71, 41, 92),
  C2 = c(56, 47, 92),
  C3 = c(72, 37, 87),
  C4 = c(18, 19, 44),
  C5 = c(41, 23, 55),
  C6 = c(67, 47, 102),
  C7 = c(32, 32, 57),
  C8 = c(93, 47, 78),
  C9 = c(66, 55, 106),
  C10 = c(67, 51, 86),
  C11 = c(78, 50, 101),
  C12 = c(62, 55, 109),
  C13 = c(48, 42, 100),
  C14 = c(52, 64, 115),
  C15 = c(54, 60, 109),
  C16 = c(58, 50, 103),
  C17 = c(6, 2, 11),
  C18 = c(57, 19, 54),
  C19 = c(37, 30, 65),
  C20 = c(70, 47, 98),
  C21 = c(76, 42, 99),
  C22 = c(72, 51, 84),
  C23 = c(61, 42, 97),
  C24 = c(25, 20, 40),
  C25 = c(28, 38, 47),
  stringsAsFactors = FALSE
)

# Reshape the comparison data to long format
data_long <- data %>%
  pivot_longer(cols = C1:C25, names_to = "Cluster", values_to = "Count")

# Create a data frame for the total cell counts, aligning clusters
total_cells_data <- data.frame(
  Gene = "Total",
  Comparison = "Total Cells",
  Cluster = paste0("C", 1:25),
  Count = total_cells
)

# Combine the comparison and total cell count data
data_combined <- bind_rows(data_long, total_cells_data)

# Ensure proper cluster ordering
data_combined$Cluster <- factor(data_combined$Cluster,
                                levels = paste0("C", 1:25))

# Convert Comparison to a factor to control legend order
data_combined$Comparison <- factor(data_combined$Comparison,
                                   levels = c("2N/2N+ vs. Ts/Ts+", "2N+/Ts+ vs. 2N/Ts",
                                              "Ts vs. Ts+", "Total Cells"))

# Define custom colors for the comparisons and total cells
colors <- c("2N/2N+ vs. Ts/Ts+" = "#56B4E9",
            "2N+/Ts+ vs. 2N/Ts" = "#F0E442",
            "Ts vs. Ts+" = "#009E73",
            "Total Cells" = "#D55E00")

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
    title = "DA Genes and Total Cell Counts per Cluster",
    x = "Cluster",
    fill = "Comparison"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
