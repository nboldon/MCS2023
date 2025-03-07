# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_Boxplot_2024-04-11.csv")

library(ggplot2)

########################################

ggplot(boxData, aes(x=treatment, y=cellAbund, shape=treatment, color=treatment)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex=1.2) +
  scale_shape_manual(values=c(15,16,17,18)) +
  facet_wrap(~cluster, scales = "free")

#Saved as cellAbund_Boxplots_2024-04-11

########################################

ggplot(boxData, aes(x=treatment, y=cellAbund, shape=treatment, color=treatment)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_jitter(position=position_dodge(0.8), cex=1.2) +
  scale_shape_manual(values=c(15,16,17,18)) +
  facet_wrap(~cluster, scales = "free")+
  geom_boxplot(alpha=0.05) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black") +
  theme(panel.grid.major.y = element_blank(),  
        panel.grid.major.x = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adjust panel borders
        panel.background = element_rect(fill = "white", color = NA),  # Adjust panel background
        strip.background = element_rect(fill = "lightgray", color = "black"),  # Adjust facet strip background
        strip.text = element_text(color = "black", size = 10),  # Adjust facet strip text size and color
        axis.text = element_text(size = 10),  # Adjust axis text size
        axis.title = element_text(size = 12)  # Adjust axis title text size
  )

#Saved as cellAbund_Boxplots_wMean_2024-04-11

##########

library(ggplot2)
library(ggpubr)

# Remove NA values
boxData <- na.omit(boxData)

# Compute pairwise t-test comparisons between treatment groups
pairwise_comparisons <- pairwise.t.test(boxData$cellAbund, boxData$treatment, p.adjust.method = "bonferroni")

# Extract p-values and adjust for NA values
p_values <- pairwise_comparisons$p.value[!is.na(pairwise_comparisons$p.value)]
treatments <- pairwise_comparisons$comparison[!is.na(pairwise_comparisons$p.value)]

# Check if there are significant pairwise comparisons
if(length(treatments) > 0) {
  # Create a dataframe with significant p-values and corresponding treatments
  signif_df <- data.frame(treatment = treatments, p_value = p_values)
  significant_pairs <- signif_df[signif_df$p_value < 0.05, ]
  
  # Generate the plot
  ggplot(boxData, aes(x = treatment, y = cellAbund)) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_jitter(position = position_dodge(0.8), cex = 1.2) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    facet_wrap(~ cluster, scales = "free_y") +
    geom_boxplot(alpha = 0.05) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black") +
    stat_pvalue_manual(
      data = significant_pairs,
      aes(label = paste0("p = ", round(p_value, 3))),
      step.increase = 0.4,
      hide.ns = TRUE,
      tip.length = 0
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
} else {
  # If no significant comparisons, print a message or take alternative action
  print("No significant pairwise comparisons found.")
}

## No significant pairwise comparisons found

#####################

library(ggpubr)
library(rstatix)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Working Docs/Cell-Abund_Boxplot_2024-04-11.csv")

## Data Preparation

# Transform `dose` into factor variable
df <- boxData
df$treatment <- as.factor(df$treatment)
# Add a random grouping variable
df$group <- factor(rep(c("grp1", "grp2"), 700))
# Add some extremely high values in column 1 at rows c(1, 3, 5).
df[c(1, 3, 5),  1] <- c(500, 495, 505)
head(df, 3)

## Multipanel plots containing two groups by panel
# Facet wrap by one variable

# Statistical tests
# Facet by the cluster variable and compare the levels of the treatment variable on the x-axis.

stat.test <- boxData %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

#########

## Multipanel plots containing three or more groups by panel
# Facet wrap by one variable
# Perform all pairwise comparisons

#Group by the cluster variable and then perform pairwise comparisons between the levels of treatment variable.

#Statistical test:
  
stat.test <- boxData %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment)
stat.test

#Add the p-values onto the plots. 
#The ggplot2 function scale_y_continuous(expand = expansion(mult = c(0, 0.1))) is used to add more spaces between labels and the plot top border

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
ggboxplot(boxData, x = "treatment", y = "cellAbund", fill = "#FC4E07", facet.by = "cluster") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#*!* C8 shows sig differences

######

## Pairwise comparisons against a reference group

#Statistical test:
  
stat.test <- boxData %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment, ref.group = "2N")
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
ggboxplot(boxData, x = "treatment", y = "cellAbund", fill = "#FC4E07", facet.by = "cluster") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#*!* C20 shows sig differences

#########################

## Box plots

# Create a box plot
bxp <- ggboxplot(
  boxData, x = "treatment", y = "cellAbund", fill = "treatment", 
  facet.by = "cluster"
)

# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(x = "treatment")
bxp + stat_pvalue_manual(stat.test)

###########

stat.test <- boxData %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

# Make the facet scale free and add jitter points
# Move down the bracket using `bracket.nudge.y`
# Hide ns (non-significant)
# Show adjusted p-values and significance levels
# Add 10% spaces between the p-value labels and the plot border
bxp <- ggboxplot(
  boxData, x = "treatment", y = "cellAbund", fill = "#00AFBB", 
  facet.by = "cluster", scales = "free_y", add = "jitter"
)
bxp +  
  stat_pvalue_manual(
    stat.test, y.position = 0.5, bracket.nudge.y = -2, hide.ns = TRUE,
    label = "{p.adj}{p.adj.signif}"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

##################################

library(ggpubr)
library(rstatix)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_Boxplot_2N-Ts_2024-04-11.csv")
# C1 removed to match the number of samples in the random grouping variable
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_TsvsTs+_Boxplot_2024-04-22.csv")
# C11 removed from data due to 4 samples missing values
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_2NvsTs+_Boxplot_2024-04-22.csv")
# C5 removed from data due to samples missing values
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_2Nvs2N+_Boxplot_2024-04-22.csv")
# C5 removed from data due to samples missing values
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_Boxplot_MCSvsCON_2024-04-22.csv")
boxData <- read.csv("/Volumes/DataBox/MCS2023/Cell Abundance/Cell-Abund_Boxplot_MCS-ALLvsTs-ALL_2024-04-22.csv")

## Data Preparation

# Transform `dose` into factor variable
df <- boxData
df$treatment <- as.factor(df$treatment)
# Add a random grouping variable
df$group <- factor(rep(c("grp1", "grp2"), 700))
# Add some extremely high values in column 1 at rows c(1, 3, 5).
df[c(1, 3, 5),  1] <- c(500, 495, 505)
head(df, 3)

## Multipanel plots containing two groups by panel
# Facet wrap by one variable

# Statistical tests
# Facet by the cluster variable and compare the levels of the treatment variable on the x-axis.

stat.test <- df %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

## Box plots

# Create a box plot
bxp <- ggboxplot(
  df, x = "treatment", y = "cellAbund", fill = "treatment", 
  facet.by = "cluster"
)

# Make facet and add p-values
stat.test <- stat.test %>% add_xy_position(x = "treatment")
bxp + stat_pvalue_manual(stat.test)

# Make the facet scale free and add jitter points
# Move down the bracket using `bracket.nudge.y`
# Hide ns (non-significant)
# Show adjusted p-values and significance levels
# Add 10% spaces between the p-value labels and the plot border
bxp <- ggboxplot(
  df, x = "treatment", y = "cellAbund", fill = "treatment", 
  facet.by = "cluster", scales = "free", add = "jitter"
)
bxp +  
  stat_pvalue_manual(
    stat.test, bracket.nudge.y = -2, hide.ns = TRUE,
    label = "{p.adj}{p.adj.signif}"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#Saved as cellAbund_Boxplots_2N-Ts_2024-04-11.pdf
#Saved as cellAbund_Boxplots_Ts-Ts+_2024-04-22.pdf
#Saved as cell-Abund_2NvsTs+_Boxplot_2024-04-22.pdf
#Saved as cell-Abund_2Nvs2N+_Boxplot_2024-04-22.pdf
#Saved as cell-Abund_MCSvsCON_Boxplot_2024-04-22.pdf
#Saved as cell-Abund_MCS-ALLvsTs-ALL_Boxplot_2024-04-22.pdf

########################################

## Compare more than 2 groups

# Global test
compare_means(cellAbund ~ treatment,  data = boxData, method = "anova")

# Default method = "kruskal.test" for multiple groups
ggboxplot(boxData, x = "treatment", y = "cellAbund",
          color = "treatment", palette = "jco")+
  stat_compare_means()
# Change method to anova
ggboxplot(boxData, x = "treatment", y = "cellAbund",
          color = "treatment", palette = "jco")+
  stat_compare_means(method = "anova")

##########

# Multiple pairwise tests against all (base-mean)
# Comparison of each group against base-mean
compare_means(cellAbund ~ treatment,  data = boxData, ref.group = ".all.",
              method = "t.test")
# Visualize
ggboxplot(boxData, x = "treatment", y = "cellAbund",
          color = "treatment", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 0.55)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                  # Pairwise comparison against all

########################################

library(ggpubr)
library(rstatix)

stat.test <- boxData %>%
  group_by(cluster) %>%
  tukey_hsd(cellAbund ~ treatment) 
stat.test 

#####

stat.test <- stat.test %>% 
  add_xy_position(x = "treatment", fun = "mean_se", scales = "free")
bp <- ggbarplot(
  boxData, x = "treatment", y = "cellAbund", fill = "#00AFBB",
  add = "mean_se", facet.by = "cluster", scales = "free"
) 
bp +
  stat_pvalue_manual(
    stat.test, hide.ns = TRUE, tip.length = 0,
    step.increase = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#Saved as: CellAbund_Tukey_hsd.pdf

########################################

stat.test <- boxData %>%
  group_by(cluster) %>%
  t_test(cellAbund ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(
  boxData, x = "treatment", y = "cellAbund", 
  color = "treatment", palette = "dark2"
)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "treatment", dodge = 0.8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0, text_color = "black"
)

# Add 10% spaces between the p-value labels and the plot border
bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

########################################

#stat_summary(fun.data="mean_sdl", mult=1, 
                 geom="crossbar", width=0.5)


########################################

## Grouping boxplots

# Plot in one graph
ggplot(boxData, aes(x = cluster, y = cellAbund,
                fill = treatment)) +
  geom_boxplot()

# Plot by cluster
ggplot(boxData, aes(x = cluster, y = cellAbund, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~cluster, scales = "free")

# To add stats
ggplot(boxData, aes(x = cluster, y = cellAbund,
                    fill = treatment)) +
  geom_boxplot(alpha=0.05) +
    stat_summary(fun.y = mean, geom = "point", shape = 20, size = 1, color = "yellow") +
  scale_fill_brewer(palette = "Set1") +
  stat_boxplot(geom = "errorbar", width = 0.25)

# To add points to the boxplot and plot by cluster

ggplot(boxData, aes(x = cluster, y = cellAbund,
                    fill = treatment)) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.4, position=position_jitter(0.2)) +
  geom_boxplot(alpha = 0.5) +
  theme(
    plot.title = element_text(size = 11)
  ) +
  facet_wrap(~cluster, scales = "free") +
  ggtitle("Cell Abundance by Cluster") +
  xlab("Cluster") + ylab("Cell Abundance") +
  stat_summary(fun.y = mean, geom = "point", shape = 20, size = 1, color = "yellow",
               position = position_dodge(width = 0.75))

# ggplot2.tidyverse.org

library(Hmisc)

# A set of useful summary functions is provided from the Hmisc package:
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, colour = "red", geom = geom, width = 0.2, ...)
}
d <- ggplot(boxData, aes(cluster, cellAbund)) + geom_point()
# The crossbar geom needs grouping to be specified when used with a continuous x axis.
d + stat_sum_df("mean_cl_boot", mapping = aes(group = cluster))

d + stat_sum_df("mean_sdl", mapping = aes(group = cluster))

d + stat_sum_df("mean_sdl", fun.args = list(mult = 1), mapping = aes(group = cluster))

###########################
###########################

# r-graph-gallery.com

## Violin plot with ggstatsplot
### Could not get to run

library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)

# Load the CSV file into a dataframe
boxData <- read.csv("/Volumes/DataBox/MCS2023/Working Docs/Cell-Abund_Boxplot_2024-04-11.csv")

plt <- ggbetweenstats(
    data = boxData,
    x = treatment,
    y = cellAbund,
)

plt <- plt +
  #Add labels and title
  labs(
    x = "Treatment",
    y = "Cell Abundance",
    title = "Cell Abundance by Treatment Group"
  ) +
  # Customizations
  theme(
    #This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
        family = "Lobster Two",
        size = 20,
        face = "bold",
        color = "#2a475e"
    ),
    # Statistical annoations below the main title
    plot.subtitle = element_text(
      family = "Roboto",
      size = 15,
      face = "bold",
      color = "#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
)

# 1. Remove axis ticks
# 2. Change default color of the axis lines with a lighter one
# 3. Remove most reference lines, only keep the major horizontal ones.
#    This reduces clutter, while keeping the reference for the variable being compared.
# 4. Set the panel and the background fill to the same light color.

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

ggsave(
  filename = ("./Boxplot_CellAbund-byCluster.png"),
  plot = plt,
  width = 8,
  height = 8,
  device = "png"
)

########################
