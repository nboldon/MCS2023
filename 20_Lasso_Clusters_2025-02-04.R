


## Complete for select clusters: 3, 8, 11, 14, 18, 22

## Complete for each behavioral test: TaskA-45, TaskB-59, TaskC-48, TaskD-149


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_Cluster_zscore-Stats")

# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(glmnet)  # For LASSO regression
library(caret)   # For data preprocessing
library(patchwork)  # For combining plots

# Read the data
df <- read.csv("C18_TaskD-149_StatsCombined_2025-02-03.csv")

# Define behavioral metrics
behavioral_metrics <- c("Correct_Response", "Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")

# Get gene columns
gene_columns <- setdiff(names(df), c("Sample", "Treatment", behavioral_metrics))

# Ensure valid column names
colnames(df)[colnames(df) %in% gene_columns] <- make.names(gene_columns)

# Function to perform LASSO regression for each behavioral metric
perform_lasso <- function(data, response_var, predictor_vars) {
  # Prepare the data
  X <- as.matrix(data[, predictor_vars, drop = FALSE])  # Ensure matrix format
  y <- data[[response_var]]
  
  # Remove rows with NA
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  # Check if there is enough data
  if (nrow(X) < 5 || length(unique(y)) < 2) {
    return(NULL)  # Not enough data to fit LASSO
  }
  
  # Standardize predictors
  X <- scale(X)
  
  # Fit LASSO with cross-validation
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  
  # Fit model with optimal lambda
  best_lambda <- cv_fit$lambda.min
  final_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  
  # Get coefficients (excluding intercept)
  coef_matrix <- as.matrix(coef(final_model))
  selected_vars <- rownames(coef_matrix)[which(coef_matrix != 0 & rownames(coef_matrix) != "(Intercept)")]
  
  # Calculate R-squared for the model
  y_pred <- predict(final_model, newx = X)
  rsq <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
  
  # Return results
  return(list(
    response = response_var,
    lambda = best_lambda,
    coefficients = coef_matrix[selected_vars, , drop = FALSE],  # Ensure correct format
    r_squared = rsq,
    model = final_model,
    cv_fit = cv_fit
  ))
}

# Function to plot LASSO results
plot_lasso_results <- function(lasso_result) {
  # If no coefficients are selected, return an empty plot with a message
  if (is.null(lasso_result) || length(lasso_result$coefficients) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No variables selected by LASSO\n(all coefficients = 0)") +
             theme_minimal() +
             theme(axis.text = element_blank(),
                   axis.title = element_blank(),
                   axis.ticks = element_blank()) +
             labs(title = paste("LASSO Coefficients for", lasso_result$response),
                  subtitle = sprintf("R² = %.3f", ifelse(is.null(lasso_result), NA, lasso_result$r_squared))))
  }
  
  # Otherwise, create a coefficient plot
  coef_data <- data.frame(
    Variable = rownames(lasso_result$coefficients),
    Coefficient = as.numeric(lasso_result$coefficients)
  )
  
  ggplot(coef_data, aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("LASSO Coefficients for", lasso_result$response),
      x = "Variables",
      y = "Coefficient Value",
      subtitle = sprintf("R² = %.3f", lasso_result$r_squared)
    )
}

# Perform LASSO for each behavioral metric
lasso_results <- list()
lasso_plots <- list()

for (metric in behavioral_metrics) {
  lasso_results[[metric]] <- perform_lasso(df, metric, gene_columns)
  
  # Only plot if the result is not NULL
  if (!is.null(lasso_results[[metric]])) {
    lasso_plots[[metric]] <- plot_lasso_results(lasso_results[[metric]])
  }
}

# Combine plots
if (length(lasso_plots) > 0) {
  final_lasso_plot <- wrap_plots(lasso_plots, ncol = 2) +
    plot_annotation(title = "LASSO Analysis Results for Cluster 18 Task D-149 Behavioral Metrics")
  
  # Save the plot
  ggsave("C18-TaskD-149_lasso_analysis_2025-02-03.pdf", final_lasso_plot, 
         width = 15, height = 15, limitsize = FALSE)
}

# Create summary table of results
summary_table <- do.call(rbind, lapply(names(lasso_results), function(metric) {
  if (!is.null(lasso_results[[metric]])) {
    data.frame(
      Metric = metric,
      R_squared = lasso_results[[metric]]$r_squared,
      Selected_Genes = ifelse(length(lasso_results[[metric]]$coefficients) > 0, 
                              paste(rownames(lasso_results[[metric]]$coefficients), collapse = ", "), 
                              "None"),
      Lambda = lasso_results[[metric]]$lambda
    )
  } else {
    data.frame(
      Metric = metric,
      R_squared = NA,
      Selected_Genes = "Not enough data",
      Lambda = NA
    )
  }
}))

# Save summary table
write.csv(summary_table, "C18_TaskD-149_LASSO_Summary_2025-02-04.csv", row.names = FALSE)

print("LASSO analysis completed successfully.")





########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################



## Complete for select clusters: 3, 8, 11, 14, 18, 22

## Complete for each "Wrong" behavioral test: TaskA-45, TaskB-59, TaskC-48, TaskD-149


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_Cluster_zscore-Stats")

# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(glmnet)  # For LASSO regression
library(caret)   # For data preprocessing
library(patchwork)  # For combining plots

# Read the data
df <- read.csv("C18_TaskC-48_WrongStatsCombined_2025-02-03.csv")

# Define behavioral metrics
behavioral_metrics <- c("Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")

# Get gene columns
gene_columns <- setdiff(names(df), c("Sample", "Treatment", behavioral_metrics))

# Ensure valid column names
colnames(df)[colnames(df) %in% gene_columns] <- make.names(gene_columns)

# Function to perform LASSO regression for each behavioral metric
perform_lasso <- function(data, response_var, predictor_vars) {
  # Prepare the data
  X <- as.matrix(data[, predictor_vars, drop = FALSE])  # Ensure matrix format
  y <- data[[response_var]]
  
  # Remove rows with NA
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  # Check if there is enough data
  if (nrow(X) < 5 || length(unique(y)) < 2) {
    return(NULL)  # Not enough data to fit LASSO
  }
  
  # Standardize predictors
  X <- scale(X)
  
  # Fit LASSO with cross-validation
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  
  # Fit model with optimal lambda
  best_lambda <- cv_fit$lambda.min
  final_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  
  # Get coefficients (excluding intercept)
  coef_matrix <- as.matrix(coef(final_model))
  selected_vars <- rownames(coef_matrix)[which(coef_matrix != 0 & rownames(coef_matrix) != "(Intercept)")]
  
  # Calculate R-squared for the model
  y_pred <- predict(final_model, newx = X)
  rsq <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
  
  # Return results
  return(list(
    response = response_var,
    lambda = best_lambda,
    coefficients = coef_matrix[selected_vars, , drop = FALSE],  # Ensure correct format
    r_squared = rsq,
    model = final_model,
    cv_fit = cv_fit
  ))
}

# Function to plot LASSO results
plot_lasso_results <- function(lasso_result) {
  # If no coefficients are selected, return an empty plot with a message
  if (is.null(lasso_result) || length(lasso_result$coefficients) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, 
                      label = "No variables selected by LASSO\n(all coefficients = 0)") +
             theme_minimal() +
             theme(axis.text = element_blank(),
                   axis.title = element_blank(),
                   axis.ticks = element_blank()) +
             labs(title = paste("LASSO Coefficients for", lasso_result$response),
                  subtitle = sprintf("R² = %.3f", ifelse(is.null(lasso_result), NA, lasso_result$r_squared))))
  }
  
  # Otherwise, create a coefficient plot
  coef_data <- data.frame(
    Variable = rownames(lasso_result$coefficients),
    Coefficient = as.numeric(lasso_result$coefficients)
  )
  
  ggplot(coef_data, aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("LASSO Coefficients for", lasso_result$response),
      x = "Variables",
      y = "Coefficient Value",
      subtitle = sprintf("R² = %.3f", lasso_result$r_squared)
    )
}

# Perform LASSO for each behavioral metric
lasso_results <- list()
lasso_plots <- list()

for (metric in behavioral_metrics) {
  lasso_results[[metric]] <- perform_lasso(df, metric, gene_columns)
  
  # Only plot if the result is not NULL
  if (!is.null(lasso_results[[metric]])) {
    lasso_plots[[metric]] <- plot_lasso_results(lasso_results[[metric]])
  }
}

# Combine plots
if (length(lasso_plots) > 0) {
  final_lasso_plot <- wrap_plots(lasso_plots, ncol = 2) +
    plot_annotation(title = "LASSO Analysis Results for Cluster 18 Task C-48 Wrong Behavioral Metrics")
  
  # Save the plot
  ggsave("C18-TaskC-48-Wrong_lasso_analysis_2025-02-03.pdf", final_lasso_plot, 
         width = 15, height = 15, limitsize = FALSE)
}

# Create summary table of results
summary_table <- do.call(rbind, lapply(names(lasso_results), function(metric) {
  if (!is.null(lasso_results[[metric]])) {
    data.frame(
      Metric = metric,
      R_squared = lasso_results[[metric]]$r_squared,
      Selected_Genes = ifelse(length(lasso_results[[metric]]$coefficients) > 0, 
                              paste(rownames(lasso_results[[metric]]$coefficients), collapse = ", "), 
                              "None"),
      Lambda = lasso_results[[metric]]$lambda
    )
  } else {
    data.frame(
      Metric = metric,
      R_squared = NA,
      Selected_Genes = "Not enough data",
      Lambda = NA
    )
  }
}))

# Save summary table
write.csv(summary_table, "C18_TaskC-48-Wrong_LASSO_Summary_2025-02-04.csv", row.names = FALSE)

print("LASSO analysis completed successfully.")

