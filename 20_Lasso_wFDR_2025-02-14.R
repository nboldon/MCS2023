

## Complete for each cell type: astro, gaba, glut, oligo, endo, microglia

## Complete for each behavioral test: TaskA-45, TaskB-59, TaskC-48, TaskD-149


# Set working directory
setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_zscore-Stats")



## Enhanced version of the LASSO analysis with FDR correction and additional features

# Load additional required libraries
library(glmnet)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(caret)
library(patchwork) 


# Your data preparation code
df <- read.csv("astro_TaskA-45_StatsCombined_2025-02-03.csv")
behavioral_metrics <- c("Correct_Response", "Premature_Response", 
                        "Missed_Response_Window", "Wrong_Choice")
gene_columns <- setdiff(names(df), c("Sample", "Treatment", behavioral_metrics))
colnames(df)[colnames(df) %in% gene_columns] <- make.names(gene_columns)

# Basic LASSO function
perform_lasso <- function(data, response_var, predictor_vars) {
  # Prepare the data
  X <- as.matrix(data[, predictor_vars, drop = FALSE])
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
  
  # Calculate R-squared
  y_pred <- predict(final_model, newx = X)
  rsq <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
  
  # Return results
  return(list(
    response = response_var,
    lambda = best_lambda,
    coefficients = coef_matrix[selected_vars, , drop = FALSE],
    r_squared = rsq,
    model = final_model,
    cv_fit = cv_fit
  ))
}



# Function to calculate p-values for LASSO coefficients using bootstrapping
calculate_lasso_pvalues <- function(data, response_var, predictor_vars, n_bootstrap = 1000) {
  # Original LASSO fit
  X <- as.matrix(data[, predictor_vars, drop = FALSE])
  y <- data[[response_var]]
  
  # Remove NA values
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  # Standardize predictors
  X <- scale(X)
  
  # Original fit
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  original_coef <- coef(glmnet(X, y, alpha = 1, lambda = cv_fit$lambda.min))
  
  # Bootstrap procedure
  boot_coefs <- matrix(0, nrow = n_bootstrap, ncol = ncol(X))
  
  for(i in 1:n_bootstrap) {
    # Sample with replacement
    boot_indices <- sample(nrow(X), replace = TRUE)
    X_boot <- X[boot_indices, ]
    y_boot <- y[boot_indices]
    
    # Fit LASSO on bootstrap sample
    cv_boot <- cv.glmnet(X_boot, y_boot, alpha = 1, nfolds = 5)
    boot_fit <- glmnet(X_boot, y_boot, alpha = 1, lambda = cv_boot$lambda.min)
    boot_coefs[i,] <- as.vector(coef(boot_fit))[-1]  # Exclude intercept
  }
  
  # Calculate p-values
  p_values <- sapply(1:ncol(X), function(j) {
    coef_distribution <- boot_coefs[, j]
    if(as.vector(original_coef)[j + 1] == 0) {
      return(1)
    } else {
      # Two-tailed test
      2 * min(
        mean(coef_distribution >= 0),
        mean(coef_distribution <= 0)
      )
    }
  })
  
  return(data.frame(
    Variable = predictor_vars,
    Coefficient = as.vector(original_coef)[-1],
    P_value = p_values
  ))
}

# Enhanced LASSO function with p-values
perform_lasso_enhanced <- function(data, response_var, predictor_vars) {
  # Original LASSO analysis
  basic_results <- perform_lasso(data, response_var, predictor_vars)
  
  if(is.null(basic_results)) {
    return(NULL)
  }
  
  # Calculate p-values
  pvalue_results <- calculate_lasso_pvalues(data, response_var, predictor_vars)
  
  # Combine results
  results <- list(
    response = response_var,
    lambda = basic_results$lambda,
    coefficients = basic_results$coefficients,
    r_squared = basic_results$r_squared,
    model = basic_results$model,
    cv_fit = basic_results$cv_fit,
    p_values = pvalue_results
  )
  
  return(results)
}

# Function to apply FDR correction across all results
apply_fdr_correction <- function(lasso_results) {
  # Collect all p-values and their corresponding information
  all_results <- do.call(rbind, lapply(names(lasso_results), function(metric) {
    if(!is.null(lasso_results[[metric]])) {
      data.frame(
        Metric = metric,
        lasso_results[[metric]]$p_values,
        stringsAsFactors = FALSE
      )
    }
  }))
  
  # Apply FDR correction
  all_results$FDR <- p.adjust(all_results$P_value, method = "BH")
  all_results$Significant <- all_results$FDR < 0.05
  
  return(all_results)
}

# Enhanced plotting function with FDR results
plot_lasso_results_enhanced <- function(lasso_result, fdr_results) {
  if(is.null(lasso_result) || length(lasso_result$coefficients) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No variables selected by LASSO\n(all coefficients = 0)") +
             theme_minimal() +
             theme(axis.text = element_blank(),
                   axis.title = element_blank(),
                   axis.ticks = element_blank()))
  }
  
  # Combine coefficient and significance information
  coef_data <- data.frame(
    Variable = rownames(lasso_result$coefficients),
    Coefficient = as.numeric(lasso_result$coefficients)
  )
  
  # Add FDR information
  coef_data <- merge(
    coef_data,
    fdr_results[, c("Variable", "FDR", "Significant")],
    by = "Variable"
  )
  
  # Create enhanced plot
  ggplot(coef_data, 
         aes(x = reorder(Variable, abs(Coefficient)), 
             y = Coefficient,
             fill = Significant)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("LASSO Coefficients for", lasso_result$response),
      subtitle = sprintf("R² = %.3f", lasso_result$r_squared),
      x = "Variables",
      y = "Coefficient Value",
      fill = "FDR < 0.05"
    ) +
    theme(legend.position = "bottom")
}

# Main analysis pipeline
main_analysis <- function(df, behavioral_metrics, gene_columns) {
  # Perform enhanced LASSO for each behavioral metric
  lasso_results <- list()
  for(metric in behavioral_metrics) {
    lasso_results[[metric]] <- perform_lasso_enhanced(df, metric, gene_columns)
  }
  
  # Apply FDR correction
  fdr_results <- apply_fdr_correction(lasso_results)
  
  # Create enhanced plots
  lasso_plots <- list()
  for(metric in behavioral_metrics) {
    if(!is.null(lasso_results[[metric]])) {
      metric_fdr <- subset(fdr_results, Metric == metric)
      lasso_plots[[metric]] <- plot_lasso_results_enhanced(
        lasso_results[[metric]], 
        metric_fdr
      )
    }
  }
  
  # Combine plots
  if(length(lasso_plots) > 0) {
    final_plot <- wrap_plots(lasso_plots, ncol = 2) +
      plot_annotation(
        title = "LASSO Analysis Results with FDR Correction",
        subtitle = "Significant associations (FDR < 0.05) shown in blue"
      )
    
    # Save the plot
    ggsave("lasso_analysis_with_fdr.pdf", final_plot,
           width = 15, height = 15, limitsize = FALSE)
  }
  
  # Create comprehensive results table
  results_table <- fdr_results %>%
    arrange(Metric, FDR) %>%
    select(Metric, Variable, Coefficient, P_value, FDR, Significant)
  
  # Save results
  write.csv(results_table, "lasso_analysis_results_with_fdr.csv", row.names = FALSE)
  
  return(list(
    results = results_table,
    plots = lasso_plots,
    full_results = lasso_results
  ))
}


# Then simply call the main_analysis function from the enhanced code
results <- main_analysis(df, behavioral_metrics, gene_columns)

# The results object will contain:
# - results$results: The comprehensive results table with FDR correction
# - results$plots: The visualization plots
# - results$full_results: Complete LASSO analysis results


# View complete results
View(results$results)

# Look at significant associations only
significant_results <- results$results %>%
  filter(Significant == TRUE)
View(significant_results)



# View individual plots
results$plots$Correct_Response
results$plots$Premature_Response


# Get detailed results for a specific metric
results$full_results$Correct_Response




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



## Loop to run through all files with the same file name ending


## First make sure all the previous LASSO functions are loaded (perform_lasso, calculate_lasso_pvalues, etc.)


# Basic LASSO function
perform_lasso <- function(data, response_var, predictor_vars) {
  # Prepare the data
  X <- as.matrix(data[, predictor_vars, drop = FALSE])
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
  
  # Calculate R-squared
  y_pred <- predict(final_model, newx = X)
  rsq <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
  
  # Return results
  return(list(
    response = response_var,
    lambda = best_lambda,
    coefficients = coef_matrix[selected_vars, , drop = FALSE],
    r_squared = rsq,
    model = final_model,
    cv_fit = cv_fit
  ))
}



# Modified p-value calculation for LASSO coefficients
calculate_lasso_pvalues <- function(data, response_var, predictor_vars, n_bootstrap = 1000) {
  # Original LASSO fit
  X <- as.matrix(data[, predictor_vars, drop = FALSE])
  y <- data[[response_var]]
  
  # Remove NA values
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  # Standardize predictors
  X <- scale(X)
  
  # Original fit
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
  original_coef <- coef(glmnet(X, y, alpha = 1, lambda = cv_fit$lambda.min))
  
  # Bootstrap procedure
  boot_coefs <- matrix(0, nrow = n_bootstrap, ncol = ncol(X))
  
  for(i in 1:n_bootstrap) {
    # Sample with replacement
    boot_indices <- sample(nrow(X), replace = TRUE)
    X_boot <- X[boot_indices, ]
    y_boot <- y[boot_indices]
    
    # Fit LASSO on bootstrap sample
    cv_boot <- cv.glmnet(X_boot, y_boot, alpha = 1, nfolds = 5)
    boot_fit <- glmnet(X_boot, y_boot, alpha = 1, lambda = cv_boot$lambda.min)
    boot_coefs[i,] <- as.vector(coef(boot_fit))[-1]  # Exclude intercept
  }
  
  # Calculate p-values - MODIFIED to be less stringent for zero coefficients
  p_values <- sapply(1:ncol(X), function(j) {
    original_val <- as.vector(original_coef)[j + 1]  # +1 because of intercept
    coef_distribution <- boot_coefs[, j]
    
    if(original_val == 0) {
      # Check if coefficient is consistently zero in bootstrap samples
      prop_zero <- mean(coef_distribution == 0)
      if(prop_zero > 0.95) {
        return(1)  # Truly zero coefficient
      } else {
        # Calculate p-value based on bootstrap distribution
        if(mean(coef_distribution) > 0) {
          return(2 * mean(coef_distribution <= 0))
        } else {
          return(2 * mean(coef_distribution >= 0))
        }
      }
    } else {
      # For non-zero coefficients, use two-tailed test
      if(original_val > 0) {
        return(2 * mean(coef_distribution <= 0))
      } else {
        return(2 * mean(coef_distribution >= 0))
      }
    }
  })
  
  return(data.frame(
    Variable = predictor_vars,
    Coefficient = as.vector(original_coef)[-1],
    P_value = p_values
  ))
}



# Enhanced LASSO function with p-values
perform_lasso_enhanced <- function(data, response_var, predictor_vars) {
  # Original LASSO analysis
  basic_results <- perform_lasso(data, response_var, predictor_vars)
  
  if(is.null(basic_results)) {
    return(NULL)
  }
  
  # Calculate p-values
  pvalue_results <- calculate_lasso_pvalues(data, response_var, predictor_vars)
  
  # Combine results
  results <- list(
    response = response_var,
    lambda = basic_results$lambda,
    coefficients = basic_results$coefficients,
    r_squared = basic_results$r_squared,
    model = basic_results$model,
    cv_fit = basic_results$cv_fit,
    p_values = pvalue_results
  )
  
  return(results)
}



# Modified FDR correction function with adjustable threshold
apply_fdr_correction <- function(lasso_results, fdr_threshold = 0.1) {
  # Collect all p-values and their corresponding information
  all_results <- do.call(rbind, lapply(names(lasso_results), function(metric) {
    if(!is.null(lasso_results[[metric]])) {
      data.frame(
        Metric = metric,
        lasso_results[[metric]]$p_values,
        stringsAsFactors = FALSE
      )
    }
  }))
  
  # Apply FDR correction
  all_results$FDR <- p.adjust(all_results$P_value, method = "BH")
  all_results$Significant <- all_results$FDR < fdr_threshold
  
  return(all_results)
}



# Enhanced plotting function with FDR results
plot_lasso_results_enhanced <- function(lasso_result, fdr_results) {
  if(is.null(lasso_result) || length(lasso_result$coefficients) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No variables selected by LASSO\n(all coefficients = 0)") +
             theme_minimal() +
             theme(axis.text = element_blank(),
                   axis.title = element_blank(),
                   axis.ticks = element_blank()))
  }
  
  # Combine coefficient and significance information
  coef_data <- data.frame(
    Variable = rownames(lasso_result$coefficients),
    Coefficient = as.numeric(lasso_result$coefficients)
  )
  
  # Add FDR information
  coef_data <- merge(
    coef_data,
    fdr_results[, c("Variable", "FDR", "Significant")],
    by = "Variable"
  )
  
  # Create enhanced plot
  ggplot(coef_data, 
         aes(x = reorder(Variable, abs(Coefficient)), 
             y = Coefficient,
             fill = Significant)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue")) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste("LASSO Coefficients for", lasso_result$response),
      subtitle = sprintf("R² = %.3f", lasso_result$r_squared),
      x = "Variables",
      y = "Coefficient Value",
      fill = "FDR < 0.05"
    ) +
    theme(legend.position = "bottom")
}




# Modified main analysis to use less stringent parameters
main_analysis <- function(df, behavioral_metrics, gene_columns, fdr_threshold = 0.1) {
  # Perform enhanced LASSO for each behavioral metric
  lasso_results <- list()
  for(metric in behavioral_metrics) {
    lasso_results[[metric]] <- perform_lasso_enhanced(df, metric, gene_columns)
  }
  
  # Apply FDR correction with custom threshold
  fdr_results <- apply_fdr_correction(lasso_results, fdr_threshold)
  
  # Create enhanced plots
  lasso_plots <- list()
  for(metric in behavioral_metrics) {
    if(!is.null(lasso_results[[metric]])) {
      metric_fdr <- subset(fdr_results, Metric == metric)
      lasso_plots[[metric]] <- plot_lasso_results_enhanced(
        lasso_results[[metric]], 
        metric_fdr
      )
    }
  }
  
  # Combine plots
  if(length(lasso_plots) > 0) {
    final_plot <- wrap_plots(lasso_plots, ncol = 2) +
      plot_annotation(
        title = "LASSO Analysis Results with FDR Correction",
        subtitle = paste0("Significant associations (FDR < ", fdr_threshold, ") shown in blue")
      )
    
    # Save the plot
    ggsave("lasso_analysis_with_fdr.pdf", final_plot,
           width = 15, height = 15, limitsize = FALSE)
  }
  
  # Create comprehensive results table
  results_table <- fdr_results %>%
    arrange(Metric, FDR) %>%
    select(Metric, Variable, Coefficient, P_value, FDR, Significant)
  
  # Save results
  write.csv(results_table, "lasso_analysis_results_with_fdr.csv", row.names = FALSE)
  
  return(list(
    results = results_table,
    plots = lasso_plots,
    full_results = lasso_results
  ))
}




# Set working directory
#setwd("/Volumes/DataBox/MCS2023/Stats/Treatment_CellType_zscore-Stats")
setwd("/Volumes/DataBox/MCS2023/Stats/Pearson_Treatment_Cluster_zscore-Stats")


# Load all required packages
library(glmnet)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(caret)
library(patchwork)
library(stringr)  


# Modified main analysis to use less stringent parameters
main_analysis <- function(df, behavioral_metrics, gene_columns, fdr_threshold = 0.1) {
  # Perform enhanced LASSO for each behavioral metric
  lasso_results <- list()
  for(metric in behavioral_metrics) {
    lasso_results[[metric]] <- perform_lasso_enhanced(df, metric, gene_columns)
  }
  
  # Apply FDR correction with custom threshold
  fdr_results <- apply_fdr_correction(lasso_results, fdr_threshold)
  
  # Create enhanced plots
  lasso_plots <- list()
  for(metric in behavioral_metrics) {
    if(!is.null(lasso_results[[metric]])) {
      metric_fdr <- subset(fdr_results, Metric == metric)
      lasso_plots[[metric]] <- plot_lasso_results_enhanced(
        lasso_results[[metric]], 
        metric_fdr
      )
    }
  }
  
  # Combine plots
  if(length(lasso_plots) > 0) {
    final_plot <- wrap_plots(lasso_plots, ncol = 2) +
      plot_annotation(
        title = "LASSO Analysis Results with FDR Correction",
        subtitle = paste0("Significant associations (FDR < ", fdr_threshold, ") shown in blue")
      )
    
    # Save the plot
    ggsave("lasso_analysis_with_fdr.pdf", final_plot,
           width = 15, height = 15, limitsize = FALSE)
  }
  
  # Create comprehensive results table
  results_table <- fdr_results %>%
    arrange(Metric, FDR) %>%
    select(Metric, Variable, Coefficient, P_value, FDR, Significant)
  
  # Save results
  write.csv(results_table, "lasso_analysis_results_with_fdr.csv", row.names = FALSE)
  
  return(list(
    results = results_table,
    plots = lasso_plots,
    full_results = lasso_results
  ))
}

# Modified run_lasso_analyses function to pass the FDR threshold
run_lasso_analyses <- function(pattern = "*_StatsCombined_2025-02-03.csv", fdr_threshold = 0.1) {
  # [Keep all existing code...]
  
  # Get list of files with more verbose output
  files <- list.files(pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found matching pattern:", pattern, "\n")
    cat("Current working directory:", getwd(), "\n")
    stop("No matching files found. Check file pattern and working directory.")
  }
  
  cat("Found", length(files), "files to process:\n")
  for (file in files) {
    cat("  -", basename(file), "\n")
  }
  
  # Create main results directory
  main_dir <- "LASSO_Analysis_Results"
  dir.create(main_dir, showWarnings = FALSE)
  
  # Process each file
  all_results <- list()
  for (file in files) {
    tryCatch({
      cat("\nAttempting to process:", basename(file), "\n")
      results <- process_single_file(file, fdr_threshold)  # Pass the threshold
      all_results[[basename(file)]] <- results
    }, error = function(e) {
      cat("Error processing", basename(file), ":", conditionMessage(e), "\n")
      cat("Full error:\n")
      print(e)
    })
  }
  
  # [Keep all existing summary code...]
  if (length(all_results) > 0) {
    overall_summary <- do.call(rbind, lapply(all_results, function(x) {
      if (!is.null(x$results)) {
        data.frame(
          Cell_Type = x$cell_type,
          Task = x$task,
          Total_Significant = sum(x$results$results$Significant)
        )
      }
    }))
    
    # Only save if we have data
    if (!is.null(overall_summary) && nrow(overall_summary) > 0) {
      write.csv(overall_summary,
                file.path(main_dir, "LASSO_overall_summary.csv"),
                row.names = FALSE)
      
      # Create summary plot
      summary_plot <- ggplot(overall_summary,
                             aes(x = Cell_Type, y = Total_Significant, fill = Task)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        labs(title = "Summary of Significant Associations Across All Analyses",
             subtitle = paste0("FDR threshold: ", fdr_threshold),
             y = "Number of Significant Associations",
             x = "Cell Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(main_dir, "LASSO_overall_summary_plot.pdf"),
             summary_plot, width = 12, height = 8)
    } else {
      cat("No valid results to summarize\n")
    }
  } else {
    cat("No results were successfully processed\n")
  }
  
  return(all_results)
}

# Modified process_single_file to accept the FDR threshold
process_single_file <- function(file_path, fdr_threshold = 0.1) {
  tryCatch({
    # [Keep existing file reading code...]
    
    # Extract cell type and task from filename
    file_info <- str_match(basename(file_path), "(.+?)_(.+?)_StatsCombined")[c(2,3)]
    
    if (is.na(file_info[1]) || is.na(file_info[2])) {
      stop("Couldn't extract cell type and task from filename: ", basename(file_path))
    }
    
    cell_type <- file_info[1]
    task <- file_info[2]
    
    # Create output directory for this analysis
    output_dir <- paste0(cell_type, "_", task, "_LASSO_results")
    dir.create(output_dir, showWarnings = FALSE)
    
    # Read the data with error checking
    cat("Reading file:", basename(file_path), "\n")
    df <- tryCatch({
      read.csv(file_path)
    }, error = function(e) {
      stop("Failed to read CSV file: ", conditionMessage(e))
    })
    
    cat("File loaded successfully. Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
    
    # Define behavioral metrics
    behavioral_metrics <- c("Correct_Response", "Premature_Response",
                            "Missed_Response_Window", "Wrong_Choice")
    
    # Check if behavioral metrics exist in the data
    missing_metrics <- behavioral_metrics[!behavioral_metrics %in% names(df)]
    if (length(missing_metrics) > 0) {
      cat("Warning: Some behavioral metrics are missing from the data:", 
          paste(missing_metrics, collapse=", "), "\n")
      behavioral_metrics <- behavioral_metrics[behavioral_metrics %in% names(df)]
      if (length(behavioral_metrics) == 0) {
        stop("None of the expected behavioral metrics found in the data")
      }
    }
    
    # Get gene columns
    gene_columns <- setdiff(names(df), c("Sample", "Treatment", behavioral_metrics))
    if (length(gene_columns) == 0) {
      stop("No gene columns identified in the data")
    }
    
    # Make valid column names
    colnames(df)[colnames(df) %in% gene_columns] <- make.names(gene_columns)
    
    # Update gene_columns with new names
    gene_columns <- setdiff(names(df), c("Sample", "Treatment", behavioral_metrics))
    
    cat("Found", length(behavioral_metrics), "behavioral metrics and", 
        length(gene_columns), "gene columns\n")
    
    # Run main analysis with custom FDR threshold
    results <- main_analysis(df, behavioral_metrics, gene_columns, fdr_threshold)
    
    # [Keep existing results saving code...]
    
    # Save results with specific filenames
    if (!is.null(results)) {
      # Save results table with specific filename
      results_filename <- file.path(output_dir,
                                    paste0(cell_type, "_", task, "_LASSO_results_with_fdr.csv"))
      write.csv(results$results, results_filename, row.names = FALSE)
      
      # Save plot with specific filename
      plot_filename <- file.path(output_dir,
                                 paste0(cell_type, "_", task, "_LASSO_plot_with_fdr.pdf"))
      if (length(results$plots) > 0) {
        final_plot <- wrap_plots(results$plots, ncol = 2) +
          plot_annotation(
            title = paste(cell_type, task, "LASSO Analysis Results with FDR Correction"),
            subtitle = paste0("Significant associations (FDR < ", fdr_threshold, ") shown in blue")
          )
        ggsave(plot_filename, final_plot, width = 15, height = 15, limitsize = FALSE)
      }
      
      # Create and save summary statistics
      summary_stats <- data.frame(
        Cell_Type = cell_type,
        Task = task,
        Total_Genes = length(gene_columns),
        Total_Significant = sum(results$results$Significant)
      )
      
      # Add per-metric statistics if available
      if (!is.null(results$results$Metric)) {
        sig_by_metric <- tapply(results$results$Significant,
                                results$results$Metric,
                                sum)
        for (metric in names(sig_by_metric)) {
          summary_stats[[paste0("Significant_", metric)]] <- sig_by_metric[[metric]]
        }
      }
      
      summary_filename <- file.path(output_dir,
                                    paste0(cell_type, "_", task, "_LASSO_summary.csv"))
      write.csv(summary_stats, summary_filename, row.names = FALSE)
      
      # Print progress
      cat("Analysis completed for:", basename(file_path), "\n")
      cat("Total significant associations:", sum(results$results$Significant), "\n")
      cat("Results saved in:", output_dir, "\n\n")
    } else {
      cat("No results returned from main_analysis for:", basename(file_path), "\n")
    }
    
    return(list(
      cell_type = cell_type,
      task = task,
      results = results
    ))
  }, error = function(e) {
    cat("Error in process_single_file for", basename(file_path), ":\n")
    cat(conditionMessage(e), "\n")
    return(NULL)
  })
}


# Run with a more general pattern
all_results <- run_lasso_analyses(pattern = "*_StatsCombined_2025-02-03.csv")


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



sessionInfo()

R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Denver
tzcode source: internal

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] stringr_1.5.1     boot_1.3-31       patchwork_1.3.0   caret_7.0-1       lattice_0.22-6   
[6] glmnet_4.1-8      Matrix_1.7-2      ggpubr_0.6.0      viridis_0.6.5     viridisLite_0.4.2
[11] ggplot2_3.5.1     dplyr_1.1.4       tidyr_1.3.1      

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_1.8.9          shape_1.4.6.1          
[5] magrittr_2.0.3          farver_2.1.2            rmarkdown_2.29          ragg_1.3.3             
[9] GlobalOptions_0.1.2     fs_1.6.5                zlibbioc_1.50.0         vctrs_0.6.5            
[13] memoise_2.0.1           ggtree_3.12.0           rstatix_0.7.2           htmltools_0.5.8.1      
[17] broom_1.0.7             pROC_1.18.5             Formula_1.2-5           gridGraphics_0.5-1     
[21] parallelly_1.42.0       plyr_1.8.9              httr2_1.1.0             lubridate_1.9.4        
[25] cachem_1.1.0            igraph_2.1.4            lifecycle_1.0.4         iterators_1.0.14       
[29] pkgconfig_2.0.3         R6_2.6.0                fastmap_1.2.0           gson_0.1.0             
[33] GenomeInfoDbData_1.2.12 future_1.34.0           clue_0.3-66             digest_0.6.37          
[37] aplot_0.2.4             enrichplot_1.24.4       colorspace_2.1-1        AnnotationDbi_1.66.0   
[41] S4Vectors_0.42.1        textshaping_1.0.0       RSQLite_2.3.9           labeling_0.4.3         
[45] timechange_0.3.0        abind_1.4-8             httr_1.4.7              polyclip_1.10-7        
[49] compiler_4.4.0          bit64_4.6.0-1           withr_3.0.2             doParallel_1.0.17      
[53] backports_1.5.0         BiocParallel_1.38.0     carData_3.0-5           DBI_1.2.3              
[57] ggforce_0.4.2           R.utils_2.12.3          ggsignif_0.6.4          lava_1.8.1             
[61] MASS_7.3-64             rappdirs_0.3.3          rjson_0.2.23            ModelMetrics_1.2.2.2   
[65] tools_4.4.0             ape_5.8-1               scatterpie_0.2.4        future.apply_1.11.3    
[69] nnet_7.3-20             R.oo_1.27.0             glue_1.8.0              nlme_3.1-167           
[73] GOSemSim_2.30.2         grid_4.4.0              shadowtext_0.1.4        cluster_2.1.8          
[77] reshape2_1.4.4          recipes_1.1.1           fgsea_1.30.0            generics_0.1.3         
[81] gtable_0.3.6            class_7.3-23            R.methodsS3_1.8.2       data.table_1.16.4      
[85] car_3.1-3               tidygraph_1.3.1         XVector_0.44.0          BiocGenerics_0.50.0    
[89] ggrepel_0.9.6           foreach_1.5.2           pillar_1.10.1           yulab.utils_0.2.0      
[93] circlize_0.4.16         splines_4.4.0           tweenr_2.0.3            treeio_1.28.0          
[97] survival_3.8-3          bit_4.5.0.1             tidyselect_1.2.1        GO.db_3.19.1           
[101] ComplexHeatmap_2.20.0   Biostrings_2.72.1       knitr_1.49              gridExtra_2.3          
[105] IRanges_2.38.1          stats4_4.4.0            xfun_0.50               graphlayouts_1.2.2     
[109] Biobase_2.64.0          hardhat_1.4.1           timeDate_4041.110       matrixStats_1.5.0      
[113] pheatmap_1.0.12         stringi_1.8.4           UCSC.utils_1.0.0        lazyeval_0.2.2         
[117] ggfun_0.1.8             yaml_2.3.10             pacman_0.5.1            evaluate_1.0.3         
[121] codetools_0.2-20        ggraph_2.2.1            tibble_3.2.1            qvalue_2.36.0          
[125] ggplotify_0.1.2         cli_3.6.3               rpart_4.1.24            systemfonts_1.2.1      
[129] munsell_0.5.1           Rcpp_1.0.14             GenomeInfoDb_1.40.1     globals_0.16.3         
[133] png_0.1-8               parallel_4.4.0          gower_1.0.2             blob_1.2.4             
[137] clusterProfiler_4.12.6  DOSE_3.30.5             listenv_0.9.1           tidytree_0.4.6         
[141] ipred_0.9-15            prodlim_2024.06.25      scales_1.3.0            purrr_1.0.4            
[145] crayon_1.5.3            GetoptLong_1.0.5        rlang_1.1.5             cowplot_1.1.3          
[149] fastmatch_1.1-6         KEGGREST_1.44.1        

