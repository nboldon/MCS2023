# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)

setwd("/Volumes/DataBox/MCS2023/TxComp")

#############################################

## Combine GO results with treatment group accessibility levels - for one cluster

# Read in both CSV files
file1 <- read.csv("./C1_significant_t1-t2_t3-t4_zscores_2024-10-18.csv", stringsAsFactors = FALSE)
file2 <- read.csv("./C1_GO_t1-t2_t3-t4_2024-10-18.csv", stringsAsFactors = FALSE)

# Split the 'geneID' column by "/" in file2 and unnest the list to create individual gene rows
file2 <- file2 %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest(geneID)

# Join the two data frames based on the gene names
merged_data <- file1 %>%
  inner_join(file2, by = c("Gene" = "geneID")) %>%
  select(Gene, t1, t2, t3, t4, pvalue, ID, Description)

# Save the result as a new CSV file
write.csv(merged_data, "C1_GO_t1-t2_t3-t4_matched_genes_2024-10-29.csv", row.names = FALSE)

#############################################
#############################################
#############################################


## Combine GO BP results with treatment group accessibility levels 
## Loop for all clusters - 2N/2N+ vs Ts/Ts+


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t1-t2_t3-t4_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_t1-t2_t3-t4_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_t1-t2_t3-t4_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}

#############################################

## Combine GO BP results with treatment group accessibility levels 
## Loop for all clusters - 2N+/Ts+ vs 2N/Ts


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t2-t4_t1-t3_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_t2-t4_t1-t3_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_t2-t4_t1-t3_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}

#############################################

## Combine GO BP results with treatment group accessibility levels 
## Loop for all clusters - Ts vs Ts+


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t3_t4_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_t3_t4_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_t3_t4_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}


#############################################
#############################################
#############################################


## Combine GO MF results with treatment group accessibility levels 
## Loop for all clusters - 2N/2N+ vs Ts/Ts+


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t1-t2_t3-t4_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_MF_t1-t2_t3-t4_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_MF_t1-t2_t3-t4_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}

#############################################

## Combine GO MF results with treatment group accessibility levels 
## Loop for all clusters - 2N+/Ts+ vs 2N/Ts


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t2-t4_t1-t3_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_MF_t2-t4_t1-t3_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_MF_t2-t4_t1-t3_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}

#############################################

## Combine GO MF results with treatment group accessibility levels 
## Loop for all clusters - Ts vs Ts+


# Set the base path for the files
base_path <- "/Volumes/DataBox/MCS2023/Tx_Comp/"

# Loop through each cluster (C1 to C25)
for (i in 1:25) {
  # Create file names for each cluster
  file1_name <- paste0(base_path, "C", i, "_significant_t3_t4_zscores_2024-10-18.csv")
  file2_name <- paste0(base_path, "C", i, "_GO_MF_t3_t4_2024-10-18.csv")
  
  # Print file paths to confirm they are correct
  message("Attempting to read files for cluster C", i)
  message("File 1: ", file1_name)
  message("File 2: ", file2_name)
  
  # Check if files exist before reading them
  if (!file.exists(file1_name) || !file.exists(file2_name)) {
    message(paste("Skipping cluster C", i, " due to missing file(s)", sep = ""))
    next
  }
  
  # Read the files since they exist
  file1 <- read.csv(file1_name, stringsAsFactors = FALSE)
  file2 <- read.csv(file2_name, stringsAsFactors = FALSE)
  
  # Process the file if both files are successfully read
  file2 <- file2 %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(geneID)
  
  # Join the two data frames based on the gene names
  merged_data <- file1 %>%
    inner_join(file2, by = c("Gene" = "geneID")) %>%
    select(Gene, t1, t2, t3, t4, pvalue, ID, Description)
  
  # Save the result as a new CSV file, dynamically naming it by cluster
  output_file <- paste0(base_path, "C", i, "_GO_MF_t3_t4_matched_genes_2024-10-29.csv")
  write.csv(merged_data, output_file, row.names = FALSE)
  
  message(paste("Successfully processed cluster C", i, sep = ""))
}
