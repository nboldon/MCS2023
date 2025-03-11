
# Set working directory to where the CSV files are located
setwd("/Volumes/DataBox/ProjMCS7/fragCounts")

######################################################
######################################################
######################################################

## Combine individual files into one and normalize frag counts
# Normalize by dividing the number of gene frags by total number of frags in sample


# List of fragment count files
fragment_files <- list.files(pattern = "fragment_counts.*\\.csv")

# Read the sample statistics file
sample_stats <- read.csv("sample_statistics.csv")

# Loop through each fragment count file and merge it with sample_stats
for (file in fragment_files) {
  # Read the fragment counts file
  fragment_counts <- read.csv(file)
  
  # Get the column name for fragment counts
  colname <- names(fragment_counts)[2] # Assuming the second column holds the fragment counts
  
  # Merge the fragment counts with the sample statistics
  sample_stats <- merge(sample_stats, fragment_counts, by = "Sample", all.x = TRUE)
  
  # Calculate the ratio for the current fragment counts and add a new column
  ratio_colname <- paste0(colname, "_Ratio")
  sample_stats[[ratio_colname]] <- sample_stats[[colname]] / sample_stats$TotalFragments
}

# Save the final merged data to a new CSV file
write.csv(sample_stats, "combined_sample_statistics.csv", row.names = FALSE)


######################################################
######################################################
######################################################


