
setwd("/Volumes/DataBox/MCS2023/Working Docs/Triplicated_chr16-17/")

# Load the CSV file (adjust the file path as needed)
regions <- read.csv("/Volumes/DataBox/MCS2023/Working Docs/Triplicated_chr16-17/Chr16 84.718.526-98.183.786 genes.csv")

# Inspect the first few rows of the CSV
head(regions)

# Rename columns for clarity (optional, but helpful)
colnames(regions) <- c("chromosome", "start", "end", "name")

# Subtract 1 from the start position for BED format (0-based start positions)
regions$start <- regions$start - 1

# Replace '/path_to_directory/' with the actual path where you want to save the file
directory_path <- "/Volumes/DataBox/MCS2023/Working Docs/Triplicated_chr16-17/"

# Specify the correct path
output_file <- file.path(directory_path, "chr16_output.bed")

# Save the BED file
write.table(regions[, c("chromosome", "start", "end")], 
            file = output_file, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

##########################
##########################


# Load the CSV file (adjust the file path as needed)
regions <- read.csv("/Volumes/DataBox/MCS2023/Working Docs/Triplicated_chr16-17/Chr17 1-9.426.821 genes.csv")

# Inspect the first few rows of the CSV
head(regions)

# Rename columns for clarity (optional, but helpful)
colnames(regions) <- c("chromosome", "start", "end", "name")

# Subtract 1 from the start position for BED format (0-based start positions)
regions$start <- regions$start - 1

# Replace '/path_to_directory/' with the actual path where you want to save the file
directory_path <- "/Volumes/DataBox/MCS2023/Working Docs/Triplicated_chr16-17/"

# Specify the correct path
output_file <- file.path(directory_path, "chr17_output.bed")

# Save the BED file
write.table(regions[, c("chromosome", "start", "end")], 
            file = output_file, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

