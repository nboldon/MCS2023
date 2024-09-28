## Copy files
# Define the target directory
target_directory="/project/eon/nboldon/MCS2023/IGV_Files"

# Adjust the top-level directory path as needed
top_level_directory="/project/eon/nboldon/MCS2023/RenameFiles"

# Use 'find' to locate files ending with .bam, .bai, or .fragments within the top level directory
find "$top_level_directory" -type f \( -name "*.bam" -o -name "*.bai" -o -name "*_fragments.tsv.gz" -o name "*_fragments.tsv.gz.tbi" \) -exec cp {} "$target_directory" \;

##########################################
##########################################

## To convert one file

conda activate pandas

# In python

import pandas as pd

# Read the gzipped fragment file, skipping lines starting with "#"
df = pd.read_csv("fragments.tsv.gz", delimiter='\t', comment='#', header=None, compression='gzip')

# Assign column names manually
df.columns = ['chromosome', 'start_position', 'end_position'] + [f'column_{i}' for i in range(3, len(df.columns))]

# Create a .bed file
bed_df = df[['chromosome', 'start_position', 'end_position']]
bed_df.to_csv("fragments.bed", sep='\t', header=False, index=False)


#########################################################
#########################################################

# To loop for all files in the folder

conda activate pandas

# In python

import os
import pandas as pd
import glob

# Define the directory where your gzipped fragments files are located
input_directory = '/project/eon/nboldon/MCS2023/IGV_Files/'

files = glob.glob(os.path.join(input_directory, '*fragments.tsv.gz'))

# Iterate over the list of files
for file in files:
    # Extract the basename of the original fragments file (without extension)
    base_name = os.path.splitext(os.path.basename(file))[0]
   
    # Read the gzipped fragment file, skipping lines starting with "#"
    df = pd.read_csv(file, delimiter='\t', comment='#', header=None, compression='gzip')
   
    # Assign column names manually
    df.columns = ['chromosome', 'start_position', 'end_position'] + [f'column_{i}' for i in range(3, len(df.columns))]
   
    # Create a .bed file with the same basename
    bed_df = df[['chromosome', 'start_position', 'end_position']]
    bed_file_name = os.path.join(input_directory, f'{base_name}.bed')
    bed_df.to_csv(bed_file_name, sep='\t', header=False, index=False)

