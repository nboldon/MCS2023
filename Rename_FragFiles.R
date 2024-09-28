This will update it so that it only finds and renames files within the outs and nested directories, not the outs and other directories themselves:


# Use 'find' to locate sample directories with names ending in "_cr"
# Adjust the top-level directory path as needed
top_level_directory="/project/eon/nboldon/MCS2023/RenameFiles"
sample_directories=$(find "$top_level_directory" -mindepth 1 -type d -name "*_cr")

# Loop through each sample directory
for sample_dir in $sample_directories; do
    # Extract the sample name without the "_cr" suffix
    sample_name=$(basename "$sample_dir" | sed 's/_cr$//')

    # Use 'find' to locate all files within "outs" and its nested subdirectories
    outs_contents=$(find "$sample_dir/outs" -type f)

    # Loop through each file within "outs"
    for item in $outs_contents; do
        # Extract the directory and item name
        dir=$(dirname "$item")
        item_name=$(basename "$item")

        # Construct the new item name with the prefix
        new_item_name="${sample_name}_${item_name}"

        # Rename the item within the same directory
        mv "$item" "$dir/$new_item_name"
    done
done

###################################################
###################################################

You'll need to re-copy the clean data over again to use the above. Or, with the prefixes attached to the outs directories, you can run the following, which accounts for variable prefixes to the beginning of each of each outs directory:

# Use 'find' to locate sample directories with names ending in "_cr"
# Adjust the top-level directory path as needed
top_level_directory="/project/eon/nboldon/MCS2023/RenameFiles"
sample_directories=$(find "$top_level_directory" -mindepth 1 -type d -name "*_cr")

# Loop through each sample directory
for sample_dir in $sample_directories; do
    # Extract the sample name without the "_cr" suffix
    sample_name=$(basename "$sample_dir" | sed 's/_cr$//')

    # Use 'find' to locate all "outs" directories within the sample directory
    outs_dirs=$(find "$sample_dir" -type d -name "*_outs")

    # Loop through each "outs" directory
    for outs_dir in $outs_dirs; do
        # Extract the name of the "outs" directory (the variable prefix)
        outs_name=$(basename "$outs_dir")

        # Use 'find' to locate all files within the "outs" and its nested subdirectories
        outs_contents=$(find "$outs_dir" -type f)

        # Loop through each file within "outs"
        for item in $outs_contents; do
            # Extract the directory and item name
            dir=$(dirname "$item")
            item_name=$(basename "$item")

            # Construct the new item name with the prefix and sample name
            new_item_name="${sample_name}_${outs_name}_${item_name}"

            # Rename the item within the same directory
            mv "$item" "$dir/$new_item_name"
        done
    done
done

######################################
######################################

##To move the renamed files

# Define the target directory
target_directory="/project/eon/nboldon/MCS2023/fragAnalysis_TSS10"

# Adjust the top-level directory path as needed
top_level_directory="/project/eon/nboldon/MCS2023/RenameFiles/FroC1_out"

# Use 'find' to locate files ending with .bam, .bai, or .fragments within the top level directory
find "$top_level_directory" -type f \( -name "*fragments.tsv.gz" -o -name "*.tsv.gz.tbi" \) -exec cp {} "$target_directory" \;

#Repeat for all directories in RenameFiles



Edit the extensions as you like, or add in any others, just adding in another -o -name before new ones. This gets anything at any nested level within the top level directory.

