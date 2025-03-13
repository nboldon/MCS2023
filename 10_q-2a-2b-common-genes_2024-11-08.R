

########################################################
########################################################
########################################################


## Determine similar genes Question 2a and 2b methods for treatment comps by cluster 
## To loop through all treatment groups and all clusters
## Also adds a gene count tally for each treatment group by cluster


# Set working directory
setwd("/Volumes/DataBox/Heatmap_Comps")

# Define the treatment groups and clusters
treatment_groups <- c("T1", "T2", "T3", "T4")
clusters <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9",
              "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19",
              "C20", "C21", "C22", "C23", "C24", "C25")

# Initialize an empty data frame to store the results
final_results <- data.frame(
  Cluster = character(),
  GeneName = character(),
  Log2FC_T1 = numeric(),
  ZScore_T1 = numeric(),
  Log2FC_T2 = numeric(),
  ZScore_T2 = numeric(),
  Log2FC_T3 = numeric(),
  ZScore_T3 = numeric(),
  Log2FC_T4 = numeric(),
  ZScore_T4 = numeric(),
  stringsAsFactors = FALSE
)

# Initialize a data frame to track the tally of common genes by treatment and cluster
gene_tally <- data.frame(
  Cluster = character(),
  Treatment = character(),
  CommonGeneCount = integer(),
  stringsAsFactors = FALSE
)

# Load necessary libraries
library(dplyr)

# Loop through each cluster
for (cluster in clusters) {
  
  # Initialize `cluster_results` as NULL to avoid row constraints
  cluster_results <- NULL
  
  # Loop through each treatment group (T1-T4) for the current cluster
  for (treatment in treatment_groups) {
    
    # Construct file paths for gene markers and Z-scores specific to each cluster and treatment
    gene_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", treatment, "_", cluster, "_GeneMarkers_2024-04-10.csv")
    zscore_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", cluster, "_byTx_zscores_2024-03-21.csv")
    
    # Use tryCatch to handle any potential errors (e.g., missing files or data)
    tryCatch({
      
      # Read the gene markers for the current treatment and cluster
      T_gene_markers <- read.csv(gene_file)
      
      # Check if `name` column is present in gene marker file
      if (!"name" %in% colnames(T_gene_markers)) {
        message(paste("Column 'name' not found in gene markers file for", treatment, "in", cluster))
        next
      }
      
      # Read the Z-scores file (specific to each cluster)
      C3_byTx_zscores <- read.csv(zscore_file, row.names = 1)
      
      # Extract gene names and find common genes
      genes_T <- T_gene_markers$name
      genes_C3 <- rownames(C3_byTx_zscores)
      common_genes <- intersect(genes_T, genes_C3)
      
      # Safety check: if no common genes, skip the iteration
      if (length(common_genes) == 0) {
        message(paste("No common genes found for treatment", treatment, "in cluster", cluster))
        next
      }
      
      # Extract Log2FC and Z-scores for common genes
      T_common_genes <- T_gene_markers %>% filter(name %in% common_genes)
      treatment_col <- tolower(treatment)
      
      # Check if the treatment column exists in the z-scores data
      if (!(treatment_col %in% colnames(C3_byTx_zscores))) {
        message(paste("Column", treatment_col, "not found in Z-score file for cluster", cluster))
        next
      }
      
      C3_common_genes <- C3_byTx_zscores[common_genes, treatment_col, drop = FALSE]
      
      # Combine the data for the common genes into a single data frame
      combined_data <- data.frame(
        GeneName = T_common_genes$name,
        Cluster = cluster,
        Log2FC = T_common_genes$Log2FC,
        ZScore = C3_common_genes[[1]]
      )
      
      # Rename columns for this specific treatment
      colnames(combined_data)[3] <- paste0("Log2FC_", treatment)
      colnames(combined_data)[4] <- paste0("ZScore_", treatment)
      
      # Ensure combined_data is not empty before joining
      if (!is.null(cluster_results)) {
        cluster_results <- full_join(cluster_results, combined_data, by = c("GeneName", "Cluster"))
      } else {
        cluster_results <- combined_data
      }
      
      # Track the tally of common genes by treatment and cluster
      gene_tally <- rbind(gene_tally, data.frame(
        Cluster = cluster,
        Treatment = treatment,
        CommonGeneCount = length(common_genes)
      ))
      
    }, error = function(e) {
      # Handle any errors and continue to the next treatment if an error occurs
      message(paste("Error processing treatment", treatment, "in cluster", cluster, ":", e$message))
    })
  }
  
  # Append cluster-specific results to the final results data frame
  final_results <- bind_rows(final_results, cluster_results)
}

# View the final results
print(final_results)

# Save the final results to a .csv file
write.csv(final_results, "q-2a-2b_common-genes_2024-11-08.csv", row.names = FALSE)

# Save the gene tally to a .csv file
write.csv(gene_tally, "q-2a-2b_gene_tally_by_treatment_cluster_2024-11-08.csv", row.names = FALSE)



# No common genes found for treatment T2 in cluster C2
# No common genes found for treatment T1 in cluster C5
# No common genes found for treatment T2 in cluster C5
# No common genes found for treatment T3 in cluster C5
# No common genes found for treatment T4 in cluster C5
# No common genes found for treatment T1 in cluster C7
# No common genes found for treatment T2 in cluster C7
# No common genes found for treatment T3 in cluster C7
# No common genes found for treatment T4 in cluster C7
# No common genes found for treatment T1 in cluster C9
# No common genes found for treatment T2 in cluster C9
# No common genes found for treatment T3 in cluster C9
# No common genes found for treatment T4 in cluster C9
# No common genes found for treatment T1 in cluster C17


########################################################
########################################################
########################################################
########################################################
########################################################
########################################################


## Barplot for gene tally


# Set working directory
setwd("/Volumes/DataBox/Heatmap_Comps")

# Load necessary libraries
library(ggplot2)
library(viridis)
library(forcats)

# Read the data
gene_tally <- read.csv("q-2a-2b_gene_tally_by_treatment_cluster_2024-11-08.csv")

# Create a custom sorting function for clusters
sort_clusters <- function(x) {
  # Extract numbers from cluster names and sort based on those
  nums <- as.numeric(gsub("C", "", x))
  x[order(nums)]
}

# Convert Cluster to factor with custom ordering
gene_tally$Cluster <- factor(gene_tally$Cluster, 
                             levels = sort_clusters(unique(gene_tally$Cluster)),
                             ordered = TRUE)

# Use the viridis color palette
color_palette <- viridis(4)  # Adjust the number if you have more treatments

# Create a barplot with increased width
ggplot(gene_tally, aes(x = Cluster, y = CommonGeneCount, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Tally of Common Genes by Treatment Group and Cluster",
    x = "Cluster",
    y = "Number of Common Genes"
  ) +
  theme_bw() +  # Changed from theme_minimal() to theme_bw()
  theme(
    panel.background = element_rect(fill = "white"),  # White background
    plot.background = element_rect(fill = "white"),   # White plot background
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white"), # White legend background
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save the plot with increased width
ggsave("q-2a-2b_gene_tally_by_tx-Cluster_Barplot_2024-11-08.png", 
       width = 12,
       height = 7,
       dpi = 300,
       bg = "white")  # Added white background for saved plot



########################################################
########################################################
########################################################
########################################################
########################################################
########################################################


## To create a new spreadsheet of gene names based on results from above comparisons


# Initialize a data frame to store the common genes across clusters by treatment
common_genes_across_clusters <- data.frame(
  Treatment = character(),
  GeneName = character(),
  ClustersSharedIn = character(),
  stringsAsFactors = FALSE
)

# Loop through each treatment group to find common genes between any two clusters
for (treatment in treatment_groups) {
  
  # Initialize a list to store genes by cluster for the current treatment
  genes_by_cluster <- list()
  
  # Loop through each cluster for the current treatment
  for (cluster in clusters) {
    
    # Construct file paths
    gene_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", treatment, "_", cluster, "_GeneMarkers_2024-04-10.csv")
    zscore_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", cluster, "_byTx_zscores_2024-03-21.csv")
    
    # Read gene markers and Z-scores files, handle errors
    tryCatch({
      if (file.exists(gene_file) && file.exists(zscore_file)) {
        # Read the gene markers and Z-scores files
        T_gene_markers <- read.csv(gene_file)
        C3_byTx_zscores <- read.csv(zscore_file, row.names = 1)
        
        # Extract gene names and find common genes
        genes_T <- T_gene_markers$name
        genes_C3 <- rownames(C3_byTx_zscores)
        
        # Debugging: Check gene names in both files
        message(paste("Cluster:", cluster, "| Treatment:", treatment))
        message("Genes in gene file:", paste(genes_T[1:5], collapse = ", "))
        message("Genes in zscore file:", paste(genes_C3[1:5], collapse = ", "))
        
        # Find common genes between gene markers and Z-scores
        common_genes <- intersect(genes_T, genes_C3)
        
        # Debugging: Display common genes found
        if (length(common_genes) > 0) {
          message("Common genes found:", paste(common_genes[1:5], collapse = ", "))
          genes_by_cluster[[cluster]] <- common_genes
        } else {
          message(paste("No common genes found for treatment", treatment, "in cluster", cluster))
        }
        
      } else {
        message(paste("File missing for treatment", treatment, "in cluster", cluster))
      }
      
    }, error = function(e) {
      message(paste("Error processing treatment", treatment, "in cluster", cluster, ":", e$message))
    })
  }
  
  # Check if genes_by_cluster has at least 2 clusters with data for the current treatment
  if (length(genes_by_cluster) > 1) {
    
    # Compare genes between every pair of clusters
    for (cluster1 in names(genes_by_cluster)) {
      for (cluster2 in names(genes_by_cluster)) {
        
        # Skip comparing the same cluster with itself
        if (cluster1 == cluster2) next
        
        # Find common genes between the two clusters
        common_in_pair <- intersect(genes_by_cluster[[cluster1]], genes_by_cluster[[cluster2]])
        
        # If there are common genes, store the result
        if (length(common_in_pair) > 0) {
          common_genes_across_clusters <- rbind(
            common_genes_across_clusters,
            data.frame(
              Treatment = treatment,
              GeneName = common_in_pair,
              ClustersSharedIn = paste(cluster1, cluster2, sep = ", "),
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  } else {
    message(paste("Insufficient clusters with data for treatment", treatment))
  }
}

# Save the common genes across clusters to a .csv file
write.csv(common_genes_across_clusters, "common_genes_between_clusters_by_treatment_2024-11-08.csv", row.names = FALSE)

# View the final results
print(common_genes_across_clusters)


############################################################


## To create a heatmap of the above outputs:

# Set working directory
setwd("/Volumes/DataBox/Heatmap_Comps")

# Define the treatment groups and clusters
treatment_groups <- c("T1", "T2", "T3", "T4")
clusters <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", 
              "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", 
              "C20", "C21", "C22", "C23", "C24", "C25")

# Initialize an empty data frame to store the results
final_results <- data.frame(
  Cluster = character(),
  GeneName = character(),
  Log2FC_T1 = numeric(),
  ZScore_T1 = numeric(),
  Log2FC_T2 = numeric(),
  ZScore_T2 = numeric(),
  Log2FC_T3 = numeric(),
  ZScore_T3 = numeric(),
  Log2FC_T4 = numeric(),
  ZScore_T4 = numeric(),
  stringsAsFactors = FALSE
)

# Initialize a data frame to track the tally of common genes by treatment and cluster
gene_tally <- data.frame(
  Cluster = character(),
  Treatment = character(),
  CommonGeneCount = integer(),
  stringsAsFactors = FALSE
)

# Initialize an empty list to store the common genes across clusters for later use in the gene matrix
common_genes_across_clusters <- list()

# Loop through each cluster
for (cluster in clusters) {
  
  # Initialize `cluster_results` as NULL to avoid row constraints
  cluster_results <- NULL
  
  # Loop through each treatment group (T1-T4) for the current cluster
  for (treatment in treatment_groups) {
    
    # Construct file paths for gene markers and Z-scores specific to each cluster and treatment
    gene_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", treatment, "_", cluster, "_GeneMarkers_2024-04-10.csv")
    zscore_file <- paste0("/Volumes/DataBox/MCS2023/Tx_Comp/", cluster, "_byTx_zscores_2024-03-21.csv")
    
    # Use tryCatch to handle any potential errors (e.g., missing files or data)
    tryCatch({
      
      # Read the gene markers for the current treatment and cluster
      T_gene_markers <- read.csv(gene_file)
      
      # Check if `name` column is present in gene marker file
      if (!"name" %in% colnames(T_gene_markers)) {
        message(paste("Column 'name' not found in gene markers file for", treatment, "in", cluster))
        next
      }
      
      # Read the Z-scores file (specific to each cluster)
      C3_byTx_zscores <- read.csv(zscore_file, row.names = 1)
      
      # Extract gene names and find common genes
      genes_T <- T_gene_markers$name
      genes_C3 <- rownames(C3_byTx_zscores)
      common_genes <- intersect(genes_T, genes_C3)
      
      # Safety check: if no common genes, skip the iteration
      if (length(common_genes) == 0) {
        message(paste("No common genes found for treatment", treatment, "in cluster", cluster))
        next
      }
      
      # Extract Log2FC and Z-scores for common genes
      T_common_genes <- T_gene_markers %>% filter(name %in% common_genes)
      treatment_col <- tolower(treatment)
      
      # Check if the treatment column exists in the z-scores data
      if (!(treatment_col %in% colnames(C3_byTx_zscores))) {
        message(paste("Column", treatment_col, "not found in Z-score file for cluster", cluster))
        next
      }
      
      C3_common_genes <- C3_byTx_zscores[common_genes, treatment_col, drop = FALSE]
      
      # Combine the data for the common genes into a single data frame
      combined_data <- data.frame(
        GeneName = T_common_genes$name,
        Cluster = cluster,
        Log2FC = T_common_genes$Log2FC,
        ZScore = C3_common_genes[[1]]
      )
      
      # Rename columns for this specific treatment
      colnames(combined_data)[3] <- paste0("Log2FC_", treatment)
      colnames(combined_data)[4] <- paste0("ZScore_", treatment)
      
      # Ensure combined_data is not empty before joining
      if (!is.null(cluster_results)) {
        cluster_results <- full_join(cluster_results, combined_data, by = c("GeneName", "Cluster"))
      } else {
        cluster_results <- combined_data
      }
      
      # Track the tally of common genes by treatment and cluster
      gene_tally <- rbind(gene_tally, data.frame(
        Cluster = cluster,
        Treatment = treatment,
        CommonGeneCount = length(common_genes)
      ))
      
      # Store the common genes information in the list for later use in the gene matrix
      common_genes_across_clusters <- append(common_genes_across_clusters, list(data.frame(
        Treatment = treatment,
        ClustersSharedIn = cluster,  # Make sure 'ClustersSharedIn' is the correct column name
        CommonGeneCount = length(common_genes)
      )))
      
    }, error = function(e) {
      # Handle any errors and continue to the next treatment if an error occurs
      message(paste("Error processing treatment", treatment, "in cluster", cluster, ":", e$message))
    })
  }
  
  # Append cluster-specific results to the final results data frame
  final_results <- bind_rows(final_results, cluster_results)
}


# Step 1: Aggregate the common gene tally data into a matrix format
gene_matrix <- gene_tally %>%
  pivot_wider(names_from = Treatment, values_from = CommonGeneCount, values_fill = list(CommonGeneCount = 0))

# View the resulting gene matrix
head(gene_matrix)

# Step 2: You may want to clean up the data by ensuring no missing values or performing any transformations.
# For example, if you want to filter out clusters with no common genes:
gene_matrix_filtered <- gene_matrix %>%
  filter(T1 > 0 | T2 > 0 | T3 > 0 | T4 > 0)

# Step 3: Save the final results to CSV
write.csv(gene_matrix_filtered, "/Volumes/DataBox/Heatmap_Comps/q-2a-2b_final_gene_matrix_2024-11-08.csv", row.names = FALSE)

# Step 4: Optionally, if you need to create a heatmap from this data, you could do so using `pheatmap`
library(pheatmap)
library(tibble)

# Convert the 'Cluster' column to rownames
gene_matrix_heatmap <- gene_matrix_filtered %>% 
  column_to_rownames("Cluster")

# Now you can generate the heatmap
pheatmap(gene_matrix_heatmap)

# Step 5: Save the heatmap plot as an image (optional)
ggsave("/Volumes/DataBox/Heatmap_Comps/q-2a-2b_gene_tally_by_tx-Cluster_Heatmap_2024-11-08.png")



############################################################



## output from above code: 

Cluster: C1 | Treatment: T1
Genes in gene file:Fgfr2, Ppfibp2, Nrg3, Ksr2, Myo1d
Genes in zscore file:Fam107a, Celsr1, 5930403L14Rik, 4930544G11Rik, Vmn1r174
Common genes found:Fgfr2, Fam107a, Gm10863, Olig1, Prdm16
Cluster: C2 | Treatment: T1
Genes in gene file:Tns3, Gm10863, A230009B12Rik, Pacrg, Kcnb1
Genes in zscore file:Hfe, Dennd1c, Lims2, Mmp2, Prps1l1
Common genes found:Tns3, Gm10863, Sox10, NA, NA
Cluster: C3 | Treatment: T1
Genes in gene file:A230009B12Rik, Kalrn, Gm10863, Rnf220, Pacrg
Genes in zscore file:2810429I04Rik, Acod1, Mucl1, S100a16, Tshb
Common genes found:A230009B12Rik, Gm10863, Rnf220, Fgfr2, Fa2h
Cluster: C4 | Treatment: T1
Genes in gene file:Cldn11, Mal, Nkain1, Prima1, Gpr37
Genes in zscore file:Olfr1040, Fa2h, Fcrl6, Olfr1214, Gm7978
Common genes found:Cldn11, A530053G22Rik, Ppp1r14a, Gng11, Folh1
Cluster: C5 | Treatment: T1
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Lbh, Snai1, Epas1, 4933440J02Rik, Foxn3
No common genes found for treatment T1 in cluster C5
Cluster: C6 | Treatment: T1
Genes in gene file:Cacng4, Olig2, Olig1, Xylt1, Lhfpl3
Genes in zscore file:Olfr324, Tmem100, Gm6961, Unc13c, Plac1
Common genes found:Cacng4, Olig2, Olig1, Lhfpl3, Lims2
Cluster: C7 | Treatment: T1
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Tcf7l1, Mir26b, Abca8b, Bmp6, Ppara
No common genes found for treatment T1 in cluster C7
Cluster: C8 | Treatment: T1
Genes in gene file:Prdm16, Msi2, Slc1a2, Fgfr3, Gm13872
Genes in zscore file:Gpr179, Ranbp3l, Prodh, Frmpd1os, Lfng
Common genes found:Prdm16, Msi2, Slc1a2, Fgfr3, Gm13872
Cluster: C9 | Treatment: T1
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Gm15441, Loxl3, Dok1, Mir7040, Gm21057
No common genes found for treatment T1 in cluster C9
Cluster: C10 | Treatment: T1
Genes in gene file:Csf1r, Zfp36l2, Fli1, Cx3cr1, Rbm47
Genes in zscore file:Olfr433, Lgals9, Slfn3, Olfr1242, B020014A21Rik
Common genes found:Csf1r, Fli1, Cx3cr1, Rbm47, Vsir
Cluster: C11 | Treatment: T1
Genes in gene file:Cdh23, Cx3cr1, Vsir, Csf1r, 4933433H22Rik
Genes in zscore file:Olfr784, Susd3, Ccr5, Dusp27, Havcr2
Common genes found:Cdh23, Cx3cr1, Vsir, Csf1r, E230029C05Rik
Cluster: C12 | Treatment: T1
Genes in gene file:Ebf1, Tns1, Rgs5, Zic4, Cald1
Genes in zscore file:Adap2, Slc38a11, Olfr1000, Mylk2, Rarres2
Common genes found:Ebf1, Rgs5, Zic4, Cald1, Pald1
Cluster: C13 | Treatment: T1
Genes in gene file:2310005E17Rik, Zic1, Nr2f2, Mylk4, Slc6a13
Genes in zscore file:Iigp1, Pcolce, Mrgprb2, Jaml, Vmn2r114
Common genes found:2310005E17Rik, Zic1, Nr2f2, Mylk4, Slc6a13
Cluster: C14 | Treatment: T1
Genes in gene file:2310005E17Rik, Zic2, Cgnl1, Tgfbr3, Rbfox3
Genes in zscore file:Olfr1490, Tas2r107, Olfr508, Samt4, Olfr384
Common genes found:2310005E17Rik, Cgnl1, Zic4, Eya2, 1700018B08Rik
Cluster: C15 | Treatment: T1
Genes in gene file:Cdh23, Vsir, Cx3cr1, Cd33, E230029C05Rik
Genes in zscore file:Ccl6, Tnfrsf17, Clec4a3, Sash3, Batf3
Common genes found:Vsir, Cx3cr1, Cd33, E230029C05Rik, Ccr1
Cluster: C16 | Treatment: T1
Genes in gene file:Tshz2, AY702102, Spon1, Gm8630, Abca12
Genes in zscore file:Pla2g2d, Vsig1, Gm22650, Trim58, A630019I02Rik
Common genes found:Tshz2, AY702102, Spon1, Gm8630, Abca12
Cluster: C17 | Treatment: T1
Genes in gene file:Gm14204, Slc32a1, Olig1, 8430430B14Rik, Kcnj10
Genes in zscore file:Olfr328, Hba-a2, Fam240b, Dlx6os2, Omt2a
No common genes found for treatment T1 in cluster C17
Cluster: C18 | Treatment: T1
Genes in gene file:Lingo1, Unc5d, Car10, Vwa5b1, Ubxn10
Genes in zscore file:Krtap1-5, Krtap1-4, Krtap4-2, Retnlb, Dact2
Common genes found:Lingo1, Unc5d, Vwa5b1, Ubxn10, Rtn4r
Cluster: C19 | Treatment: T1
Genes in gene file:Tshz2, Itprid1, Gm8267, Adcy10, Htr2c
Genes in zscore file:Aard, Htr4, Htr7, B430306N03Rik, 4930546C10Rik
Common genes found:Tshz2, Itprid1, Gm8267, Adcy10, Htr2c
Cluster: C20 | Treatment: T1
Genes in gene file:4933412E24Rik, 1700121N20Rik, Fam19a1, Bcl11b, Ldb2
Genes in zscore file:Dnah9, Acrv1, Olfr887, Olfr1395, Vmn1r-ps103
Common genes found:4933412E24Rik, 1700121N20Rik, Fam19a1, Bcl11b, Vat1l
Cluster: C21 | Treatment: T1
Genes in gene file:Hs3st4, Foxp2, Zfpm2, Gm8630, Syt6
Genes in zscore file:Igfbp4, Lrrc15, Crym, Daw1, Gal3st2
Common genes found:Hs3st4, Foxp2, Zfpm2, Gm8630, Syt6
Cluster: C22 | Treatment: T1
Genes in gene file:Sox6, Dlgap2, Nxph1, Grip1, Ptprm
Genes in zscore file:Olfr282, Klhl14, Gm572, Bend4, Sox6
Common genes found:Sox6, Nxph1, Grip1, Ptprm, Kcnmb2
Cluster: C23 | Treatment: T1
Genes in gene file:Adarb2, Erbb4, Dlx6os1, Zfp536, Kcnip1
Genes in zscore file:4930444F02Rik, Taar7f, Olfr1384, Olfr1383, Spp2
Common genes found:Adarb2, Erbb4, Dlx6os1, Zfp536, Kcnip1
Cluster: C24 | Treatment: T1
Genes in gene file:Prdm16, Slc1a2, Gm35978, Notch1, Fgfr3
Genes in zscore file:Gfap, 4930539C22Rik, Vmn1r64, Vmn1r78, Usp17le
Common genes found:Prdm16, Fgfr3, Gm13872, 2900052N01Rik, 4930550C17Rik
Cluster: C25 | Treatment: T1
Genes in gene file:Gpr17, Pdgfra, Lims2, Cacng4, Lhfpl3
Genes in zscore file:Ttr, Rpl13, Snord68, Inava, C1ql1
Common genes found:Gpr17, Pdgfra, Lims2, Cacng4, Lhfpl3
Cluster: C1 | Treatment: T2
Genes in gene file:Nat8, Nat8f2, Erc2, A230009B12Rik, Slc1a2
Genes in zscore file:Fam107a, Celsr1, 5930403L14Rik, 4930544G11Rik, Vmn1r174
Common genes found:Nat8, Nat8f2, Fgfr2, Prdm16, Kcnj10
Cluster: C2 | Treatment: T2
Genes in gene file:Syt1, NA, NA, NA, NA
Genes in zscore file:Hfe, Dennd1c, Lims2, Mmp2, Prps1l1
No common genes found for treatment T2 in cluster C2
Cluster: C3 | Treatment: T2
Genes in gene file:A230009B12Rik, Rnf220, Gm10863, Kalrn, Pacrg
Genes in zscore file:2810429I04Rik, Acod1, Mucl1, S100a16, Tshb
Common genes found:A230009B12Rik, Rnf220, Gm10863, Fgfr2, Mal
Cluster: C4 | Treatment: T2
Genes in gene file:Sox10, Tmem63a, Mal, Gng11, Lims2
Genes in zscore file:Olfr1040, Fa2h, Fcrl6, Olfr1214, Gm7978
Common genes found:Gng11, Rnase1, Gm1979, Rpl10l, NA
Cluster: C5 | Treatment: T2
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Lbh, Snai1, Epas1, 4933440J02Rik, Foxn3
No common genes found for treatment T2 in cluster C5
Cluster: C6 | Treatment: T2
Genes in gene file:Olig2, Olig1, Cacng4, Xylt1, Kank1
Genes in zscore file:Olfr324, Tmem100, Gm6961, Unc13c, Plac1
Common genes found:Olig2, Olig1, Cacng4, Kank1, Pdgfra
Cluster: C7 | Treatment: T2
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Tcf7l1, Mir26b, Abca8b, Bmp6, Ppara
No common genes found for treatment T2 in cluster C7
Cluster: C8 | Treatment: T2
Genes in gene file:Prdm16, Msi2, Slc1a2, Fgfr3, 0610039H22Rik
Genes in zscore file:Gpr179, Ranbp3l, Prodh, Frmpd1os, Lfng
Common genes found:Prdm16, Msi2, Slc1a2, Fgfr3, 0610039H22Rik
Cluster: C9 | Treatment: T2
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Gm15441, Loxl3, Dok1, Mir7040, Gm21057
No common genes found for treatment T2 in cluster C9
Cluster: C10 | Treatment: T2
Genes in gene file:Runx1, Lyn, Fli1, Cdh23, Cx3cr1
Genes in zscore file:Olfr433, Lgals9, Slfn3, Olfr1242, B020014A21Rik
Common genes found:Fli1, Cx3cr1, Ccr6, Vsir, Tnfrsf1b
Cluster: C11 | Treatment: T2
Genes in gene file:Cdh23, Vsir, Cx3cr1, Csf1r, E230029C05Rik
Genes in zscore file:Olfr784, Susd3, Ccr5, Dusp27, Havcr2
Common genes found:Cdh23, Vsir, Cx3cr1, Csf1r, E230029C05Rik
Cluster: C12 | Treatment: T2
Genes in gene file:Atp13a5, Pald1, Pdgfrb, Rbpms, Ets1
Genes in zscore file:Adap2, Slc38a11, Olfr1000, Mylk2, Rarres2
Common genes found:Atp13a5, Pald1, Pdgfrb, Rbpms, Ebf1
Cluster: C13 | Treatment: T2
Genes in gene file:Uaca, Angptl8, Zic4, 2310005E17Rik, Zic1
Genes in zscore file:Iigp1, Pcolce, Mrgprb2, Jaml, Vmn2r114
Common genes found:Uaca, Zic4, 2310005E17Rik, Zic1, Cgnl1
Cluster: C14 | Treatment: T2
Genes in gene file:2310005E17Rik, Atp2b2, Cgnl1, Zic1, Tgfbr3
Genes in zscore file:Olfr1490, Tas2r107, Olfr508, Samt4, Olfr384
Common genes found:2310005E17Rik, Cgnl1, Zic1, Uaca, Zic4
Cluster: C15 | Treatment: T2
Genes in gene file:Csf1r, Cx3cr1, Fli1, Vsir, E230029C05Rik
Genes in zscore file:Ccl6, Tnfrsf17, Clec4a3, Sash3, Batf3
Common genes found:Csf1r, Cx3cr1, Fli1, Vsir, E230029C05Rik
Cluster: C16 | Treatment: T2
Genes in gene file:Tshz2, AY702102, Hs3st4, Abca12, Abcc12
Genes in zscore file:Pla2g2d, Vsig1, Gm22650, Trim58, A630019I02Rik
Common genes found:Tshz2, AY702102, Hs3st4, Abca12, Abcc12
Cluster: C17 | Treatment: T2
Genes in gene file:Slc32a1, Fam240b, Olig1, Gad2, Kcnj10
Genes in zscore file:Olfr328, Hba-a2, Fam240b, Dlx6os2, Omt2a
Common genes found:Fam240b, Gad2, NA, NA, NA
Cluster: C18 | Treatment: T2
Genes in gene file:Unc5d, D430041D05Rik, Nell2, Mir466k, Lingo1
Genes in zscore file:Krtap1-5, Krtap1-4, Krtap4-2, Retnlb, Dact2
Common genes found:Unc5d, Lingo1, Vwa5b1, Rtn4r, Dact2
Cluster: C19 | Treatment: T2
Genes in gene file:Tshz2, Gm8267, Adcy10, Gm32141, Itprid1
Genes in zscore file:Aard, Htr4, Htr7, B430306N03Rik, 4930546C10Rik
Common genes found:Tshz2, Gm8267, Adcy10, Gm32141, Itprid1
Cluster: C20 | Treatment: T2
Genes in gene file:4933412E24Rik, Fam19a1, 1700121N20Rik, Fras1, A330093E20Rik
Genes in zscore file:Dnah9, Acrv1, Olfr887, Olfr1395, Vmn1r-ps103
Common genes found:4933412E24Rik, Fam19a1, 1700121N20Rik, Fras1, Vat1l
Cluster: C21 | Treatment: T2
Genes in gene file:Hs3st4, Foxp2, Thsd7b, Zfpm2, Gm8630
Genes in zscore file:Igfbp4, Lrrc15, Crym, Daw1, Gal3st2
Common genes found:Hs3st4, Foxp2, Thsd7b, Zfpm2, Gm8630
Cluster: C22 | Treatment: T2
Genes in gene file:Sox6, Dlgap2, Nxph1, Kcnmb2, Ptprm
Genes in zscore file:Olfr282, Klhl14, Gm572, Bend4, Sox6
Common genes found:Sox6, Nxph1, Kcnmb2, Ptprm, Btbd11
Cluster: C23 | Treatment: T2
Genes in gene file:Adarb2, Erbb4, Zfp536, Dlx6os1, Btbd11
Genes in zscore file:4930444F02Rik, Taar7f, Olfr1384, Olfr1383, Spp2
Common genes found:Adarb2, Erbb4, Zfp536, Dlx6os1, Btbd11
Cluster: C24 | Treatment: T2
Genes in gene file:Fgfr3, 2900052N01Rik, Slc38a3, 5930403L14Rik, Grin2c
Genes in zscore file:Gfap, 4930539C22Rik, Vmn1r64, Vmn1r78, Usp17le
Common genes found:Fgfr3, 2900052N01Rik, 5930403L14Rik, AA619741, Ndrg2
Cluster: C25 | Treatment: T2
Genes in gene file:Olig2, Gpr17, Xylt1, Lims2, Cacng4
Genes in zscore file:Ttr, Rpl13, Snord68, Inava, C1ql1
Common genes found:Gpr17, Lims2, Cacng4, Myt1, C1ql1
Cluster: C1 | Treatment: T3
Genes in gene file:F630040K05Rik, Prdm16, Gab1, Erc2, Fgfr3
Genes in zscore file:Fam107a, Celsr1, 5930403L14Rik, 4930544G11Rik, Vmn1r174
Common genes found:Prdm16, Kcnj10, Plpp3, Eva1a, 3930402G23Rik
Cluster: C2 | Treatment: T3
Genes in gene file:Gabbr2, Tns3, Frmd4b, Tgfa, NA
Genes in zscore file:Hfe, Dennd1c, Lims2, Mmp2, Prps1l1
Common genes found:Tns3, Frmd4b, NA, NA, NA
Cluster: C3 | Treatment: T3
Genes in gene file:Gm10863, Kalrn, Rnf220, A230009B12Rik, Olig1
Genes in zscore file:2810429I04Rik, Acod1, Mucl1, S100a16, Tshb
Common genes found:Gm10863, Rnf220, A230009B12Rik, Olig1, Fgfr2
Cluster: C4 | Treatment: T3
Genes in gene file:Gm10863, Fa2h, 5430431A17Rik, Mal, Ppp1r14a
Genes in zscore file:Olfr1040, Fa2h, Fcrl6, Olfr1214, Gm7978
Common genes found:Fa2h, Ppp1r14a, Gjc3, 5430425K12Rik, A530053G22Rik
Cluster: C5 | Treatment: T3
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Lbh, Snai1, Epas1, 4933440J02Rik, Foxn3
No common genes found for treatment T3 in cluster C5
Cluster: C6 | Treatment: T3
Genes in gene file:Olig2, Cacng4, Pdgfra, Olig1, Xylt1
Genes in zscore file:Olfr324, Tmem100, Gm6961, Unc13c, Plac1
Common genes found:Olig2, Cacng4, Pdgfra, Olig1, Lims2
Cluster: C7 | Treatment: T3
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Tcf7l1, Mir26b, Abca8b, Bmp6, Ppara
No common genes found for treatment T3 in cluster C7
Cluster: C8 | Treatment: T3
Genes in gene file:Prdm16, Msi2, Slc1a2, Fgfr3, 0610039H22Rik
Genes in zscore file:Gpr179, Ranbp3l, Prodh, Frmpd1os, Lfng
Common genes found:Prdm16, Msi2, Slc1a2, Fgfr3, 0610039H22Rik
Cluster: C9 | Treatment: T3
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Gm15441, Loxl3, Dok1, Mir7040, Gm21057
No common genes found for treatment T3 in cluster C9
Cluster: C10 | Treatment: T3
Genes in gene file:Cx3cr1, Cdh23, Tmem119, Rgs10, Slco2b1
Genes in zscore file:Olfr433, Lgals9, Slfn3, Olfr1242, B020014A21Rik
Common genes found:Cx3cr1, Tmem119, Rgs10, Vsir, Fli1
Cluster: C11 | Treatment: T3
Genes in gene file:Cdh23, E230029C05Rik, Cx3cr1, Kalrn, 4933433H22Rik
Genes in zscore file:Olfr784, Susd3, Ccr5, Dusp27, Havcr2
Common genes found:Cdh23, E230029C05Rik, Cx3cr1, Ccr6, Vsir
Cluster: C12 | Treatment: T3
Genes in gene file:Pdgfrb, Carmn, Pald1, Ebf1, Rbpms
Genes in zscore file:Adap2, Slc38a11, Olfr1000, Mylk2, Rarres2
Common genes found:Pdgfrb, Carmn, Pald1, Ebf1, Rbpms
Cluster: C13 | Treatment: T3
Genes in gene file:Zic4, 2310005E17Rik, Zic1, Rbpms, Phldb2
Genes in zscore file:Iigp1, Pcolce, Mrgprb2, Jaml, Vmn2r114
Common genes found:Zic4, 2310005E17Rik, Zic1, Rbpms, Foxc1
Cluster: C14 | Treatment: T3
Genes in gene file:Tgfbr3, 2310005E17Rik, Zic2, Rbfox3, Zic4
Genes in zscore file:Olfr1490, Tas2r107, Olfr508, Samt4, Olfr384
Common genes found:2310005E17Rik, Zic4, Uaca, Eya2, Cgnl1
Cluster: C15 | Treatment: T3
Genes in gene file:Ccr6, Cx3cr1, Vsir, Siglech, Csf1r
Genes in zscore file:Ccl6, Tnfrsf17, Clec4a3, Sash3, Batf3
Common genes found:Ccr6, Cx3cr1, Vsir, Siglech, Csf1r
Cluster: C16 | Treatment: T3
Genes in gene file:Tshz2, AY702102, Hs3st4, Etl4, Grik3
Genes in zscore file:Pla2g2d, Vsig1, Gm22650, Trim58, A630019I02Rik
Common genes found:Tshz2, AY702102, Hs3st4, Grik3, Abcc12
Cluster: C17 | Treatment: T3
Genes in gene file:Dlx1as, Dlx6os2, Fam240b, Dlx6, Rassf2
Genes in zscore file:Olfr328, Hba-a2, Fam240b, Dlx6os2, Omt2a
Common genes found:Dlx6os2, Fam240b, NA, NA, NA
Cluster: C18 | Treatment: T3
Genes in gene file:Lingo1, Unc5d, Hs6st3, Vwa5b1, Satb2
Genes in zscore file:Krtap1-5, Krtap1-4, Krtap4-2, Retnlb, Dact2
Common genes found:Lingo1, Unc5d, Vwa5b1, Ubxn10, Dact2
Cluster: C19 | Treatment: T3
Genes in gene file:Tshz2, Itprid1, Adcy10, 4930556N09Rik, Chrm2
Genes in zscore file:Aard, Htr4, Htr7, B430306N03Rik, 4930546C10Rik
Common genes found:Tshz2, Itprid1, Adcy10, 4930556N09Rik, Htr2c
Cluster: C20 | Treatment: T3
Genes in gene file:4933412E24Rik, Fam19a1, 1700121N20Rik, Tarm1, Rab3c
Genes in zscore file:Dnah9, Acrv1, Olfr887, Olfr1395, Vmn1r-ps103
Common genes found:4933412E24Rik, Fam19a1, 1700121N20Rik, Tarm1, Rab3c
Cluster: C21 | Treatment: T3
Genes in gene file:Hs3st4, Foxp2, Etl4, Gm8630, Igsf21
Genes in zscore file:Igfbp4, Lrrc15, Crym, Daw1, Gal3st2
Common genes found:Hs3st4, Foxp2, Gm8630, Zfpm2, Grik3
Cluster: C22 | Treatment: T3
Genes in gene file:Nfix, Nxph1, Sox6, Kcnip1, Rbms3
Genes in zscore file:Olfr282, Klhl14, Gm572, Bend4, Sox6
Common genes found:Nxph1, Sox6, Kcnip1, Rbms3, Btbd11
Cluster: C23 | Treatment: T3
Genes in gene file:Adarb2, Erbb4, Dlx6os1, Zfp536, Gad1
Genes in zscore file:4930444F02Rik, Taar7f, Olfr1384, Olfr1383, Spp2
Common genes found:Adarb2, Erbb4, Dlx6os1, Zfp536, Gad1
Cluster: C24 | Treatment: T3
Genes in gene file:Prdm16, F630040K05Rik, Fgfr3, 0610039H22Rik, Sparcl1
Genes in zscore file:Gfap, 4930539C22Rik, Vmn1r64, Vmn1r78, Usp17le
Common genes found:Prdm16, F630040K05Rik, Fgfr3, 0610039H22Rik, Gm13872
Cluster: C25 | Treatment: T3
Genes in gene file:Olig2, Gpr17, Cspg4, Sox6, Kank1
Genes in zscore file:Ttr, Rpl13, Snord68, Inava, C1ql1
Common genes found:Gpr17, Cspg4, Cacng4, Lims2, Matn4
Cluster: C1 | Treatment: T4
Genes in gene file:Kcnh1, Gm35978, Shank2, A830018L16Rik, AU022754
Genes in zscore file:Fam107a, Celsr1, 5930403L14Rik, 4930544G11Rik, Vmn1r174
Common genes found:Prdm16, 4930447J18Rik, Fgfr2, Gm10863, NA
Cluster: C2 | Treatment: T4
Genes in gene file:Gria1, Tns3, Pitpnm2, Gab1, Olig1
Genes in zscore file:Hfe, Dennd1c, Lims2, Mmp2, Prps1l1
Common genes found:Tns3, Olig1, NA, NA, NA
Cluster: C3 | Treatment: T4
Genes in gene file:Kalrn, Gm10863, Rnf220, A230009B12Rik, Olig1
Genes in zscore file:2810429I04Rik, Acod1, Mucl1, S100a16, Tshb
Common genes found:Gm10863, Rnf220, A230009B12Rik, Olig1, Fgfr2
Cluster: C4 | Treatment: T4
Genes in gene file:Gm10863, Ppp1r14a, Cldn11, Trf, Nkx2-2
Genes in zscore file:Olfr1040, Fa2h, Fcrl6, Olfr1214, Gm7978
Common genes found:Ppp1r14a, Cldn11, Trf, Gjc3, Gng11
Cluster: C5 | Treatment: T4
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Lbh, Snai1, Epas1, 4933440J02Rik, Foxn3
No common genes found for treatment T4 in cluster C5
Cluster: C6 | Treatment: T4
Genes in gene file:Olig2, Olig1, Cacng4, Cspg4, Lims2
Genes in zscore file:Olfr324, Tmem100, Gm6961, Unc13c, Plac1
Common genes found:Olig2, Olig1, Cacng4, Cspg4, Lims2
Cluster: C7 | Treatment: T4
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Tcf7l1, Mir26b, Abca8b, Bmp6, Ppara
No common genes found for treatment T4 in cluster C7
Cluster: C8 | Treatment: T4
Genes in gene file:Prdm16, Msi2, Slc1a2, Sox1ot, Mertk
Genes in zscore file:Gpr179, Ranbp3l, Prodh, Frmpd1os, Lfng
Common genes found:Prdm16, Msi2, Slc1a2, Sox1ot, Mertk
Cluster: C9 | Treatment: T4
Genes in gene file:NA, NA, NA, NA, NA
Genes in zscore file:Gm15441, Loxl3, Dok1, Mir7040, Gm21057
No common genes found for treatment T4 in cluster C9
Cluster: C10 | Treatment: T4
Genes in gene file:Csf1r, Runx1, Cdh23, Irf8, Mertk
Genes in zscore file:Olfr433, Lgals9, Slfn3, Olfr1242, B020014A21Rik
Common genes found:Csf1r, Irf8, Fli1, Ccr6, Vsir
Cluster: C11 | Treatment: T4
Genes in gene file:Cx3cr1, Cdh23, Kalrn, Itgb5, E230029C05Rik
Genes in zscore file:Olfr784, Susd3, Ccr5, Dusp27, Havcr2
Common genes found:Cx3cr1, Cdh23, E230029C05Rik, Vsir, Csf1r
Cluster: C12 | Treatment: T4
Genes in gene file:Pdgfrb, Ebf1, Tns1, Zic4, Pitpnc1
Genes in zscore file:Adap2, Slc38a11, Olfr1000, Mylk2, Rarres2
Common genes found:Pdgfrb, Ebf1, Zic4, Pald1, Atp13a5
Cluster: C13 | Treatment: T4
Genes in gene file:Zic1, Zic4, 2610035F20Rik, Zic2, Colec12
Genes in zscore file:Iigp1, Pcolce, Mrgprb2, Jaml, Vmn2r114
Common genes found:Zic1, Zic4, 2610035F20Rik, Zic2, Colec12
Cluster: C14 | Treatment: T4
Genes in gene file:2310005E17Rik, Atp2b2, Zic2, Zic4, Rbfox3
Genes in zscore file:Olfr1490, Tas2r107, Olfr508, Samt4, Olfr384
Common genes found:2310005E17Rik, Zic4, Uaca, Cgnl1, Zic1
Cluster: C15 | Treatment: T4
Genes in gene file:Csf1r, Vsir, C1qc, Adgre1, Cd33
Genes in zscore file:Ccl6, Tnfrsf17, Clec4a3, Sash3, Batf3
Common genes found:Csf1r, Vsir, C1qc, Adgre1, Cd33
Cluster: C16 | Treatment: T4
Genes in gene file:Tshz2, AY702102, Cdh18, Man1c1, Spon1
Genes in zscore file:Pla2g2d, Vsig1, Gm22650, Trim58, A630019I02Rik
Common genes found:Tshz2, AY702102, Spon1, Hs3st4, Myzap
Cluster: C17 | Treatment: T4
Genes in gene file:Gad2, Dlx6os2, Daam2, Cntnap3, Mir486
Genes in zscore file:Olfr328, Hba-a2, Fam240b, Dlx6os2, Omt2a
Common genes found:Gad2, Dlx6os2, Cntnap3, Crhbp, NA
Cluster: C18 | Treatment: T4
Genes in gene file:Lingo1, Unc5d, D430041D05Rik, Dact2, Vwa5b1
Genes in zscore file:Krtap1-5, Krtap1-4, Krtap4-2, Retnlb, Dact2
Common genes found:Lingo1, Unc5d, Dact2, Vwa5b1, Ubxn10
Cluster: C19 | Treatment: T4
Genes in gene file:Tshz2, Gm8267, Bmper, Neurod6, Gpr161
Genes in zscore file:Aard, Htr4, Htr7, B430306N03Rik, 4930546C10Rik
Common genes found:Tshz2, Gm8267, Neurod6, Gpr161, Itprid1
Cluster: C20 | Treatment: T4
Genes in gene file:4933412E24Rik, Fam19a1, 1700121N20Rik, Fras1, Vat1l
Genes in zscore file:Dnah9, Acrv1, Olfr887, Olfr1395, Vmn1r-ps103
Common genes found:4933412E24Rik, Fam19a1, 1700121N20Rik, Fras1, Vat1l
Cluster: C21 | Treatment: T4
Genes in gene file:Hs3st4, Foxp2, Gm8630, Syt6, Zfpm2
Genes in zscore file:Igfbp4, Lrrc15, Crym, Daw1, Gal3st2
Common genes found:Hs3st4, Foxp2, Gm8630, Syt6, Zfpm2
Cluster: C22 | Treatment: T4
Genes in gene file:Nxph1, Btbd11, Grip1, Sox6, Dlgap2
Genes in zscore file:Olfr282, Klhl14, Gm572, Bend4, Sox6
Common genes found:Nxph1, Btbd11, Grip1, Sox6, Kcnmb2
Cluster: C23 | Treatment: T4
Genes in gene file:Adarb2, Erbb4, Zfp536, Dlx6os1, Grip1
Genes in zscore file:4930444F02Rik, Taar7f, Olfr1384, Olfr1383, Spp2
Common genes found:Adarb2, Erbb4, Zfp536, Dlx6os1, Grip1
Cluster: C24 | Treatment: T4
Genes in gene file:Prdm16, Mertk, Sox1ot, 9430041J12Rik, Fgfr3
Genes in zscore file:Gfap, 4930539C22Rik, Vmn1r64, Vmn1r78, Usp17le
Common genes found:Prdm16, Fgfr3, 5930403L14Rik, Gm13872, 2900052N01Rik
Cluster: C25 | Treatment: T4
Genes in gene file:Gpr17, Pdgfra, Lims2, Ntn1, Cacng4
Genes in zscore file:Ttr, Rpl13, Snord68, Inava, C1ql1
Common genes found:Gpr17, Pdgfra, Lims2, Ntn1, Cacng4
