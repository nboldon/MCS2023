C3_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T1vsT4[abs(C3_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T1vsT4[abs(C3_T1vsT4$Log2FC) >=0.5,][,"name"]

##########################################
##########################################
##########################################

library(topGO)
library(ggplot2)

# Prepare the gene list (a named vector of scores, e.g., Log2FC)
geneList <- C3_T1vsT4$Log2FC
names(geneList) <- C3_T1vsT4$name

# Create a topGOdata object
GOdata <- new("topGOdata", 
              ontology = "BP",  # BP = Biological Process
              allGenes = geneList, 
              geneSel = function(x) abs(x) >= 0.5,  # Gene selection criterion
              nodeSize = 10,  # Minimum number of annotated genes per GO term
              annot = annFUN.org, 
              mapping = "org.Mm.eg.db", 
              ID = "SYMBOL")

# Run enrichment analysis with Fisher's exact test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Summarize results
GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)

# Get the top significant terms (e.g., top 10)
top_results <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)

# Bar plot using ggplot2
ggplot(top_results, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)))) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Top 10 GO Biological Process Terms", x = "GO Term", y = "-log10(Fisher p-value)") +
  theme_minimal()

# Dot plot
ggplot(top_results, aes(x = reorder(Term, -as.numeric(classicFisher)), y = -log10(as.numeric(classicFisher)), size = Significant)) +
  geom_point() +
  coord_flip() + 
  labs(title = "Top 10 GO Terms", x = "GO Term", y = "-log10(Fisher p-value)", size = "Gene Count") +
  theme_minimal()

##########################################
##########################################
##########################################

library(goseq)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Create a binary vector indicating which genes are differentially expressed
geneVector <- as.integer(abs(C3_T1vsT4$Log2FC) >= 0.5)
names(geneVector) <- C3_T1vsT4$name

# Perform GO analysis
pwf <- nullp(geneVector, "mm10", "geneSymbol")  # Adjust for length bias
GO_results <- goseq(pwf, "mm10", "geneSymbol")

# Filter significant results
significantGO <- GO_results[GO_results$over_represented_pvalue < 0.05, ]
head(significantGO)

# Select top 10 terms for visualization
top_results <- head(significantGO[order(significantGO$over_represented_pvalue), ], 10)

# Bar plot of top 10 GO terms
ggplot(top_results, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  labs(title = "Top 10 GO Biological Process Terms", x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()

# Bubble plot for GO terms
ggplot(top_results, aes(x = reorder(term, -log10(over_represented_pvalue)), y = -log10(over_represented_pvalue), size = numDEInCat)) +
  geom_point() +
  coord_flip() +
  labs(title = "Top 10 GO Terms", x = "GO Term", y = "-log10(p-value)", size = "Gene Count") +
  theme_minimal()


##########################################
##########################################
##########################################

library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve GO terms for the genes
go_data <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), 
                 filters = "external_gene_name", 
                 values = genes_to_test, 
                 mart = ensembl)
head(go_data)

# Count the number of genes associated with each GO term
go_counts <- as.data.frame(table(go_data$name_1006))

# Plot the top 10 GO terms by gene count
top_go <- head(go_counts[order(go_counts$Freq, decreasing = TRUE), ], 10)

ggplot(top_go, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 GO Terms by Gene Count", x = "GO Term", y = "Gene Count") +
  theme_minimal()

##########################################
##########################################
##########################################

# Custom bubble plot
ggplot(data = top_results, aes(x = reorder(Term, -log10(classicFisher)), y = -log10(classicFisher), size = Significant)) +
  geom_point(aes(color = -log10(classicFisher))) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "GO Term Bubble Plot", x = "GO Term", y = "-log10(p-value)", size = "Gene Count", color = "-log10(p-value)") +
  theme_minimal()
