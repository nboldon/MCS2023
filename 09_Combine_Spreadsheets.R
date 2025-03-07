# Set working directory
setwd("/Volumes/DataBox/MCS2023/Tx_Comp/")

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("T1_C1_GeneMarkers_2024-04-10.csv"), select = -c(X, MeanDiff))
df2 = subset(read.csv("T2_C1_GeneMarkers_2024-04-10.csv"), select = -c(X, MeanDiff))
df3 = subset(read.csv("T3_C1_GeneMarkers_2024-04-10.csv"), select = -c(X, MeanDiff))
df4 = subset(read.csv("T4_C1_GeneMarkers_2024-04-10.csv"), select = -c(X, MeanDiff))
names(df1) = c("seqnames", "start", "end", "strand", "name", "idx", "T1.Log2FC", "T1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "idx", "T2.Log2FC", "T2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "idx", "T3.Log2FC", "T3.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "idx", "T4.Log2FC", "T4.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C1_TxMarkers_Combined.csv", row.names = FALSE)

