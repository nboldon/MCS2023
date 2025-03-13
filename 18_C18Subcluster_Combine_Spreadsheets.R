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

#######################

# Set working directory
setwd("/Volumes/DataBox/MCS2023/C18_Subset/")

# Read spreadsheets into data frames and add a source column
df1 = subset(read.csv("C1_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df2 = subset(read.csv("C2_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df3 = subset(read.csv("C3_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df4 = subset(read.csv("C4_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df5 = subset(read.csv("C5_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df6 = subset(read.csv("C6_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df7 = subset(read.csv("C7_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df8 = subset(read.csv("C8_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df9 = subset(read.csv("C9_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df10 = subset(read.csv("C10_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df11 = subset(read.csv("C11_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df12 = subset(read.csv("C12_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df13 = subset(read.csv("C13_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df14 = subset(read.csv("C14_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df15 = subset(read.csv("C15_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df16 = subset(read.csv("C16_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df17 = subset(read.csv("C17_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
df18 = subset(read.csv("C18_C18ArchRSubset_2024-06-19.csv"), select = -c(X, MeanDiff))
names(df1) = c("seqnames", "start", "end", "strand", "name", "idx", "C1.Log2FC", "C1.FDR")
names(df2) = c("seqnames", "start", "end", "strand", "name", "idx", "C2.Log2FC", "C2.FDR")
names(df3) = c("seqnames", "start", "end", "strand", "name", "idx", "C3.Log2FC", "C3.FDR")
names(df4) = c("seqnames", "start", "end", "strand", "name", "idx", "C4.Log2FC", "C4.FDR")
names(df5) = c("seqnames", "start", "end", "strand", "name", "idx", "C5.Log2FC", "C5.FDR")
names(df6) = c("seqnames", "start", "end", "strand", "name", "idx", "C6.Log2FC", "C6.FDR")
names(df7) = c("seqnames", "start", "end", "strand", "name", "idx", "C7.Log2FC", "C7.FDR")
names(df8) = c("seqnames", "start", "end", "strand", "name", "idx", "C8.Log2FC", "C8.FDR")
names(df9) = c("seqnames", "start", "end", "strand", "name", "idx", "C9.Log2FC", "C9.FDR")
names(df10) = c("seqnames", "start", "end", "strand", "name", "idx", "C10.Log2FC", "C10.FDR")
names(df11) = c("seqnames", "start", "end", "strand", "name", "idx", "C11.Log2FC", "C11.FDR")
names(df12) = c("seqnames", "start", "end", "strand", "name", "idx", "C12.Log2FC", "C12.FDR")
names(df13) = c("seqnames", "start", "end", "strand", "name", "idx", "C13.Log2FC", "C13.FDR")
names(df14) = c("seqnames", "start", "end", "strand", "name", "idx", "C14.Log2FC", "C14.FDR")
names(df15) = c("seqnames", "start", "end", "strand", "name", "idx", "C15.Log2FC", "C15.FDR")
names(df16) = c("seqnames", "start", "end", "strand", "name", "idx", "C16.Log2FC", "C16.FDR")
names(df17) = c("seqnames", "start", "end", "strand", "name", "idx", "C17.Log2FC", "C17.FDR")
names(df18) = c("seqnames", "start", "end", "strand", "name", "idx", "C18.Log2FC", "C18.FDR")

combined_df = Reduce(function(a, b) merge(a, b, all = TRUE), list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18))

# Export combined data frame to a new spreadsheet
write.csv(combined_df, "C18Subcluster_TxMarkers_Combined_2024-06-19.csv", row.names = FALSE)
