library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

setwd("/project/eon/nboldon/MCS2023/Subset/C18_Subset")

####################################
####################################
####################################

C18.3_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.3_T124vsT3[abs(C18.3_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.3_T124vsT3[abs(C18.3_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.3_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.4_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.4_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.4_T124vsT3[abs(C18.4_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.4_T124vsT3[abs(C18.4_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.4_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.5_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.5_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.5_T124vsT3[abs(C18.5_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.5_T124vsT3[abs(C18.5_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.5_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.6_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.6_T124vsT3[abs(C18.6_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.6_T124vsT3[abs(C18.6_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.6_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.7_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.7_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.7_T124vsT3[abs(C18.7_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.7_T124vsT3[abs(C18.7_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.7_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.8_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.8_T124vsT3[abs(C18.8_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.8_T124vsT3[abs(C18.8_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C18.10_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.10_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.10_T124vsT3[abs(C18.10_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.10_T124vsT3[abs(C18.10_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.10_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.11_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.11_T124vsT3[abs(C18.11_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.11_T124vsT3[abs(C18.11_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.11_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.12_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.12_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.12_T124vsT3[abs(C18.12_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.12_T124vsT3[abs(C18.12_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.12_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.13_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.13_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.13_T124vsT3[abs(C18.13_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.13_T124vsT3[abs(C18.13_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.13_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.14_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.14_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.14_T124vsT3[abs(C18.14_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.14_T124vsT3[abs(C18.14_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.14_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.15_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.15_T124vsT3[abs(C18.15_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.15_T124vsT3[abs(C18.15_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.15_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.16_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.16_T124vsT3[abs(C18.16_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.16_T124vsT3[abs(C18.16_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18.16_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18.17_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.17_T124vsT3[abs(C18.17_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.17_T124vsT3[abs(C18.17_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C18.18_T124vsT3 <-read.csv("/project/eon/nboldon/MCS2023/Subset/C18_Subset/C18.18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18.18_T124vsT3[abs(C18.18_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18.18_T124vsT3[abs(C18.18_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

