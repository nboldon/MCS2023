library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

setwd("/Volumes/DataBox/GO_Analysis")

####################################
####################################
####################################

C3_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

#C3_T1vsT3[C3_T1vsT3$Log2FC >=0.5,][,"name"]

C3_T1vsT3[abs(C3_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T1vsT3[abs(C3_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

#BP = Biological Process
#MP = Molecular function
#CC = Cellular component

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C3_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T2vsT4[abs(C3_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T2vsT4[abs(C3_T2vsT4$Log2FC) >=0.5,][,"name"]
genes_to_test

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned 

####################################
####################################
####################################

C3_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C3_T124vsT3[abs(C3_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T124vsT3[abs(C3_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C4_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T1vsT3[abs(C4_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T1vsT3[abs(C4_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C4_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T2vsT4[abs(C4_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T2vsT4[abs(C4_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C6_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T1vsT3[abs(C6_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T1vsT3[abs(C6_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C6_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T2vsT4[abs(C6_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T2vsT4[abs(C6_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C6_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C6_T124vsT3[abs(C6_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T124vsT3[abs(C6_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C8_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T1vsT3[abs(C8_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T1vsT3[abs(C8_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C8_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T2vsT4[abs(C8_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T2vsT4[abs(C8_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C8_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C8_T124vsT3[abs(C8_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T124vsT3[abs(C8_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C10_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T1vsT3[abs(C10_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T1vsT3[abs(C10_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C10_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T2vsT4[abs(C10_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T2vsT4[abs(C10_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C10_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C10_T3vsT4[abs(C10_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T3vsT4[abs(C10_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T3vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C11_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T1vsT3[abs(C11_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T1vsT3[abs(C11_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C11_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T2vsT4[abs(C11_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T2vsT4[abs(C11_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C11_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C11_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C11_T124vsT3[abs(C11_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T124vsT3[abs(C11_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C11_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C12_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T1vsT3[abs(C12_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T1vsT3[abs(C12_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C12_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C12_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T2vsT4[abs(C12_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T2vsT4[abs(C12_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C12_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C13_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T1vsT3[abs(C13_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T1vsT3[abs(C13_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C13_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T2vsT4[abs(C13_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T2vsT4[abs(C13_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C14_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T1vsT3[abs(C14_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T1vsT3[abs(C14_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C14_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T2vsT4[abs(C14_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T2vsT4[abs(C14_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C15_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C15_T1vsT3[abs(C15_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T1vsT3[abs(C15_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C15_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C15_T2vsT4[abs(C15_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T2vsT4[abs(C15_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C15_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C15_T3vsT4[abs(C15_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T3vsT4[abs(C15_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T3vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C15_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C15_T124vsT3[abs(C15_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T124vsT3[abs(C15_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C16_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C16_T1vsT3[abs(C16_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T1vsT3[abs(C16_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C16_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C16_T2vsT4[abs(C16_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T2vsT4[abs(C16_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C16_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C16_T3vsT4[abs(C16_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T3vsT4[abs(C16_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T3vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C16_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C16_T124vsT3[abs(C16_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T124vsT3[abs(C16_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C17_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C17_T1vsT3[abs(C17_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T1vsT3[abs(C17_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C17_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C17_T2vsT4[abs(C17_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T2vsT4[abs(C17_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C17_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T3vsT4[abs(C17_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T3vsT4[abs(C17_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
#No GO results returned

####################################
####################################
####################################

C17_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C17_T124vsT3[abs(C17_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T124vsT3[abs(C17_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T1vsT3[abs(C18_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T1vsT3[abs(C18_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C18_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T2vsT4[abs(C18_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T2vsT4[abs(C18_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C18_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C18_T124vsT3[abs(C18_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T124vsT3[abs(C18_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C19_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C19_T1vsT3[abs(C19_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T1vsT3[abs(C19_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C19_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C19_T2vsT4[abs(C19_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T2vsT4[abs(C19_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C19_T1vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C19_T1vsT2[abs(C19_T1vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T1vsT2[abs(C19_T1vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T1vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C19_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C19_T124vsT3[abs(C19_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T124vsT3[abs(C19_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C20_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T1vsT3[abs(C20_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T1vsT3[abs(C20_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C20_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T2vsT4[abs(C20_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T2vsT4[abs(C20_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C20_T1vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T1vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T1vsT2[abs(C20_T1vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T1vsT2[abs(C20_T1vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C20_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T3vsT4[abs(C20_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T3vsT4[abs(C20_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T3vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C20_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C20_T124vsT3[abs(C20_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T124vsT3[abs(C20_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C21_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T1vsT3[abs(C21_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T1vsT3[abs(C21_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C21_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T2vsT4[abs(C21_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T2vsT4[abs(C21_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C21_T3vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T3vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T3vsT4[abs(C21_T3vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T3vsT4[abs(C21_T3vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T3vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C21_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C21_T124vsT3[abs(C21_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T124vsT3[abs(C21_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C22_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T1vsT3[abs(C22_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T1vsT3[abs(C22_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C22_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T2vsT4[abs(C22_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T2vsT4[abs(C22_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C22_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C22_T124vsT3[abs(C22_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T124vsT3[abs(C22_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C23_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T1vsT3[abs(C23_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T1vsT3[abs(C23_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C23_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C23_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T2vsT4[abs(C23_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T2vsT4[abs(C23_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No GO results returned

####################################
####################################
####################################

C23_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C23_T124vsT3[abs(C23_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T124vsT3[abs(C23_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C23_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C24_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T1vsT3[abs(C24_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T1vsT3[abs(C24_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C24_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T2vsT4[abs(C24_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T2vsT4[abs(C24_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C25_T1vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T1vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T1vsT3[abs(C25_T1vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T1vsT3[abs(C25_T1vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C25_T1vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C25_T2vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T2vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T2vsT4[abs(C25_T2vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T2vsT4[abs(C25_T2vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C25_T2vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C25_T124vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T1-2-4vsT3_FDR-0-1_Log2FC-0-5_2024-07-31.csv")

C25_T124vsT3[abs(C25_T124vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T124vsT3[abs(C25_T124vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C25_T124v3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################
####################################
####################################
####################################

C3_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C3_T3vsT124[abs(C3_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T3vsT124[abs(C3_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C11_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C11_T3vsT124[abs(C11_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T3vsT124[abs(C11_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C11_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C17_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C17_T3vsT124[abs(C17_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T3vsT124[abs(C17_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)
# No GO results returned

####################################
####################################
####################################

C18_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C18_T3vsT124[abs(C18_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T3vsT124[abs(C18_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C19_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C19_T3vsT124[abs(C19_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T3vsT124[abs(C19_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C20_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C20_T3vsT124[abs(C20_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T3vsT124[abs(C20_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)
# No GO results returned

####################################
####################################
####################################

C21_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C21_T3vsT124[abs(C21_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T3vsT124[abs(C21_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)
# No GO results returned

####################################
####################################
####################################

C22_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C22_T3vsT124[abs(C22_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T3vsT124[abs(C22_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

C23_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C23_T3vsT124[abs(C23_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T3vsT124[abs(C23_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C23_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################

