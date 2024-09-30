library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

setwd("/Volumes/DataBox/GO_Analysis")

####################################
####################################
####################################

C3_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T1vsT4[abs(C3_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T1vsT4[abs(C3_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

#BP = Biological Process
#MP = Molecular function
#CC = Cellular component

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C3_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C3_T2vsT1[abs(C3_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T2vsT1[abs(C3_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C3_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T2vsT3[abs(C3_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T2vsT3[abs(C3_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C3_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C3_T3vsT1[abs(C3_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T3vsT1[abs(C3_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C3_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C3_T3vsT2[abs(C3_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T3vsT2[abs(C3_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C3_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T4vsT1[abs(C3_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T4vsT1[abs(C3_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C3_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C3_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T4vsT2[abs(C3_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T4vsT2[abs(C3_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C3_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C3_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C3_T4vsT3[abs(C3_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C3_T4vsT3[abs(C3_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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
####################################
####################################
####################################


C4_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T1vsT4[abs(C4_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T1vsT4[abs(C4_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C4_T2vsT1[abs(C4_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T2vsT1[abs(C4_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T2vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T2vsT3[abs(C4_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T2vsT3[abs(C4_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C4_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C4_T3vsT1[abs(C4_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T3vsT1[abs(C4_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C4_T3vsT2[abs(C4_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T3vsT2[abs(C4_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T4vsT1[abs(C4_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T4vsT1[abs(C4_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T4vsT2[abs(C4_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T4vsT2[abs(C4_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C4_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C4_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C4_T4vsT3[abs(C4_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T4vsT3[abs(C4_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C4_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C4_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C4_T3vsT124[abs(C4_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C4_T3vsT124[abs(C4_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C6_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T1vsT4[abs(C6_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T1vsT4[abs(C6_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C6_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C6_T2vsT1[abs(C6_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T2vsT1[abs(C6_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# no return

####################################
####################################
####################################

C6_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T2vsT3[abs(C6_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T2vsT3[abs(C6_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C6_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C6_T3vsT1[abs(C6_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T3vsT1[abs(C6_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C6_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C6_T3vsT2[abs(C6_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T3vsT2[abs(C6_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C6_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T4vsT1[abs(C6_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T4vsT1[abs(C6_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C6_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T4vsT2[abs(C6_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T4vsT2[abs(C6_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C6_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C6_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C6_T4vsT3[abs(C6_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T4vsT3[abs(C6_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C6_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C6_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C6_T3vsT124[abs(C6_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C6_T3vsT124[abs(C6_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C7_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C7_T1vsT4[abs(C7_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T1vsT4[abs(C7_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C7_T2vsT1[abs(C7_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T2vsT1[abs(C7_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C7_T2vsT3[abs(C7_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T2vsT3[abs(C7_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C7_T3vsT1[abs(C7_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T3vsT1[abs(C7_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C7_T3vsT2[abs(C7_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T3vsT2[abs(C7_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C7_T4vsT1[abs(C7_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T4vsT1[abs(C7_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C7_T4vsT2[abs(C7_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T4vsT2[abs(C7_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C7_T4vsT3[abs(C7_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T4vsT3[abs(C7_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C7_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C7_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C7_T3vsT124[abs(C7_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C7_T3vsT124[abs(C7_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C8_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T1vsT4[abs(C8_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T1vsT4[abs(C8_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C8_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C8_T2vsT1[abs(C8_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T2vsT1[abs(C8_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C8_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T2vsT3[abs(C8_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T2vsT3[abs(C8_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C8_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C8_T3vsT1[abs(C8_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T3vsT1[abs(C8_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C8_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C8_T3vsT2[abs(C8_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T3vsT2[abs(C8_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C8_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T4vsT1[abs(C8_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T4vsT1[abs(C8_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C8_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T4vsT2[abs(C8_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T4vsT2[abs(C8_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C8_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C8_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C8_T4vsT3[abs(C8_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T4vsT3[abs(C8_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C8_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C8_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C8_T3vsT124[abs(C8_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C8_T3vsT124[abs(C8_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C9_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C9_T1vsT4[abs(C9_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T1vsT4[abs(C9_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C9_T2vsT1[abs(C9_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T2vsT1[abs(C9_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C9_T2vsT3[abs(C9_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T2vsT3[abs(C9_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C9_T3vsT1[abs(C9_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T3vsT1[abs(C9_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C9_T3vsT2[abs(C9_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T3vsT2[abs(C9_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C9_T4vsT1[abs(C9_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T4vsT1[abs(C9_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C9_T4vsT2[abs(C9_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T4vsT2[abs(C9_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C9_T4vsT3[abs(C9_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T4vsT3[abs(C9_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C9_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C9_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C9_T3vsT124[abs(C9_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C9_T3vsT124[abs(C9_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C10_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T1vsT4[abs(C10_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T1vsT4[abs(C10_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C10_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C10_T2vsT1[abs(C10_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T2vsT1[abs(C10_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T2vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T2vsT3[abs(C10_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T2vsT3[abs(C10_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C10_T3vsT1[abs(C10_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T3vsT1[abs(C10_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C10_T3vsT2[abs(C10_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T3vsT2[abs(C10_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T4vsT1[abs(C10_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T4vsT1[abs(C10_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T4vsT2[abs(C10_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T4vsT2[abs(C10_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C10_T4vsT3[abs(C10_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T4vsT3[abs(C10_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C10_T4vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C10_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C10_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C10_T3vsT124[abs(C10_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C10_T3vsT124[abs(C10_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C11_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T1vsT4[abs(C11_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T1vsT4[abs(C11_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C11_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C11_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C11_T2vsT1[abs(C11_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T2vsT1[abs(C11_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C11_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T2vsT3[abs(C11_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T2vsT3[abs(C11_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C11_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C11_T3vsT1[abs(C11_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T3vsT1[abs(C11_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C11_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C11_T3vsT2[abs(C11_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T3vsT2[abs(C11_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C11_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T4vsT1[abs(C11_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T4vsT1[abs(C11_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C11_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C11_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T4vsT2[abs(C11_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T4vsT2[abs(C11_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C11_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C11_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C11_T4vsT3[abs(C11_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C11_T4vsT3[abs(C11_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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
####################################
####################################
####################################


C12_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T1vsT4[abs(C12_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T1vsT4[abs(C12_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C12_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C12_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C12_T2vsT1[abs(C12_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T2vsT1[abs(C12_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C12_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T2vsT3[abs(C12_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T2vsT3[abs(C12_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C12_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C12_T3vsT1[abs(C12_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T3vsT1[abs(C12_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C12_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C12_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C12_T3vsT2[abs(C12_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T3vsT2[abs(C12_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C12_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T4vsT1[abs(C12_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T4vsT1[abs(C12_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C12_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C12_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T4vsT2[abs(C12_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T4vsT2[abs(C12_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C12_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C12_T4vsT3[abs(C12_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T4vsT3[abs(C12_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C12_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C12_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C12_T3vsT124[abs(C12_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C12_T3vsT124[abs(C12_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C13_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T1vsT4[abs(C13_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T1vsT4[abs(C13_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C13_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C13_T2vsT1[abs(C13_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T2vsT1[abs(C13_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C13_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T2vsT3[abs(C13_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T2vsT3[abs(C13_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C13_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C13_T3vsT1[abs(C13_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T3vsT1[abs(C13_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C13_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C13_T3vsT2[abs(C13_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T3vsT2[abs(C13_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C13_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T4vsT1[abs(C13_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T4vsT1[abs(C13_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C13_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T4vsT2[abs(C13_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T4vsT2[abs(C13_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C13_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C13_T4vsT3[abs(C13_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T4vsT3[abs(C13_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C13_T4vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C13_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C13_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C13_T3vsT124[abs(C13_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C13_T3vsT124[abs(C13_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C14_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T1vsT4[abs(C14_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T1vsT4[abs(C14_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C14_T2vsT1[abs(C14_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T2vsT1[abs(C14_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C14_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T2vsT3[abs(C14_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T2vsT3[abs(C14_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C14_T3vsT1[abs(C14_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T3vsT1[abs(C14_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C14_T3vsT2[abs(C14_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T3vsT2[abs(C14_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T4vsT1[abs(C14_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T4vsT1[abs(C14_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T4vsT2[abs(C14_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T4vsT2[abs(C14_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C14_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C14_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C14_T4vsT3[abs(C14_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T4vsT3[abs(C14_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C14_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C14_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C14_T3vsT124[abs(C14_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C14_T3vsT124[abs(C14_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################

C15_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C15_T1vsT4[abs(C15_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T1vsT4[abs(C15_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C15_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C15_T2vsT1[abs(C15_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T2vsT1[abs(C15_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C15_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C15_T2vsT3[abs(C15_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T2vsT3[abs(C15_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C15_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C15_T3vsT1[abs(C15_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T3vsT1[abs(C15_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C15_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C15_T3vsT2[abs(C15_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T3vsT2[abs(C15_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C15_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C15_T4vsT1[abs(C15_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T4vsT1[abs(C15_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C15_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C15_T4vsT2[abs(C15_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T4vsT2[abs(C15_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C15_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C15_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C15_T4vsT3[abs(C15_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T4vsT3[abs(C15_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C15_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C15_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C15_T3vsT124[abs(C15_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C15_T3vsT124[abs(C15_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C16_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C16_T1vsT4[abs(C16_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T1vsT4[abs(C16_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C16_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C16_T2vsT1[abs(C16_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T2vsT1[abs(C16_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C16_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C16_T2vsT3[abs(C16_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T2vsT3[abs(C16_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C16_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C16_T3vsT1[abs(C16_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T3vsT1[abs(C16_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C16_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C16_T3vsT2[abs(C16_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T3vsT2[abs(C16_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C16_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C16_T4vsT1[abs(C16_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T4vsT1[abs(C16_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C16_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C16_T4vsT2[abs(C16_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T4vsT2[abs(C16_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C16_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C16_T4vsT3[abs(C16_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T4vsT3[abs(C16_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C16_T4vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C16_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C16_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C16_T3vsT124[abs(C16_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C16_T3vsT124[abs(C16_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C17_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C17_T1vsT4[abs(C17_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T1vsT4[abs(C17_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C17_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T2vsT1[abs(C17_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T2vsT1[abs(C17_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C17_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C17_T2vsT3[abs(C17_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T2vsT3[abs(C17_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C17_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C17_T3vsT1[abs(C17_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T3vsT1[abs(C17_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C17_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T3vsT2[abs(C17_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T3vsT2[abs(C17_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C17_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T4vsT1[abs(C17_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T4vsT1[abs(C17_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C17_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T4vsT2[abs(C17_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T4vsT2[abs(C17_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C17_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C17_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C17_T4vsT3[abs(C17_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T4vsT3[abs(C17_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C17_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C17_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C17_T3vsT124[abs(C17_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C17_T3vsT124[abs(C17_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C18_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T1vsT4[abs(C18_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T1vsT4[abs(C18_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C18_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T2vsT1[abs(C18_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T2vsT1[abs(C18_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C18_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T2vsT3[abs(C18_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T2vsT3[abs(C18_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C18_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C18_T3vsT1[abs(C18_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T3vsT1[abs(C18_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C18_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C18_T3vsT2[abs(C18_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T3vsT2[abs(C18_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C18_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T4vsT1[abs(C18_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T4vsT1[abs(C18_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C18_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T4vsT2[abs(C18_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T4vsT2[abs(C18_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C18_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C18_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C18_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C18_T4vsT3[abs(C18_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C18_T4vsT3[abs(C18_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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
####################################
####################################
####################################

C19_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C19_T1vsT4[abs(C19_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T1vsT4[abs(C19_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C19_T2vsT1[abs(C19_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T2vsT1[abs(C19_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-02.csv")

C19_T2vsT3[abs(C19_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T2vsT3[abs(C19_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-06-26.csv")

C19_T3vsT1[abs(C19_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T3vsT1[abs(C19_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C19_T3vsT2[abs(C19_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T3vsT2[abs(C19_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T3vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C19_T4vsT1[abs(C19_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T4vsT1[abs(C19_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C19_T4vsT2[abs(C19_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T4vsT2[abs(C19_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C19_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C19_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C19_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C19_T4vsT3[abs(C19_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C19_T4vsT3[abs(C19_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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
####################################
####################################
####################################


C20_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T1vsT4[abs(C20_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T1vsT4[abs(C20_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C20_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T2vsT1[abs(C20_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T2vsT1[abs(C20_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C20_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T2vsT3[abs(C20_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T2vsT3[abs(C20_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C20_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T3vsT1[abs(C20_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T3vsT1[abs(C20_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C20_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T3vsT2[abs(C20_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T3vsT2[abs(C20_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C20_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T4vsT1[abs(C20_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T4vsT1[abs(C20_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C20_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T4vsT2[abs(C20_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T4vsT2[abs(C20_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C20_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C20_T4vsT3[abs(C20_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T4vsT3[abs(C20_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C20_T4vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C20_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C20_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C20_T3vsT124[abs(C20_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C20_T3vsT124[abs(C20_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C21_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T1vsT4[abs(C21_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T1vsT4[abs(C21_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C21_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T2vsT1[abs(C21_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T2vsT1[abs(C21_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T2vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C21_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T2vsT3[abs(C21_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T2vsT3[abs(C21_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C21_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T3vsT1[abs(C21_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T3vsT1[abs(C21_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C21_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T3vsT2[abs(C21_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T3vsT2[abs(C21_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C21_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T4vsT1[abs(C21_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T4vsT1[abs(C21_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C21_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C21_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T4vsT2[abs(C21_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T4vsT2[abs(C21_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C21_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C21_T4vsT3[abs(C21_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T4vsT3[abs(C21_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C21_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C21_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C21_T3vsT124[abs(C21_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C21_T3vsT124[abs(C21_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C22_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T1vsT4[abs(C22_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T1vsT4[abs(C22_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C22_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T2vsT1[abs(C22_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T2vsT1[abs(C22_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C22_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T2vsT3[abs(C22_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T2vsT3[abs(C22_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C22_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T3vsT1[abs(C22_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T3vsT1[abs(C22_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C22_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T3vsT2[abs(C22_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T3vsT2[abs(C22_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C22_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T4vsT1[abs(C22_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T4vsT1[abs(C22_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C22_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C22_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T4vsT2[abs(C22_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T4vsT2[abs(C22_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C22_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C22_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C22_T4vsT3[abs(C22_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C22_T4vsT3[abs(C22_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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

png("C3_T3vsT124_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

####################################
####################################
####################################
####################################
####################################
####################################


C23_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T1vsT4[abs(C23_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T1vsT4[abs(C23_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T2vsT1[abs(C23_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T2vsT1[abs(C23_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T2vsT3[abs(C23_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T2vsT3[abs(C23_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T3vsT1[abs(C23_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T3vsT1[abs(C23_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C23_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C23_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T3vsT2[abs(C23_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T3vsT2[abs(C23_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T4vsT1[abs(C23_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T4vsT1[abs(C23_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T4vsT2[abs(C23_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T4vsT2[abs(C23_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C23_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C23_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C23_T4vsT3[abs(C23_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C23_T4vsT3[abs(C23_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

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
####################################
####################################
####################################


C24_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T1vsT4[abs(C24_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T1vsT4[abs(C24_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T1vsT4_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C24_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T2vsT1[abs(C24_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T2vsT1[abs(C24_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C24_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T2vsT3[abs(C24_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T2vsT3[abs(C24_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C24_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T3vsT1[abs(C24_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T3vsT1[abs(C24_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C24_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T3vsT2[abs(C24_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T3vsT2[abs(C24_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C24_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T4vsT1[abs(C24_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T4vsT1[abs(C24_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T4vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C24_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T4vsT2[abs(C24_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T4vsT2[abs(C24_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C24_T4vsT2_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C24_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C24_T4vsT3[abs(C24_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T4vsT3[abs(C24_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C24_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C24_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C24_T3vsT124[abs(C24_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C24_T3vsT124[abs(C24_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################


C25_T1vsT4 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T1vsT4Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T1vsT4[abs(C25_T1vsT4$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T1vsT4[abs(C25_T1vsT4$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T2vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T2vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T2vsT1[abs(C25_T2vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T2vsT1[abs(C25_T2vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T2vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T2vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T2vsT3[abs(C25_T2vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T2vsT3[abs(C25_T2vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C25_T2vsT3_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C25_T3vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T3vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T3vsT1[abs(C25_T3vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T3vsT1[abs(C25_T3vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20))
#Displays top 20 

png("C25_T3vsT1_top20.png", res = 250, width = 1200, height = 2750)
print(fit)
dev.off()

fit

####################################
####################################
####################################

C25_T3vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T3vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T3vsT2[abs(C25_T3vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T3vsT2[abs(C25_T3vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T4vsT1 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T4vsT1Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T4vsT1[abs(C25_T4vsT1$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T4vsT1[abs(C25_T4vsT1$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T4vsT2 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T4vsT2Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T4vsT2[abs(C25_T4vsT2$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T4vsT2[abs(C25_T4vsT2$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T4vsT3 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T4vsT3Comp_FDR-0-1_Log2FC-0-5_2024-07-05.csv")

C25_T4vsT3[abs(C25_T4vsT3$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T4vsT3[abs(C25_T4vsT3$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################

C25_T3vsT124 <-read.csv("/Volumes/DataBox/ProjMCS6/ResearchQuestions/C25_T3vsT1-2-4_FDR-0-1_Log2FC-0-5_2024-08-09.csv")

C25_T3vsT124[abs(C25_T3vsT124$Log2FC) >=0.5,][,"name"]

genes_to_test <- C25_T3vsT124[abs(C25_T3vsT124$Log2FC) >=0.5,][,"name"]

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
# No return

####################################
####################################
####################################
####################################
####################################
####################################

