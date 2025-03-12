#Setup an interactive session
salloc --account=eon -t 0-16:00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Updated conda env 12-2023
module load miniconda3/23.1.0

conda activate archr2023_12

#Load libraries
R
library(ArchR)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(Seurat)
library(Signac)
library(BiocManager)
library(BiocGenerics)

#Additional setup
setwd("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7")
addArchRGenome("mm10")
addArchRThreads(threads = 16)

#Load project
projMCS5 <- loadArchRProject(path = "/project/eon/nboldon/MCS2023/Save-ProjMCS5", force = FALSE, showLogo = FALSE)

##########################################  

# t1 = 2N
# t2 = 2N+
# t3 = Ts
# t4 = Ts+

############################################
############################################

### This code does not appear to work to subset specific samples

# Subset the project by desired sample(s) of interest
C25ArchRSubset <- projMCS5[projMCS5$Clusters==("C25"),]

C25ArchRSubset

# Get peak marker features
C25ArchRSubset <- getMarkerFeatures(
	ArchRProj = C25ArchRSubset, 
	useMatrix = "PeakMatrix",
	groupBy = "Sample",
	bias = c("TSSEnrichment","log10(nFrags)"),
	testMethod = "wilcoxon"
)

C24ArchRSubset 

C24_PeakMarker_List <- getMarkers(
	C24ArchRSubset,
	cutOff = "FDR<=0.01 & abs(Log2FC)>=1.25")
	
C24_PeakMarker_List

print(C24_PeakMarker_List)

write.csv(C24_PeakMarker_List, file="/project/eon/nboldon/MCS2023/fragAnalysis_TSS7/C24_PeakMarker_List.csv", row.names=FALSE)


#########################

### Returns for getMarkerFeatures

C1ArchRSubset
#2024-01-23 20:17:51.164971 : Matching Known Biases, 0.009 mins elapsed.
#2024-01-23 20:17:51.233476 : Found less than 100 cells for background matching, Lowering k to 15
#2024-01-23 20:17:51.251183 : Found less than 100 cells for background matching, Lowering k to 3
#2024-01-23 20:17:51.267151 : Found less than 100 cells for background matching, Lowering k to 14
#2024-01-23 20:17:51.284092 : Found less than 100 cells for background matching, Lowering k to 20
#2024-01-23 20:17:51.30279 : Found less than 100 cells for background matching, Lowering k to 22
#2024-01-23 20:17:51.32009 : Found less than 100 cells for background matching, Lowering k to 12
#2024-01-23 20:17:51.3366 : Found less than 100 cells for background matching, Lowering k to 3
#2024-01-23 20:17:51.352582 : Found less than 100 cells for background matching, Lowering k to 12
#2024-01-23 20:17:51.369076 : Found less than 100 cells for background matching, Lowering k to 13
#2024-01-23 20:17:51.385702 : Found less than 100 cells for background matching, Lowering k to 28
#2024-01-23 20:17:51.402452 : Found less than 100 cells for background matching, Lowering k to 20
#2024-01-23 20:17:51.419199 : Found less than 100 cells for background matching, Lowering k to 51
#2024-01-23 20:17:51.437197 : Found less than 100 cells for background matching, Lowering k to 18
#2024-01-23 20:17:51.454166 : Found less than 100 cells for background matching, Lowering k to 13
#2024-01-23 20:17:51.472334 : Found less than 100 cells for background matching, Lowering k to 0
#Error: Cloud has no points

C2ArchRSubset
#2024-01-23 20:20:49.526524 : Found less than 100 cells for background matching, Lowering k to 11
#2024-01-23 20:20:49.544098 : Found less than 100 cells for background matching, Lowering k to 0
#Error: Cloud has no points

C4ArchRSubset
2024-01-24 13:02:51.588151 : Computing Pairwise Tests (56 of 56), 2.371 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 1 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C301_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C301_>
Error Found Iteration 3 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C303_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C303_>
Error Found Iteration 5 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C305_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C305_>
Error Found Iteration 33 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C338_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C338_>
Error Found Iteration 46 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background 
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  6 function calls resulted in an error

C5ArchRSubset
2024-01-24 13:04:14.155222 : Found less than 100 cells for background matching, Lowering k to 0
Error: Cloud has no points

C6ArchRSubset
2024-01-24 13:09:45.028582 : Computing Pairwise Tests (56 of 56), 2.369 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 7 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C307_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C307_>
Error Found Iteration 10 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C310_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C310_>
Error Found Iteration 14 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C315_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C315_>
Error Found Iteration 17 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C319_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C319_>
Error Found Iteration 18 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and backgroun
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  19 function calls resulted in an error

C7ArchRSubset
2024-01-24 13:08:43.573214 : Matching Known Biases, 0.024 mins elapsed.
2024-01-24 13:08:43.635333 : Found less than 100 cells for background matching, Lowering k to 0
Error: Cloud has no points

C9ArchRSubset
2024-01-24 13:12:52.589709 : Found less than 100 cells for background matching, Lowering k to 3
2024-01-24 13:12:52.607722 : Found less than 100 cells for background matching, Lowering k to 0
Error: Cloud has no points

C10ArchRSubset
2024-01-24 13:16:32.102094 : Computing Pairwise Tests (56 of 56), 2.078 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 25 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C328_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C328_>
Error Found Iteration 27 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C332_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C332_>
Error Found Iteration 32 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C337_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C337_>
Error Found Iteration 52 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C360_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C360_>
Error Found Iteration 53 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and backgrou
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  5 function calls resulted in an error

C11ArchRSubset
2024-01-24 13:24:01.095443 : Found less than 100 cells for background matching, Lowering k to 0
Error: Cloud has no points

C14ArchRSubset
2024-01-24 13:36:33.848919 : Computing Pairwise Tests (56 of 56), 1.295 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 4 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C304_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C304_>
Error Found Iteration 5 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C305_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C305_>
Error Found Iteration 6 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C306_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C306_>
Error Found Iteration 9 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C309_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C309_>
Error Found Iteration 14 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background a
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  23 function calls resulted in an error

C17ArchRSubset
2024-01-24 13:45:13.814755 : Computing Pairwise Tests (56 of 56), 2.102 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 1 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C301_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C301_>
Error Found Iteration 12 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C313_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C313_>
Error Found Iteration 15 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C316_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C316_>
Error Found Iteration 17 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C319_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C319_>
Error Found Iteration 18 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and backgroun
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  16 function calls resulted in an error

C19ArchRSubset
2024-01-24 13:48:13.925806 : Computing Pairwise Tests (56 of 56), 1.893 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 7 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C307_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C307_>
Error Found Iteration 8 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C308_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C308_>
Error Found Iteration 14 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C315_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C315_>
Error Found Iteration 18 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C320_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C320_>
Error Found Iteration 19 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  8 function calls resulted in an error

C20ArchRSubset
024-01-24 13:50:49.732185 : Computing Pairwise Tests (56 of 56), 2.802 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 3 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C303_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C303_>
Error Found Iteration 4 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C304_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C304_>
Error Found Iteration 24 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C327_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C327_>
Error Found Iteration 32 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C337_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C337_>
Error Found Iteration 34 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  7 function calls resulted in an error

C23ArchRSubset
2024-01-24 14:03:53.363902 : Computing Pairwise Tests (56 of 56), 2.962 mins elapsed.
Error in .safelapply(seq_along(matchObj[[1]]), function(x) { : 
Error Found Iteration 23 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C325_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C325_>
Error Found Iteration 33 : 
	[1] "Error in FUN(X[[i]], ...) : \n  Cells in foreground and background are 0 for group = C338_\n"
	<simpleError in FUN(X[[i]], ...): Cells in foreground and background are 0 for group = C338_>
In addition: Warning message:
In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
  2 function calls resulted in an error

C25ArchRSubset
Error: Cloud has no points
