## 03_projMCS1.R
- Adds doublet scores
- Creates dataframes for nFrags and TSS Enrichment
- Plots QC scores 
    - Density plot for nFrags vs TSS Enrichment
    - Ridge and violin plots for nFrags and TSS Enrichment individually plotted)



#Create Project 1

ArrowFiles <- list.files("/project/eon/nboldon/MCS2023/fragAnalysis_TSS7", pattern=".arrow$")
ArrowFiles

projMCS1 <- ArchRProject(
	ArrowFiles = ArrowFiles,
	outputDirectory = "/project/eon/nboldon/MCS2023/fragAnalysis_TSS7",
	copyArrows = TRUE
	)

#Add Doublet Scores
projMCS1 <- addDoubletScores(
	input = projMCS1,
	k = 10, #Refers to how many cells near a "pseudo-doublet" to count
	knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection
	LSIMethod = 1
	)

projMCS1$Sample
#The sampleNames of each sample were obtained from Arrow files
#A matrix called sampleColData contains data associated with each sample
#A matrix called cellColData contains data associated with each cell

#To check how much memory size is used to store the ArchRProject in memory wihtin R
paste0("Memory Size = ", round(object.size(projMCS1) / 10^6, 3), " MB")

#To access the cell names associated with each cell
head(projMCS1$cellNames)

#To access the sample names associated with each cell
head(projMCS1$Sample)

getAvailableMatrices(projMCS1)

#To access TSS Enrichment Scores for each cell
quantile(projMCS1$TSSEnrichment)

#To subset the project numerically (ex: 1st 100 cells in project)
# projMCS1[1:100]

#To subset the project based on certain cell names
# projMCS1[projMCS1$cellNames[1:100, ]

#To subset the project to keep all cells corresponding to a specific sample
# idxSample <- BiocGenerics::which(projMCS1$Sample %in% "scATAC_BMMC_R1")
# cellsSample <- projMCS1$cellNames[idxSample]
# projMCS1[cellsSample, ]

#To subset the project to only keep cells that meet a specific cutoff for the TSS enrichment score
# idxPass <- which(projMCS1$TSSEnrichment >= 10)
# cellsPass <- projMCS1$cellNames[idxPass]
# projMCS1[cellsPass, ]

#Plot QC metrics: Log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(projMCS1, select = c("log10(nFrags)", "TSSEnrichment"))
df

#Plot the Log10(Unique Fragments) by TSS enrichment score to identify high quality cells
p <- ggPoint(
	x = df[,1],
	y = df[,2],
	colorDensity = TRUE,
	continuousSet = "sambaNight",
	xlabel = "Log10 Unique Fragments",
	ylabel = "TSS Enrichment",
	xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
	ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p

#To save an editable vectorized version of the plot
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projMCS1, addDOC = FALSE)

#Create ridge plots for each sample using TSS enrichment scores
p1 <- plotGroups(
	ArchRProj = projMCS1, 
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "ridges"
	)
p1

#Create violin plots for each sample using TSS enrichment scores
p2 <- plotGroups(
	ArchRProj = projMCS1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
	)
p2

#Create ridge plots for each sample using log10(Unique Fragments)
p3 <- plotGroups(
	ArchRProj = projMCS1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "ridges"
	)
p3

#Create violin plots for each sample using log10(Unique Fragments)
p4 <- plotGroups(
	ArchRProj = projMCS1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
	)
p4

#Save editable vectorized versions of the plots
plotPDF(p1,p2,p3,p4, name = "TSS7_QC-MCS1.pdf", ArchRProj = projMCS1, addDOC = FALSE, width = 8, height = 8)

################################################################
##Does not run

#Plot fragment size distributions
p5 <- plotFragmentSizes(ArchRProj = projMCS5)
p5

#Plot TSS enrichment profiles
p6 <- plotTSSEnrichment(ArchRProj = projMCS5)
p6

#Save an editable vectorized versions of the plots
plotPDF(p5,p6, name = "QC_FragSize-Distro_2024-02-29.pdf", ArchRProj = projMCS5, addDOC = FALSE, width = 5, height = 5)
#########################################################################################################################

#Save the ArchRProject
saveArchRProject(ArchRProj = projMCS1, outputDirectory = "/project/eon/nboldon/MCS2023/Save-ProjMCS1", load = FALSE)

#This process does not automatically update the ArchRProject object that is active in your current R session. Specifically, the object named ProjMCS1 in the current R session will still point to the original location of the arrow files, not the copied arrow files that reside in the specified outputDirectory. If we wanted to do this, we would specify load = TRUE which causes the saveArchRProject() function to return the saved ArchRProject object which you can assign to overwrite the original ArchRProject object using <-. This effectively saves and loads the ArchRProject from its new location. 
