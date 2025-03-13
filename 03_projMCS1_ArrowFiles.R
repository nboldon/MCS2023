## 03_projMCS1_ArrowFiles.R
- Creates arrow files in ArchR to create project



#Setup an interactive session
salloc --account=eon -t 1-00:00 --mem=128G --nodes=2 --ntasks-per-node=16

#Load required dependencies
module load miniconda3/4.12.0
conda activate /pfs/tc1/project/eon/archr_env

#Load libraries
R
library(ArchR)

#Create arrow files
setwd("/project/eon/nboldon/MCS2023/Analysis")

inputFiles <- list.files("/project/eon/nboldon/MCS2023/Analysis", pattern=".bam$")
namesFiles <- gsub(".bam", "", inputFiles)

addArchRGenome("mm10")
addArchRThreads(threads = 8)

#If chromosomes don't have chr prefixes as expected by ArchR, run the following:
# addArchRChrPrefix(chrPrefix = FALSE)

ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles,
	sampleNames = namesFiles,
	minTSS = 10,
	minFrags = 1000,
	addTileMat = TRUE,
	addGeneScoreMat = TRUE,
	bcTag = 'CB'
)
