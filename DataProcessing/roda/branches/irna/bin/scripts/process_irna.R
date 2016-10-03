#!/usr/bin/Rscript

###########################################
# Parse command line arguments
###########################################

args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Usage:./processHTseqs.R analysisDir outputDir")
} else {
        inputDir <- args[1]
	outputDir <- args[2]
}

###########################################
# Collect the read counts for all samples
###########################################

hts <- list.files(
	path = inputDir,
	pattern = "\\.ssv$",
	full.names = TRUE,
	recursive = TRUE
)

print(hts)

counts <- list()
for(sample in hts){
	counts[[sample]] <- read.table(
		file = sample, 
		header = FALSE,
		sep = "\t",
		row.names = 1
	)
	name <- basename(sample)
	name <- gsub(name,pattern="\\.htSeqCounts_ensemblID\\.ssv",replacement="")
	colnames(counts[[sample]]) <- name
}
counts <- do.call("cbind",counts)

goodRows <- grep("^__",rownames(counts),invert=T)
counts <- counts[goodRows,]

# save to R-data file
dir.create(outputDir)
setwd(outputDir)
save(counts,file="allSamples_readCounts.RData")

# save to tab seperate data table
nc <- ncol(counts)
counts$id <- rownames(counts)
counts <- counts[,c(nc+1,1:nc)]
options(scipen=99)
write.table(
	x = counts,
	file = "allSamples_readCounts.tsv",
	sep = "\t",
	row.names = FALSE,
	col.names = TRUE,
	quote=FALSE
)
