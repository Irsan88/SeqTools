#!/usr/bin/Rscript

# check dependencies
require(ggplot2)
require(reshape)

# Get command line arguments
# maybe change for getopt
args<-commandArgs(TRUE)
if(length(args) != 2){
	stop("Usage: ./counts2rpkm.R inputDir windows.bed")
} else {
	input.dir <- args[1]
	windowsFile <- args[2]
  minmax <- as.numeric(args[3])
}

# collect counts for all samples in the
# provided input directory
samples <- list.files(
	input.dir,
	pattern="*.count$",
	full.names=T
)

# import regions in bed file
counts <- read.table(windowsFile,header=F,sep="\t")
colnames(counts) <- c("chr","start","end","symbol")

# import data for each .count file
for(file in samples){
	basename <- basename(file)	
	counts[,basename] <- read.table(file,header=F)
}

# assume the sample name can be inferred
# from everything before the _ character
# example: VU1456-hiseq-13Nov_Sorted.recal.count
colnames(counts) <- gsub("(.*)_.*",replacement="\\1",x=colnames(counts))

# calculate library sizes per million
# mapped reads, called rpm
lib.sizes <- apply(counts[,5:ncol(counts)],2,sum) / 1e6
rpm <- t(t(counts[,5:ncol(counts)]) / lib.sizes)

# calculate the baseline rpm based on the median
# rpm value for each windows
baseline <- apply(rpm,1,median)

# calculate the log ratios (these are your 
# copy number estimates!
lrr <- rpm / baseline
# identify regions with bad baseline
good.windows <- which(baseline>0)
lrr <- lrr[good.windows,]
lrr <- log2(lrr)
# annotate the windows with chr start end
lrr <- cbind(counts[good.windows,1:4],lrr)

# reshape the data so that you can plot
# with famous ggplot2
lrr <- melt(lrr,id.vars=colnames(lrr)[1:4])
colnames(lrr)[5:6] <- c("sample","lrr")
# apply min and max values
lrr[which(lrr$lrr > 3),"lrr"] <- minmax
lrr[which(lrr$lrr < -3),"lrr"] <- minmax

# and export the log-r-ratios data frame
write.csv(lrr,file="lrr-values.csv",row.names=F,col.names=T)





