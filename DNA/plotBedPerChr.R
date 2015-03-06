# this file takes as input a directory and searches for files 
# with the provided filename pattern. These files should be 5-column 
# bed files where the 5th column contains numeric data. For each file
# a pdf will be made containing plots per chromosome of the 5th column
# chromosome names should contain "chr"-prefix, columns tab-seperated,
# no headers. Order should be: chr, start, end, anything, value

# get command line args
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2){
  message("No command arguments detected, default settings will be used")
  inputDir <- "./"
  searchPattern <- "*.bed"
} else {
  inputDir <- as.character(args[1])
  searchPattern <- as.character(args[2])
}

# for debugging
print(inputDir)
print(searchPattern)

# load dependendies
require(ggplot2)

# define chrs
useChromosomes <- paste("chr",c(1:22,"X","Y"),sep="")

# define file inputs
inputFiles <- list.files(path=inputDir,pattern=searchPattern,full.names=T,recursive=T)
stopifnot(length(inputFiles) > 0)
print(inputFiles)

# import data
inputDataList <- list()
for(file in inputFiles){
  inputDataList[[file]] <- read.table(file,sep="\t",header=F) 
}
chromosomesDetected <- unique(sapply(inputDataList,function(x){x[,1]}))
numberChromosomesDetected <- sum(chromosomesDetected %in% useChromosomes)
stopifnot(numberChromosomesDetected > 0)

# make sure chromosomes are in nice order
inputDataList <- lapply(inputDataList,function(x){
  x$V1 <- factor(x$V1,levels=useChromosomes)
  x <- droplevels(x)
  return(x)
})

# for each bed-file, plot the value column (5th) along genomic
# coordinates for each chromosome seperately
for(file in names(inputDataList)){
  data <- inputDataList[[file]]
  name <- gsub("(.*)\\..*","\\1",file)
  # create pdf file
  pdf(paste(name,".pdf",sep=""),width=10,height=3,onefile=T)
  # loop over all chrs
  for(chr in chromosomesDetected){
    tmp <- subset(data,V1==chr)
    totalSnps <- nrow(tmp)
    p <- ggplot(tmp) + 
      geom_point(aes(x=V2,y=V5)) + 
      scale_x_continuous(breaks=seq(0,250e6,by=10e6),labels=seq(0,250,by=10)) +
      labs(x="Genomic Position (Mbps)",y="VAF",title=paste(chr,": ",totalSnps," SNPs",sep=""))
    print(p)
  }
  # save pdf
  dev.off()
}

q(save="no")
