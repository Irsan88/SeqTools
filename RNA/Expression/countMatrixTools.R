# to convert counts to (log)FPKM use rpkm()-method from edgeR
# ...

# filter genes with low counts
removeLowExpressionGenes <- function(x,percentageSamples=0.8,countThreshold=10){
  if(missing(x)){
    stop("Provide count matrix")
  }
  if(class(x)!="matrix"){
    stop("x should be a count matrix")
  }
  numberSamples <- ncol(x)
  numberSamplesThreshold <- numberSamples * percentageSamples
  bad.genes <- apply(x,1,function(x){sum(x < countThreshold)}) > numberSamplesThreshold
  good.genes <- !bad.genes
  return(x[good.genes,])
}

# Define clustering functions
clustAverage <- function(x) hclust(x, method="average") 
clustComplete <- function(x) hclust(x,method="complete")
clustWard <- function(x) hclust(x,method="ward.D")
distPearson <- function(x) as.dist(1-cor(t(x), method="pearson"))
distEuclidean <- function(x) dist(x,method = 'euclidean')
distMax <- function(x) dist(x,method = 'maximum')
distManhattan <- function(x) dist(x,method = 'manhattan')
distCanberra <- function(x) dist(x,method = 'canberra')
distBinary <- function(x) dist(x,method = 'binary')
distMinkowski <- function(x) dist(x,method = 'minkowski')

#  cluster the columns of a matrix, can be applied to normalized expression matrices
dendrogramSamples <- function(x,distanceMethod=distPearson,clusterMethod=clustWard){
  if(missing(x)){
    stop("Provide count matrix")
  }
  if(class(x)!="matrix"){
    stop("x should be a count matrix")
  }
  # by default, clustering/distances are calculated on genes of a matrix, 
  # so for samples, transpose the matrix so that columns become rows
  x <- t(x)
  distances <- distanceMethod(x)
  clusters <- clusterMethod(distances)
  dd <- as.dendrogram(clusters)
  return(dd)
}

pcaOnSamples <- function(x){
  if(missing(x)){
    stop("Provide count matrix")
  }
  if(class(x)!="matrix"){
    stop("x should be a count matrix")
  }
  x <- t(x)
  pcaResult <- prcomp(x,center=T,scale=T)
  return(pcaResult)
}

# make a 3d plot from a prcomp-object
plotPCA <- function(x,boxAngle=45,showLabels=T,labelOffset=3,labelAngle=45,labelCex=1,colors=NULL,...){
	stopifnot(class(x)=="prcomp")
	pca <- x
	if(missing(colors)){
		colors <- rep("black",nrow(pca$x))
	}	
  pc1.var <- round((pca$sdev[1]^2 / sum(pca$sdev^2)) * 100,digits=2)
	pc2.var <- round((pca$sdev[2]^2 / sum(pca$sdev^2)) * 100,digits=2)
	pc3.var <- round((pca$sdev[3]^2 / sum(pca$sdev^2)) * 100,digits=2)
	require(scatterplot3d)
	scatter <- scatterplot3d(
		pca$x[,1:3],	
		angle=boxAngle,
		color = colors,
    xlab = paste("PC1 (",pc1.var," %)",sep=""),
    ylab = paste("PC2 (",pc2.var," %)",sep=""),
		zlab = paste("PC3 (",pc3.var," %)",sep=""),
		...
	)
	if(showLabels){
		text(
			scatter$xyz.convert(pca$x[,1], pca$x[,2], pca$x[,3]), 
			labels=rownames(pca$x),
			pos = labelOffset,
			srt = labelAngle,
			cex = labelCex,
			col = colors
		) 	
	}	
	return(scatter)
}

# get the N most variable genes
getMostVar <- function(x,n=250,stat=sd){
	# calculate std dev for all genes
	sds <- apply(x,1,stat)
	# sort descending
	sds <- sds[order(sds,decreasing=T)]
	# get the names of the N most variable genes
	most.var.genes <- names(sds)[1:n]
	return(most.var.genes)
}

compareCountMatrix <- function(x,y,plotDir="correlationPlots",labels=c("m1","m2")){
  require(reshape)
  require(ggplot2)
  m1 <- x
  m2 <- y
  stopifnot(colnames(m1)==colnames(m2))
  # get correlations
  cors <- numeric()
  dir.create(plotDir,showWarnings=F)
  labels <- gsub("\\s","\\.",labels)
  for(sample in colnames(m1)){
    print(sample)
    d <- data.frame(m1=m1[,sample],m2=m2[,sample])
    colnames(d) <- labels
    d <- log2(d+0.01)
    cors[sample] <- round(cor(d[,1],d[,2],method="pearson"),digits=3)
    png(paste(plotDir,"/",sample,"-correlationPlot.png",sep=""),700,700)
    p <- ggplot(d) + 
          geom_point(aes_string(x=colnames(d)[1],y=colnames(d)[2]),pch=21) + 
          labs(title=paste("Correlation",cors[sample]))
    print(p)
    dev.off()
  }  
  return(cors)
}



