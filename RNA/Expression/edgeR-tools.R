# description: wrapper for anova on N groups and all pairwise comparisons

# to-do: when the group factor only contains two levels, no
# one-way anova is needed (is the same as two-samples t-test)

compareNGroupsEdgeR <- function(y,model="group",normalization="TMM"){
	if(missing(y)){
		stop("Provide DGElist object")
	}
	stopifnot(class(y)=="DGEList")
	# the DGEList should containt y$samples with one column
	# named groups which is the factor for the analysis
	
	model <- paste("~0 + ",model,sep="")
	model <- as.formula(model)
	design <- model.matrix(model,data=y$samples)
	colnames(design)  <- gsub("^group(.*$)","\\1",colnames(design))
	# get all pairwise combinations in vector
	combinationsMatrix <- combn(levels(y$samples$group),2)
	comparisonsCharacter <- apply(combinationsMatrix,2,function(x){paste(x[c(2,1)],collapse=" - ")})
	# make all contrasts
	pairwise.contrasts <- makeContrasts(contrasts=eval(comparisonsCharacter),levels=design)

	# normalize counts and fit general model
	y <- calcNormFactors(y,method=normalization)
	y <- estimateGLMCommonDisp(y,design)
	y <- estimateGLMTrendedDisp(y,design)
	y <- estimateGLMTagwiseDisp(y,design)
	fit <- glmFit(y,design)
	
	# get results from the fit
	top.tables <- list()
	top.tables[["OneWayAnova"]] <- topTags(glmLRT(fit,contrast=pairwise.contrasts),n=nrow(y))$table
	for(contrast in colnames(pairwise.contrasts)){
		top.tables[[contrast]]  <- topTags(glmLRT(fit,contrast=pairwise.contrasts[,contrast]),n=nrow(y))$table
	}
	return(top.tables)
}

# description: inspect counts for a gene per sample
plotGenes <- function(y){
	require(ggplot2)
	require(reshape)
	x <- cpm(y,log=T,normalized.lib.sizes=T)
	x <- melt(x)
	colnames(x) <- c("gene","sample","logCPM")
	x$group <- y$samples[as.character(x$sample),"group"]
	print(head(x))
	plot <- ggplot(x) + 
		geom_jitter(aes(x=group,y=logCPM,color=group)) +
		scale_color_brewer(type="qual",palette="Set1") +
#		position_jitter(width=0.1) +
		facet_wrap(~gene,scales="free_y")
	return(plot)
}
