# description: wrapper for anova on N groups and all pairwise comparisons

# to-do: by default all pairwise comparisons are analyzed. Make it possible
# for the user to provide preset comparisons himself

compareNGroupsEdgeR <- function(y,model="group",normalization="TMM",BCV){
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
	
  if(!missing(BCV)){
    # plot BCV
    png(BCV,15,15,"cm",res=300)
    plotBCV(y)
    dev.off()    
  }

	# get results from the fit
	top.tables <- list()
  if(ncol(pairwise.contrasts) > 1){
	  top.tables[["OneWayAnova"]] <- topTags(glmLRT(fit,contrast=pairwise.contrasts),n=nrow(y))$table
  }
	for(contrast in colnames(pairwise.contrasts)){
    print(contrast)
		top.tables[[contrast]]  <- topTags(glmLRT(fit,contrast=pairwise.contrasts[,contrast]),n=nrow(y))$table
	}
	return(top.tables)
}

anovaLikeEdgeR <- function(y,model="group",normalization="TMM",BCV){
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
  
  if(!missing(BCV)){
    # plot BCV
    png(BCV,15,15,"cm",res=300)
    plotBCV(y)
    dev.off()    
  }
  
  # get results from the fit
  top.table <- topTags(glmLRT(fit,contrast=pairwise.contrasts),n=nrow(y))$table
  return(top.table)
}

# description: inspect counts for a gene per sample
plotGenes <- function(y,scales="free_y",base=12,minThreshold,extraIdColumn,numberColumns){
	require(ggplot2)
	require(reshape2)
  order <- rownames(y)
	x <- cpm(y,log=T,normalized.lib.sizes=T)
	x <- melt(x)
	colnames(x) <- c("gene","sample","logCPM")
  x$gene <- factor(x$gene,levels=order)
	x$group <- y$samples[as.character(x$sample),"group"]
  if(!missing(extraIdColumn)){
    x$symbol <- y$genes[as.character(x$gene),extraIdColumn]
  } else{
    x$symbol <- ""
  }
  if(!missing(minThreshold)){
    x$logCPM[x$logCPM < minThreshold] <- minThreshold
  }
	plot <- ggplot(x) + 
		geom_jitter(aes(x=group,y=logCPM,color=group),position=position_jitter(width=0.2)) +
		scale_color_brewer(type="qual",palette="Set1") +
		facet_wrap(~gene+symbol,scales=scales) +
    theme_gray(base_size=base) + 
    theme(legend.position="top")
  if(!missing(numberColumns)){
    plot <- plot + facet_wrap(~gene+symbol,scales=scales,ncol=numberColumns)
  }
	return(plot)
}

# description: inspect counts for a gene per sample
jitterPlotGenes <- function(y,scales="free_y",base=12,minThreshold,extraIdColumn,numberColumns){
  require(ggplot2)
  require(reshape2)
  order <- rownames(y)
  x <- cpm(y,log=T,normalized.lib.sizes=T)
  x <- melt(x)
  colnames(x) <- c("gene","sample","logCPM")
  x$gene <- factor(x$gene,levels=order)
  x$group <- y$samples[as.character(x$sample),"group"]
  if(!missing(extraIdColumn)){
    x$symbol <- y$genes[as.character(x$gene),extraIdColumn]
  } else{
    x$symbol <- ""
  }
  if(!missing(minThreshold)){
    x$logCPM[x$logCPM < minThreshold] <- minThreshold
  }
  plot <- ggplot(x) + 
    geom_jitter(aes(x=group,y=logCPM,color=group),position=position_jitter(width=0.2)) +
    scale_color_brewer(type="qual",palette="Set1") +
    facet_wrap(~gene+symbol,scales=scales) +
    theme_gray(base_size=base) + 
    theme(legend.position="top")
  if(!missing(numberColumns)){
    plot <- plot + facet_wrap(~gene+symbol,scales=scales,ncol=numberColumns)
  }
  return(plot)
}

dotPlotGenes <- function(y,scales="free_y",base=12,minThreshold,extraIdColumn,numberColumns){
  require(ggplot2)
  require(reshape2)
  order <- rownames(y)
  x <- cpm(y,log=T,normalized.lib.sizes=T)
  x <- melt(x)
  colnames(x) <- c("gene","sample","logCPM")
  x$gene <- factor(x$gene,levels=order)
  x$group <- y$samples[as.character(x$sample),"group"]
  if(!missing(extraIdColumn)){
    x$symbol <- y$genes[as.character(x$gene),extraIdColumn]
  } else{
    x$symbol <- ""
  }
  if(!missing(minThreshold)){
    x$logCPM[x$logCPM < minThreshold] <- minThreshold
  }
  plot <- ggplot(x) + 
    geom_dotplot(aes(x=group,y=logCPM,fill=group),binaxis = "y",stackdir="center") +
    scale_fill_brewer(type="qual",palette="Set1") +
    facet_wrap(~gene+symbol,scales=scales) +
    theme_gray(base_size=base) + 
    theme(legend.position="top")
  if(!missing(numberColumns)){
    plot <- plot + facet_wrap(~gene+symbol,scales=scales,ncol=numberColumns)
  }
  return(plot)
}
