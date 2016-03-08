segmentDifferential <- function(originalCNAsegments,...){
  require(DNAcopy)
  # check that input format is data frame with 6 columns names 
  # "ID"        "chrom"     "loc.start" "loc.end"   "num.mark"  "seg.mean" 
  # check classes of columns
  # get all unique sample IDs
  allSamples <- unique(originalCNAsegments$ID)
  numberAllSamples <- length(allSamples)
  currentSampleIndex <- 1
  lstSegmentedDiff <- lapply(
    allSamples,
    function(currentSample){
      currentSegments <- subset(originalCNAsegments,ID==currentSample)
      absDiffSegMean <- abs(diff(currentSegments$seg.mean))
      numberMarkers <- length(absDiffSegMean)
      # make DNAcopy input
      currentCna <- CNA(
        genomdat=absDiffSegMean,
        chrom = rep("noChrom",times = numberMarkers),
        maploc = 1:numberMarkers,
        data.type = "logratio",
        sampleid = currentSample
      )
      # segment
      currentSegmentedDiff <- segment(
        x = currentCna,
        min.width = 5,
        ...
      )
      return(currentSegmentedDiff)
    }
  )
  names(lstSegmentedDiff) <- allSamples
  return(lstSegmentedDiff)
}

collectSegmentDetailsFromListSegmentedDiff <- function(lstSegmentedDiff,withPValues=T){
  require(DNAcopy)
  stopifnot(class(lstSegmentedDiff)=="list")
  df <- lapply(lstSegmentedDiff,function(data){
    if(withPValues){
      df <- segments.p(data)
    } else {
      df <- data$output
    }
  })
  df <- do.call("rbind",df)
  rownames(df) <- NULL
  return(df)
}

getUniqueSamplesFromSegments <- function(segments){
  stopifnot("ID" %in% colnames(segments))
  uniqueSamples <- sort(unique(as.character(segments$ID)))
  return(uniqueSamples)
}

validateSegmentFormat <- function(segments){
  print("To be implemented. Better: define segments object")
}

getOriginalSegmentsFromSegmentedDiff <- function(segmentedDiff,originalSegments){
  # Check that unique sample IDS in orignal and segmentedDiff tables are identical
  samplesInOriginalSegments <- getUniqueSamplesFromSegments(originalSegments)
  samplesInSegmentedDiff <- getUniqueSamplesFromSegments(segmentedDiff)
  stopifnot(all(samplesInSegmentedDiff %in% samplesInOriginalSegments))
  # container for results
  allOriginalSegmentsPerSample <- data.frame()
  for(segmentRowId in rownames(segmentedDiff)){ # could also loop over samplesInOriginalSegments, should be same
    currentSegment <- segmentedDiff[segmentRowId,]
    currentSample <- currentSegment[,"ID"]
    currentOriginalSegments <- subset(originalSegments,ID==currentSample)
    rownames(currentOriginalSegments) <- NULL
    # get the indexes of the original segments of this particular segmentedDiff
    originalSegmentsIndexLefToRight <- currentSegment[,c("loc.start","loc.end")]
    leftMostSegmentIndex <- originalSegmentsIndexLefToRight[,"loc.start"]
    rightMostSegmentIndex <- originalSegmentsIndexLefToRight[,"loc.end"] + 1
    foundOriginalSegments <- currentOriginalSegments[leftMostSegmentIndex:rightMostSegmentIndex,]
    allOriginalSegmentsPerSample<- rbind(allOriginalSegmentsPerSample,foundOriginalSegments)
  }
  return(allOriginalSegmentsPerSample)
}

matchSegments <- function(segments.A,segments.B,whichColumns){ 
  stopifnot(all(colnames(segments.A) == colnames(segments.B)))
  # determine what columns are to be used for mathching, by default, all
  if(missing(whichColumns)){
    whichColumns <- colnames(segments.A)
  } else {
    stopifnot(all(whichColumns %in% c(colnames(segments.A),colnames(segments.B))))
  }
  # first make sure everyting is text
  segments.A <- sapply(segments.A[,whichColumns],as.character)
  segments.B <- sapply(segments.B[,whichColumns],as.character)
  # define function to merge values
  collapser <- function(segment){paste(segment,collapse = ";")}
  # create ids  
  a <- apply(segments.A,1,collapser)
  b <- apply(segments.B,1,collapser)
  # match
  return(which(a %in% b))
}

ggSegmentation <- function(segments,highlightSegments,yMinMax,scaleFreedom="free_x",spaceFreedom="free_x",...){
  require(ggplot2)
  require(grid)
  if(missing(yMinMax)){
    yMinMax <- range(segments$seg.mean) * 1.1
  }
  # reduce outliers
  segments[segments$seg.mean > max(yMinMax),"seg.mean"] <- max(yMinMax)
  segments[segments$seg.mean < min(yMinMax),"seg.mean"] <- min(yMinMax)
  # color the segments in the chromorexis region or other region of interest
  if(!missing(highlightSegments)){
    indexSegmentsToHighlight <- matchSegments(segments,highlightSegments)
    segments$highlightMe <- "no"
    segments$highlightMe[indexSegmentsToHighlight] <- "yes"
    # further make plot
    p <- ggplot(segments) +
      geom_segment(aes(x=loc.start,xend=loc.end,y=seg.mean,yend=seg.mean,col=highlightMe)) + 
      scale_color_manual(values = c("black","red"))
  } else { # do not add highlight
    p <- ggplot(segments) +
      geom_segment(aes(x=loc.start,xend=loc.end,y=seg.mean,yend=seg.mean))
  }
  p <- p + facet_grid(ID ~ chrom,scales = scaleFreedom,space = spaceFreedom) + 
    labs(x="Genomic position",y="Copy number estimate") +
    coord_cartesian(ylim = yMinMax) + 
    theme_bw() + 
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position="none",
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.margin = unit(0,units = "cm")
    )
  return(p)
}

whichAreShortAndHighLevel <- function(segments){
  df <- segments[,c("num.mark","seg.mean")]
  df$seg.mean <- abs(df$seg.mean)
  df$num.mark <- log10(df$num.mark + 0.1)
  # if there are only 2 segments in the data frame than you 
  # don't need k-means. The highest amplitude is chro
  if(nrow(df) < 2){
    stop("Fewer than 2 original segments in maximum segmented diff")
  } else if(nrow(df)==2){
    result <- abs(df$seg.mean)==max(abs(df$seg.mean))
  } else {
    kmeansResults <- kmeans(df,centers = 2)
    # find out what K is the clusters with the focal amplicons
    # the rule is that the group with highest ratio of amplitude/length
    ratioSizeAmp <- kmeansResults$centers[,"seg.mean"] / kmeansResults$centers[,"num.mark"]
    highestAndSmallest.cluster = order(ratioSizeAmp,decreasing = T)[1]
    result <- kmeansResults$cluster == highestAndSmallest.cluster
  }
  return(result)
}

getSegmentStats <- function(segments){
  # define amps and lenghts
  segmentAmplitudes <- segments$seg.mean
  segmentLength <- log10((segments$loc.end - segments$loc.start) + 0.1)
  # get stats
  numberSegments <- nrow(segments)
  numberUniqueChroms <- length(unique(segments$chrom))
  maxAmplitude <- max(segmentAmplitudes)
  minAmplitude <- min(segmentAmplitudes)
  meanAmplitude <- mean(segmentAmplitudes)
  medianAmplitude <- median(segmentAmplitudes)
  sdAmplitude <- sd(segmentAmplitudes)
  maxLength <- max(segmentLength)
  minLenght <- min(segmentLength)
  meanLength <- mean(segmentLength)
  medianLenght <- median(segmentLength)
  sdLength <- sd(segmentLength)
  results <- c(
    numberSegments,
    numberUniqueChroms,
    maxAmplitude,
    minAmplitude,
    meanAmplitude,
    medianAmplitude,
    sdAmplitude,
    maxLength,
    minLenght,
    meanLength,
    medianLenght,
    sdLength
  )
  # var names
  names(results) <- c(
    "numberSegments",
    "numberUniqueChroms",
    "maxAmplitude",
    "minAmplitude",
    "meanAmplitude",
    "medianAmplitude",
    "sdAmplitude",
    "maxLength",
    "minLenght",
    "meanLength",
    "medianLenght",
    "sdLength"
  )
  return(results)
}

getFeaturesPerSegmentedDiff <- function(segmentedDiff,originalSegments){
  # loop over all segmented diffs, extract the original CN segments 
  # and return features about the size and amplitude of both base segments and 
  # altered segments in that respective segmented diff. You and up with a matrix
  # where each row represents a segmented diff and each column with a feature. 
  # Based on particular features, each segmented diff will be tested whether
  # this looks like chromorexis (or something else useful, like chromothrypsis, firestorm, ...)
  
  # get unique samples names
  samplesInOriginalSegments <- getUniqueSamplesFromSegments(originalSegments)
  samplesInSegmentedDiff <- getUniqueSamplesFromSegments(segmentedDiff)
  stopifnot(all(samplesInSegmentedDiff %in% samplesInOriginalSegments))
  
  # define matrix of features
  featureMatrixPerSegmentedDiff <- matrix(nrow = nrow(segmentedDiff),ncol = 27) #aai, use sapply to auto dimension the feature matrix
  index <- 1
  labels <- c()
  featureLabels <- c()
  # loop over samples
  for(currentSample in samplesInSegmentedDiff){
    print(currentSample)
    currentSegmentedDiffs <- subset(segmentedDiff,ID==currentSample)
    rownames(currentSegmentedDiffs) <- NULL
    currentOriginalSegments <- subset(originalSegments,ID==currentSample)
    rownames(currentOriginalSegments) <- NULL
    # loop over segmented diffs in current samples
    for(currentSegmentName  in rownames(currentSegmentedDiffs)){
      currentSegment <- currentSegmentedDiffs[currentSegmentName,]
      currentLabel <- paste(currentSample,currentSegmentName,sep=".")
      # get basic segmented diff stats
      stopifnot(c("ID","loc.start","loc.end","seg.mean","pval") %in% colnames(currentSegment))
      segmentedDiffAmp <- currentSegment[,"seg.mean"]
      names(segmentedDiffAmp) <- "segmentedDiffAmp"
      segmentedDiffPval <- currentSegment[,"pval"]
      names(segmentedDiffPval) <- "segmentedDiffPval"
      # get original CN segments in this current segmented diff
      currentOriginalSegmentsInCurrentSegmentedDiff <- getOriginalSegmentsFromSegmentedDiff(currentSegment,currentOriginalSegments)
      # use K-means to seperate the small high-levels form long-low level segments
      binShortAmp <- whichAreShortAndHighLevel(currentOriginalSegmentsInCurrentSegmentedDiff)
      alteredSegments <- currentOriginalSegmentsInCurrentSegmentedDiff[binShortAmp,]
      baseSegments <- currentOriginalSegmentsInCurrentSegmentedDiff[!binShortAmp,]
      # extract features
      baseSegmentFeatures <- getSegmentStats(baseSegments) # add function to retrieve number of unique chroms
      alteredSegmentFeatures <- getSegmentStats(alteredSegments)
      # calculate jump height (amplitude difference bewteen base and altered segments)
      jumpHeight <- alteredSegmentFeatures["meanAmplitude"] - baseSegmentFeatures["meanAmplitude"]
      # change labels
      names(baseSegmentFeatures) <- paste("baseSegments",names(baseSegmentFeatures),sep = ".")
      names(alteredSegmentFeatures) <- paste("alteredSegments",names(alteredSegmentFeatures),sep = ".")
      names(jumpHeight) <- "jumpHeight"
      # combine features for current segmented diff
      currentFeatures <- c(segmentedDiffAmp,segmentedDiffPval,jumpHeight,baseSegmentFeatures,alteredSegmentFeatures)
      # add to feature matrix
      featureMatrixPerSegmentedDiff[index,] <- currentFeatures
      # save labels for dimnames matrix
      labels[index] <- currentLabel
      featureLabels <- names(currentFeatures)
      # increment index
      index <- index + 1
    }
  }
  # add row and column names to matrxi
  rownames(featureMatrixPerSegmentedDiff) <- labels
  colnames(featureMatrixPerSegmentedDiff) <- featureLabels
  return(featureMatrixPerSegmentedDiff)
}

