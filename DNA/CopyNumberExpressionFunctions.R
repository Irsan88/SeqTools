heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }

        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(colnames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }

    if (!missing(ColSideColors)) {

        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }

        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

# Define clustering functions
clustAverage <- function(x) hclust(x, method="average") 
clustComplete <- function(x) hclust(x,method="complete")
clustWard <- function(x) hclust(x,method="ward.D2")
distPearson <- function(x) as.dist(1-cor(t(x), method="pearson"))
distEuclidean <- function(x) dist(x,method = 'euclidean')
distMax <- function(x) dist(x,method = 'maximum')
distManhattan <- function(x) dist(x,method = 'manhattan')
distCanberra <- function(x) dist(x,method = 'canberra')
distBinary <- function(x) dist(x,method = 'binary')
distMinkowski <- function(x) dist(x,method = 'minkowski')

# colors for heatmaps
blueRedCols <- colorRampPalette(c("blue3","blue2","blue1","blue","white","red","red1","red2","red3"))(100)

rna.clustering <- function(x,scale="row",hc=clustComplete,dist=distPearson,minmaxCol=2,palette="RdBu",col,column.colors,row.colors,labRow=F,...){
	require(affy)
	require(RColorBrewer)
	require(gplots)
  br <- seq(-1*minmaxCol,minmaxCol,length.out=41)
  if(missing(col)){
    col <- colorRampPalette(rev(brewer.pal(n=11,name=palette)))(length(br)-1)
  }
	if(missing(column.colors) & missing(row.colors)){
		h <- heatmap.3(
			x = x,
			col = col,
      breaks=br,
			labRow = labRow,
			scale=scale,
      hclustfun = hc,
      distfun = dist,
			...
		)
	} else if(!missing(column.colors) & missing(row.colors)) {
		h <- heatmap.3(
			x = x,
			col = col,
			breaks=br,
			labRow = labRow,
			ColSideColors = column.colors,
			scale=scale,
			hclustfun = hc,
			distfun = dist,
			...
		)
	} else if(missing(column.colors) & !missing(row.colors)) {
		h <- heatmap.3(
			x = x,
			col = col,
			breaks=br,
			labRow = labRow,
			RowSideColors = row.colors,
			scale=scale,
			hclustfun = hc,
			distfun = dist,
			...
		)
	} 
	else if(!missing(column.colors) & !missing(row.colors)) { 
		h <- heatmap.3(
			x = x,
			col = col,
			breaks=br,
			labRow = labRow,
			ColSideColors = column.colors,
			RowSideColors = row.colors,
			scale=scale,
			hclustfun = hc,
			distfun = dist,
			...
		)
	} else {
		stop("Something went wrong ...")
	}
	return(h)
}


makeColSideColors <- function(x){
	require(affy)
	col <- data.frame(row.names=rownames(x))
	for(column in colnames(x)){
		if(is.numeric(x[,column])){
			col[,column] <- numericToColors(x[,column])
		} else {
			col[,column] <- categoriesToColors(x[,column])		
		}
	}
	return(as.matrix(col))
}

categoriesToColors <- function(x){
  require(RColorBrewer)
  if(is.numeric(x)){ 
		stop("You should provide non-numeric vector")
	}
	# determine the amount of categories
	numberColorsRequired <- length(unique(x))
	# number of individuals
	n <- length(x)
	# make default returnColors
	returnColors <- "not set yet"
	if(numberColorsRequired==1){
	  returnColors <- rep("red",n)
	}
	if(numberColorsRequired==2){
	  returnColors <- c("red","black")[factor(x)]
	}
	if(numberColorsRequired > 2 & numberColorsRequired < 10){
	  colorSet <- brewer.pal(n = numberColorsRequired,name = "Set1")
	  returnColors <- colorSet[factor(x)]
	}
	if(numberColorsRequired > 9 & numberColorsRequired < 13){
	  colorSet <- brewer.pal(n = numberColorsRequired,name = "Set3")
	  returnColors <- colorSet[factor(x)]
	}
	if(numberColorsRequired > 12){
	  colorSet <- colorRampPalette(colors = brewer.pal(n = 9,name = "Set1"))(numberColorsRequired)
	  returnColors <- colorSet[factor(x)]
	}
	returnColors[which(is.na(returnColors))] <- "grey"
	return(returnColors)
}

numericToColors <- function(x,na.color="grey"){
	require(gplots)
	if(!is.numeric(x)){ 
		stop("You should provide numeric vector")
	}
	# determine the amount of numbers	
	n <- length(x)
	colors <- rep(na.color,n)
	no.na <- which(!is.na(x))
	require(RColorBrewer)
	colors[no.na] <- colorpanel(length(no.na),low="green",mid="black",high="red")[rank(x[no.na])]
	return(colors)
}

# jitter plot of individual genes
plotProbes <- function(eset,group=NULL,scales="free_y",rows=NULL,jitterWidth=0.1){
	require(ggplot2)
	require(reshape)
	require(affy)
	# get melted expression matrix 
	probe.order <- rownames(eset)
	d <- melt(exprs(eset))
	colnames(d) <- c("probe","sample","expression")
	# preserve probe order, useful when you supply expression-set
	# that is already ordered by p-value
	d$probe <- factor(d$probe,levels=unique(as.character(d$probe)))
	if(missing(group)){
		d$group <- "all"
	}
	else {
		d$group <- pData(eset)[as.character(d$sample),group]
	}
	plot <- ggplot(d) + 
		geom_jitter(aes(x=group,y=expression,color=group),position=position_jitter(width=jitterWidth)) + 
		facet_wrap(~probe,scales=scales,nrow=rows) +
		scale_color_brewer(type="qual",palette="Set1") +
		theme(legend.position="none")
	return(plot)
}

dna.clustering <- function(x,dist=distPearson,clust=clustWard,col=T,column.colors,colors=c("purple","blue","grey","red","black"),labRow=F,...){
	require(RColorBrewer)
	if(missing(column.colors)){	
		h <- heatmap.3(
			x = x,
			col = colors,
			Rowv = F,
			Colv = col,
			keysize=1,
			labRow=labRow,
			distfun=dist,
			hclust=clust,
      dendrogram="column",
			...
		)
	} else {
		h <- heatmap.3(
			x = x,
			col = colors,
			Rowv = F,
			Colv = col,
			keysize=1,
			labRow=labRow,
			distfun=dist,
			hclust=clust,
      dendrogram="column",
			ColSideColors = column.colors,
			...
		)
	}
	return(h)
}

getDNAinstability <- function(x,divideBy=1e6){
  require(CGHregions)
  instability <- numeric()
  for(sample in sampleNames(x)){
    not.normal <- which(regions(x[,sample])!=0)
    instability[sample] <- sum(bpend(x)[not.normal] - bpstart(x)[not.normal]) / divideBy
  }
  return(instability)
}


cghToCrosstab <- function(x,factor,levels=3){
  require(CGHregions)
  f <- factor(pData(x)[,factor])
  ct <- list()
  stopifnot(levels %in% c(3,5))
  for(i in featureNames(x)){
    ct[[i]] <- data.frame(row.names=levels(f))
    if(levels==3){
      for(group in levels(f)){
        ct[[i]][group,"loss"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]<0)
        ct[[i]][group,"normal"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==0)
        ct[[i]][group,"gain"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]>0)
      }
    } else {
        for(group in levels(f)){
          ct[[i]][group,"big.loss"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==-2)
          ct[[i]][group,"loss"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==-1)
          ct[[i]][group,"normal"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==0)
          ct[[i]][group,"gain"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==1)
          ct[[i]][group,"big.gain"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==2)
        }
    }
  }
  return(ct)
}

cghFreqPlot <- function(x,group=NULL,base=12,scales="free_x",alpha=1,ylim100percent=T){
	require(CGHregions)
	require(reshape)
	require(ggplot2)
	require(grid)
	r <- as.data.frame(regions(x))
  r$index <- 1:nrow(r)
  r$chr <- chromosomes(x)
  r$start <- bpstart(x)
  r$end <- bpend(x)
	r <- melt(r,id.vars=c("index","chr","start","end"))
  colnames(r) <- c("index","chr","start","end","sample","call")
	r$call[r$call < -1] <- -1
	r$call[r$call > 1] <- 1
	if(missing(group)){
		r$group <- "all"
    r$group <- as.factor(r$group)
		groupSizes <- ncol(x)
		names(groupSizes) <- "all"
	} else {
		r$group <- pData(x)[as.character(r$sample),group]
		groupSizes <- summary(pData(x)[,group])
	}
	groupNames <- paste(names(groupSizes)," (N=",groupSizes,")",sep="")
	names(groupSizes) <- groupNames
	levels(r$group) <- groupNames
	r <- cast(r,index ~ call + group,fun.aggregate=length,value="sample")
	r <- melt(r,id.vars="index")
	r$chr <- chromosomes(x)[r$index]
	r$start <- bpstart(x)[r$index]
	r$end <- bpend(x)[r$index]
	r$midpoint <- rowMeans(r[,c("start","end")])
	r$size <- r$end - r$start
	r <- subset(r,call != 0)
	r$value <- r$value * r$call
	# get percentages
	r$value <- r$value / groupSizes[as.character(r$group)]	
	# plotting
	p <- ggplot(r) + 
	  geom_bar(aes(x=midpoint,y=value,width=size,fill=factor(call)),stat="identity",position="identity") + 
		facet_grid(group~chr,scales=scales,space=scales) +
		scale_fill_manual(values=c("blue","red")) + 
		theme_bw(base_size=base) + 
		theme(
			legend.position="none",
			axis.ticks.x = element_blank(),
			axis.text.x = element_blank(),
			axis.title = element_blank(),
			panel.grid = element_blank(),
			panel.margin = unit(0,"cm")
		)
	if(ylim100percent){
	  p <- p + ylim(c(-1,1))
	} 
	return(p)
}

cghKaryoHeatmap <-function(x,rows=2,scales="fixed",base=12,flip=F,colors=c("purple","blue","grey","red","black")){
	require(CGHregions)
	require(reshape)
	require(ggplot2)
	require(grid)
	totalSamples <- ncol(x)
	r <- as.data.frame(regions(x))
	r$chr <- chromosomes(x)
	r$size <- bpend(x) - bpstart(x)
	r<-melt(r,id.vars=c("chr","size"))
	colnames(r)[3:4] <- c("sample","call")
	if(flip){
	  p <- ggplot(r) + 
	    geom_bar(aes(x=sample,y=size,fill=factor(call)),stat="identity",width=1) +
	    scale_y_continuous(breaks=seq(0,300e6,10e6),labels=paste(seq(0,300,10),"Mb",sep=" ")) + 
      coord_flip() + 
	    scale_fill_manual(values=colors) + 
	    facet_wrap(~chr,nrow=rows,scales=scales) + 
	    theme_bw(base_size=base) + 
	    theme(
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.ticks.y = element_blank(),
	      axis.text.y = element_blank(),
	      panel.grid = element_blank(),
	      axis.title.y = element_blank()	    ) +
      labs(y=paste("Samples (N=",totalSamples,")",sep=""))
	} else {
	  p <- ggplot(r) + 
	    geom_bar(aes(x=sample,y=size,fill=factor(call)),stat="identity",width=1) +
	    scale_fill_manual(values=colors) + 
	    facet_wrap(~chr,nrow=rows,scales=scales) + 
	    theme_bw(base_size=base) + 
	    theme(
	      legend.position = "top",
	      legend.title = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.text.x = element_blank(),
	      panel.grid = element_blank(),
	      axis.title.y = element_blank()
	    ) + 
	    scale_y_continuous(breaks=seq(0,300e6,10e6),labels=paste(seq(0,300,10),"Mb",sep=" ")) + 
	    labs(x=paste("Samples (N=",totalSamples,")",sep=""))
	}
	return(p)
} 

bedHeatmap <- function(x,fill="green",scales="fixed",rows=2,base=12){
	require(ggplot2)
	colnames(x) <- c("chr","start","end","sample","value")
	x$size <- x$end - x$start
	p <- ggplot(x) +
		geom_bar(aes(x=sample,y=size),fill=fill,stat="identity",width=1) +
		facet_wrap(~chr,nrow=rows,scales=scales) +
		scale_y_continuous(breaks=seq(0,300e6,10e6),labels=paste(seq(0,300,10),"Mb",sep=" ")) +
		theme_bw(base_size=base) + 
		theme(
			legend.position = "top",
			legend.title = element_blank(),
			axis.ticks.x = element_blank(),
			axis.text.x = element_blank(),
			panel.grid = element_blank(),
			axis.title.y = element_blank()
		)
	return(p)
}

bedFreqPlot <- function(x,fill=T,color="green",base=12){
  require(ggplot2)
  require(grid)
  colnames(x) <- c("chr","start","end","group","value")
  x$midpoint <- rowMeans(x[,c("start","end")])
  x$size <- x$end - x$start
  p <- ggplot(x) + 
    geom_bar(aes(x=midpoint,y=value,width=size),fill=color,stat="identity",position="identity") + 
    facet_grid(group~chr,scales="free_x",space="free_x") + 
    theme_bw(base_size=base) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="top",
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.margin = unit(0,"cm")
    ) + 
    labs(
      x="Genomic position",
      y=paste("Frequency")
    )
  return(p)
}

cghSquareHeatmap <- function(x){
  require(CGHregions)
  require(reshape)
  require(ggplot2)
  require(grid)
  d <- as.data.frame(regions(x))
  d$region <- rownames(d)
  d$chr <- chromosomes(x)
  d$start <- bpstart(x)
  d$end <- bpend(x)
  sapply(d,class)
  d <- melt(d,id.vars=c("region","chr","start","end"))
  colnames(d)[5:6] <- c("sample","call")
  d$size <- d$end - d$start
  p <- ggplot(d) + 
    geom_bar(aes(x=sample,y=size,fill=factor(call)),position="stack",stat="identity",width=1) + 
    facet_grid(chr ~ .,scales="free",space="free") + 
    scale_y_reverse() + 
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.position = "bottom",
      panel.margin = unit(0,"cm"),
      plot.margin = unit(c(0,0,0,0),"cm")
    ) + 
    scale_fill_manual(values=c("blue","grey","red")) +
    labs(fill="Copy number call")
  return(p)
}

plotCopyNumberSegmentsandPoint <- function(lrr,segments,callThresholds,dilutionPercentage=0.1,chrList,sampleList,startMin,endMax,pointSize=1,pointColor="grey",pointShape=21,ploidy=F,segmentSize=1){
  # make sure lrr is a data frame with chr, pos, lrrSample1, lrrSample2, lrrSample3
  colnames(lrr)[1:2] <- c("chrom","pos")
  samples <- colnames(lrr)[3:ncol(lrr)]
  # subset chromosomes
  if(!missing(chrList)){
    lrr <- subset(lrr,chrom %in% chrList)
  }
  # subset samples
  if(!missing(sampleList)){
    lrr <- lrr[,c("chrom","pos",sampleList)]
  }
  # subset start
  if(!missing(startMin)){
    lrr <- subset(lrr,pos >  startMin)
  }
  # subset end
  if(!missing(endMax)){
    lrr <- subset(lrr,pos <  endMax)
  }
  # randomly diltue lrr, useful if many points
  keep <- sample(1:nrow(lrr),size=round(dilutionPercentage*nrow(lrr)),replace=F)
  lrr <- lrr[keep,]
  require(reshape)
  lrr <- melt(lrr,id.vars=c("chrom","pos"))
  colnames(lrr)[3:4] <- c("ID","seg.mean")
  # make sure the sample order of plotting is not alfabetic but the original
  # provided order in the lrr data frame
  lrr$ID <- factor(lrr$ID,levels=samples)
  require(ggplot2)
  # convert to ploidy of requested
  if(ploidy==T){
    lrr$seg.mean <- 2*2^(lrr$seg.mean)
    if(!missing(segments)){
    	segments$seg.mean <- 2*2^(segments$seg.mean)
    }
  }
  p <- ggplot(lrr) +
    geom_point(aes(x=pos,y=seg.mean),size=pointSize,col=pointColor,shape=pointShape) +
    facet_grid(ID ~ chrom,space="free_x",scales="free_x")
  if(!missing(segments)){
    # first make sure the samples are in the right order
    segments$ID <- factor(segments$ID,levels=samples)
    # subset chrs
    if(!missing(chrList)){
      segments <- subset(segments,chrom %in% chrList)
    }
    # subset samples
    if(!missing(sampleList)){
      segments <- subset(segments,ID %in% sampleList)
    }
    if(!missing(callThresholds)){
      segments$call <- "unchanged"
      segments$call[segments$seg.mean > callThresholds[1]] <- "gain"
      segments$call[segments$seg.mean < callThresholds[2]] <- "loss"
      p <- p + geom_segment(aes(x=loc.start,xend=loc.end,y=seg.mean,yend=seg.mean,col=call),data=segments,size=segmentSize) + 
        scale_color_manual(values = c("red","blue","black")) + 
        theme(legend.position="top")
    }
    else {
      p <- p + geom_segment(aes(x=loc.start,xend=loc.end,y=seg.mean,yend=seg.mean),col="red",data=segments,size=segmentSize)
    }
  }
  if(!missing(startMin) & !missing(endMax)){
    p <- p + coord_cartesian(xlim=c(startMin,endMax))
  }
  return(p)
}

# Generic copy number frequency plot tool
copyNumberFrequencyPlot <- function(calls,featureData,groups=NULL){
  # calls and featureData is required, groups is optional
  if(missing(calls) | missing(featureData)){
    stop("Provide calls and featureData")
  }
  # featureData should contain chr, start & end
  if(!all(c("chr","start","end") %in% colnames(featureData))){
    stop("featureData should at least contain columns named chr, start & end")
  }
  # the chr column of feature data should be integer or factor
  if(!class(featureData$chr) %in% c("integer","factor")){
    stop("The chr column is not an integer or a sorted character (factor)")
  }
  # the number of features should be equal to the number of rows in calls
  if(nrow(calls) != nrow(featureData)){
    stop("calls and featureData should have equal number of rows")
  }
  # if groups is defined, it should have the same lenght as the number of 
  # columns of calls
  if(!missing(groups) & length(groups)!= ncol(calls)){
    stop("Number of group members should be equal to number of columns in calls")
  }
  x <- cbind(featureData[,c("chr","start","end")],calls)
  require(reshape)
  x <- melt(c,id.vars=c("chr","start","end"))
  colnames(x) <- c("feature","sample","call")
  # add group info to each sample when available
  if(missing(groups)){
    x$group <- "all"
  } else {
    names(groups) <- colnames(calls)
    x$group <- groups[as.character(x$sample)]
  }
  return(x)
}

plotDNAcopySegmentsAll <- function(x,outputDir="segmented-plots",width=15,height=5,ylim,...){
  stopifnot(class(x)=="DNAcopy")
  require(DNAcopy)
  require(Cairo)
  if(missing(ylim)){
    max <- ceiling(max(abs(x$output$seg.mean)))
    if(max > 8){
      max <- 8
    }
    ylim <- c(-1*max,max)
  }
  # retrieve sample anc chr names
  samples <- unique(x$output$ID)
  chromosomes <- unique(x$output$chrom)
  # create output dir
  dir.create(path = file.path(outputDir))
  # export IGV file
  write.table(x$output,file=file.path(outputDir,"segments.seg"),row.names=F,col.names=T,quote=F,sep="\t")
  # save segmented DNAcopy object
  save(x,file=file.path(outputDir,"segments.RData"))
  for(currentSample in samples){
    print(currentSample)
    # create dir for each sample
    sampleDirPath <- file.path(outputDir,currentSample)
    dir.create(sampleDirPath)
    # grab current sample from data
    currentData <- subset(x,samplelist = currentSample)
    # define file name
    wholeGenomeFile <- file.path(sampleDirPath,"whole-genome.pdf")
    # render graphics for whole genome
    CairoPDF(file=wholeGenomeFile,width,height)
    plot(currentData,altCol=T,ylim=ylim,...)
    dev.off()
    # plot each chr
    for(currentChr in chromosomes){
      # define file name
      basename <- paste(currentChr,".pdf",sep="")
      chrFile <- file.path(sampleDirPath,basename)
      # grab current chr from data
      currentData <- subset(x,samplelist = currentSample,chromlist = currentChr)
      # render graphics for whole genome
      CairoPDF(file=chrFile,width,height)
      plot(currentData,altcol=F,ylim=ylim,...)
      dev.off()
    }
  }
}
