dna.clustering <- function(x,...){
  source("~/scripts/heatmap.R")
}


cghToCrosstab <- function(x,factor){
  require(CGHregions)
  f <- factor(pData(x)[,factor])
  ct <- list()
  for(i in featureNames(x)){
    ct[[i]] <- data.frame(row.names=levels(f))
    for(group in levels(f)){
      ct[[i]][group,"loss"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]<0)
      ct[[i]][group,"normal"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]==0)
      ct[[i]][group,"gain"] <- sum(regions(x)[i,which(pData(x)[,factor]==group)]>0)
    }
  }
  return(ct)
}

cghFreqPlot <-function(x,fill=T){
  require(CGHregions)
  require(reshape)
  require(ggplot2)
  totalSamples <- ncol(x)
  r <- regions(x)
  loss <- apply(r,1,function(x){sum(x<0)}) / totalSamples
  gain <- apply(r,1,function(x){sum(x>0)}) / totalSamples
  df <- data.frame(
    chr = chromosomes(x),
    midpoint = rowMeans(cbind(bpend(x),bpstart(x))),
    size = bpend(x) - bpstart(x),
    loss = loss * -1,
    gain = gain
  )
  df <- melt(df,id.vars=c("chr","midpoint","size"))
  colnames(df)[4:5] <- c("direction","percentage")
  p <- ggplot(df) + 
    geom_step(aes(x=midpoint,y=percentage,color=direction)) + 
    facet_grid(.~chr,scales="free_x",space="free_x") + 
    scale_color_manual(values=c("blue","red")) +
    theme_bw() +
    ylim(c(-1,1)) + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="top",
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) + 
    labs(
      x="Genomic position",
      y=paste("Percentage of cohort (N=",totalSamples,")",sep="")
    )
  if(fill){
    p <- p +  
      geom_bar(
        aes(x=midpoint,y=percentage,fill=direction,width=size),
        alpha=1/3,
        stat="identity",
        position="identity"
      ) +
      scale_fill_manual(values=c("blue","red"))
  }
  return(p)
} 

bedFreqPlot <- function(x,fill=T){
  require(ggplot2)
  colnames(x) <- c("chr","start","end","group","value")
  x$midpoint <- rowMeans(x[,c("start","end")])
  x$size <- x$end - x$start
  p <- ggplot(x) + 
    geom_step(aes(x=midpoint,y=value,color=group)) + 
    facet_grid(.~chr,scales="free_x",space="free_x") + 
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position="top",
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) + 
    labs(
      x="Genomic position",
      y=paste("Frequency")
    )
  if(fill){
    p <- p +  
      geom_bar(
        aes(x=midpoint,y=value,fill=group,width=size),
        alpha=1/3,
        stat="identity",
        position="identity"
      )
  }
  return(p)
}