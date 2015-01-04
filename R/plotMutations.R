#' Plot mutations
#'
#' Plots mutations given mutation data, and matched
#' @param x list of mutation positions
#' @param y list of mutation counts for each corresponding position mutated
#' @param protLen protein length; can be from retrieveProtLen
#' @param annotatePos positions to annotate, default empty list
#' @param annotateSymbol symbol to use for annotation, default *
#' @param annotateProt protein annotations, e.g. conserved domains; can be from retrieveDomains
#' @param xaxis Boolean to turn on/off x-axis
#' @param yaxis Boolean to turn on/off y-axis
#' @param xlab x-axis label, default = "Amino acid", only show if xaxis = TRUE
#' @param ylab y-axis label, default = "Number of Mutations", only show if yaxis = TRUE
#' @param main title, value '[gene]' will be replaced by actual gene name (queried using given UniProtID)
#' @return plots
#' @export
#' 
plotMutations <- function(x, y, 
                          UniProtID = NA,
                          protLen = NA,
                          annotatePos = c(), annotateSymbol = "*",
                          annotateProt = c(),
                          xaxis = TRUE, yaxis = TRUE,
                          xlab = "Amino acid", ylab = "Number of Mutations",
                          main = "[gene]"){
  
  #default colors
  mutcol <- rgb(180,0,14,220,maxColorValue=255)
  mutcolborder <- rgb(100,0,5,255,maxColorValue=255)
  stemcol <- rgb(150,150,150,255,maxColorValue=255)
  barcol <- rgb(100,100,100,255,maxColorValue=255)
  baranotcol <- rgb(0,130,5,255,maxColorValue=255)
  
  #quality controls
  if(length(protLen) > 1){
    warning("protLen should not be a vector, will only take the first element")
    protLen <- protLen[1]
  }
  if(nchar(UniProtID) < 1){
    UniProtID <- NA
  }
  
  #initialize protLen
  if(is.na(protLen)){
    if(!is.na(UniProtID)){
      #retrieve protein length if gene name (UniProtID) is provided and protein length undefined
      protLen <- retrieveProtLen(UniProtID)
      
      if(length(protLen) < 1 || protLen < 1){
        warning("Error in retrieving protein length information, defaulting to protein length as maximum of given x value")
        protLen <- max(x)
      }
    } else {
      #otherwise default to maximum value of x
      protLen <- max(x)
    }
  }
  
  #initialize annotateProt
  if(length(annotateProt)<1 && !is.na(UniProtID)){
    annotateProt <- retrieveDomains(UniProtID)
  }
  
  #derive additional parameter values
  geneSymbol <- ""
  titleName <- ""
  if(!is.na(UniProtID)){
    geneSymbol <- convertID_UniProt2HGNC(UniProtID)
    titleName <- gsub("\\[gene\\]", geneSymbol, main)
  }
  
  #derive plot parameters
  ymax <- ceiling((max(y)/10))*10
  ymin <- 0 #lowest point of y, with tickmarks
  protheight <- ymax*0.09
  protheightAdd <- (protheight*0.18)
  #protheightAdd <- 0
  ylower <- -(protheight+protheightAdd) #lower point of y, without tickmarks
  
  #initialize plot
  plot(0,0, xlim=c(0,protLen), axes=FALSE, ylim=c(ylower,ymax), xlab="", ylab="",
       cex=0, frame.plot=FALSE, main=titleName)
  
  #plot mutations
  for(idx in 1:length(x)){
    segments(x[idx],0,x[idx],y[idx], col=stemcol) #plot stems
  }
  par(new=TRUE)
  plot(x, y, xlim=c(0,protLen), axes=FALSE, ylim=c(ylower,ymax), xlab="", ylab="",
       cex=2, bg=mutcol, col=mutcolborder, pch=21, frame.plot=FALSE)
  
  if(xaxis){
    axis(1)
    title(xlab=xlab)
  }
  
  if(yaxis){
    axis(2, at=seq(ymin, ymax, by=floor((ymax-ymin)/5)))
    title(ylab=ylab)
  }
    
  #annotate positions
  if(!is.null(annotatePos)){
    y.increment <- ymax*0.15
    for(idx in 1:length(annotatePos)){
      if(annotatePos[idx] %in% x){
        y.indx <- match(annotatePos[idx], x)
        text(annotatePos[idx]-0.1, y[y.indx]+y.increment, labels=annotateSymbol, cex=2)
      }
    }
  }
  
  #plot protein
  rect(0,-protheight,protLen,0,col=barcol,border=barcol)
  
  #annotate protein
  if(!is.null(annotateProt)){
    for(idx in 1:nrow(annotateProt)){
      posMiddle <- floor((as.numeric(annotateProt$end[idx])-as.numeric(annotateProt$start[idx]))/2)
      posMiddle <- as.numeric(annotateProt$start[idx]) + posMiddle
      rect(annotateProt$start[idx],-(protheight+protheightAdd),
           annotateProt$end[idx], protheightAdd,
           col=baranotcol,border=baranotcol)
      text(posMiddle, -((protheight+protheightAdd)/2.1), adj=0.5,
           labels=annotateProt$pfamName[idx], cex=0.75)
    }
  }
}
