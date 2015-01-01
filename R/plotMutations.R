#' Plot mutations
#'
#' Plots mutations given mutation data, and matched
#' @param x list of mutation positions
#' @param y list of mutation counts for each corresponding position mutated
#' @param protLen protein length; can be from retrieveProtLen
#' @param annotatePos positions to annotate, default empty list
#' @param annotateSymbol symbol to use for annotation, default *
#' @param annotateProt protein annotations, e.g. conserved domains; can be from retrieveDomains
#' @return plots
#' @export
#' 
plotMutations <- function(x, y, 
                          UniProtID = NA,
                          protLen = NA,
                          annotatePos = c(), annotateSymbol = "*",
                          annotateProt = c()){
  
  #default colors
  mutcol <- rgb(180,0,14,220,maxColorValue=255)
  mutcolborder <- rgb(100,0,5,255,maxColorValue=255)
  stemcol <- rgb(150,150,150,255,maxColorValue=255)
  barcol <- rgb(100,100,100,255,maxColorValue=255)
  baranotcol <- rgb(0,130,5,255,maxColorValue=255)
  
  #initialize protLen
  #quality controls
  if(length(protLen) > 1){
    warning("protLen should not be a vector, will only take the first element")
    protLen <- protLen[1]
  }
  
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
  
  #initialize plot
  plot(0,0, xlim=c(0,protLen), ylim=c(-0.4,max(y)+1), ylab="Number of Mutations", xlab="Amino acid", 
       cex=0,frame.plot=FALSE)
  
  #plot mutations
  for(idx in 1:length(x)){
    segments(x[idx],0,x[idx],y[idx], col=stemcol) #plot stems
  }
  par(new=TRUE)
  plot(x, y, xlim=c(0,protLen), ylim=c(-0.4,max(y)+1), ylab="Number of Mutations", xlab="Amino acid", 
       cex=2, bg=mutcol, col=mutcolborder, pch=21, frame.plot=FALSE)
  
  #annotate positions
  if(!is.null(annotatePos)){
    y.max <- max(y)
    y.increment <- y.max*0.15
    for(idx in 1:length(annotatePos)){
      if(annotatePos[idx] %in% x){
        y.indx <- match(annotatePos[idx], x)
        text(annotatePos[idx]-0.1, y[y.indx]+y.increment, labels=annotateSymbol, cex=2)
      }
    }
  }
  
  #plot protein
  rect(0,-0.4,protLen,0,col=barcol,border=barcol)
  
  #annotate protein
  if(!is.null(annotateProt)){
    for(idx in 1:nrow(annotateProt)){
      posMiddle <- floor((as.numeric(annotateProt$end[idx])-as.numeric(annotateProt$start[idx]))/2)
      posMiddle <- as.numeric(annotateProt$start[idx]) + posMiddle
      rect(annotateProt$start[idx],-0.5,
           annotateProt$end[idx],0.1,
           col=baranotcol,border=baranotcol)
      text(posMiddle, -0.20, adj=0.5,
           labels=annotateProt$pfamName[idx], cex=0.75)
    }
  }
}
