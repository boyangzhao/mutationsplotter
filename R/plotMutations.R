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
                          annotateProt = NA){
  
  #default colors
  mutcol <- rgb(180,0,14,220,maxColorValue=255)
  mutcolborder <- rgb(100,0,5,255,maxColorValue=255)
  stemcol <- rgb(150,150,150,255,maxColorValue=255)
  barcol <- rgb(100,100,100,255,maxColorValue=255)
  
  #initialize protLen
  if(is.NA(protLen)){
    if(!is.NA(UniProtID)){
      #retrieve protein length if gene name (UniProtID) is provided and protein length undefined
      protLen <- retrieveProtLen(UniProtID)
    } else {
      #otherwise default to maximum value of x
      protLen <- max(x)
    }
  }
  
  #initialize annotateProt
  if(is.na(annotateProt) && !is.na(UniProtID)){
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
  y.max <- max(y)
  y.increment <- y.max*0.15
  for(idx in 1:length(annotatePos)){
    if(annotatePos[idx] %in% x){
      y.indx <- match(annotatePos[idx], x)
      text(annotatePos[idx]-0.1, y[y.indx]+y.increment, labels=annotateSymbol, cex=2)
    }
  }
  
  #plot protein
  rect(0,-0.4,protLen,0,col=barcol,border="gray")
}
