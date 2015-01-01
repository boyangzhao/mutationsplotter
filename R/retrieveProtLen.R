#' Retrieve protein lengths
#'
#' Retrieve protein lengths using UniProt web service
#' @param UniProtID UniProt Accession ID
#' @return protein length
#' @export
#' 
retrieveProtLen <- function(UniProtID){
  library(RCurl)
  fasta.data <- getURL(gettextf("http://www.uniprot.org/uniprot/%s.fasta", UniProtID))
  
  #gets rid of header
  idx <- regexpr("\\\n",fasta.data)
  idx <- idx[1]
  
  #remove \n's
  seq <- substr(fasta.data,idx[1]+1,nchar(fasta.data))
  seq <- gsub("\\\n","",seq)
  
  #calculate protein length
  protLen <- nchar(seq)
  
  return(protLen)
}
