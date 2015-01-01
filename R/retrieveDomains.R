#' Retrieve conserved domains
#'
#' @param UniProtID UniProt Accession ID
#' @return conserved domains
#' @export
#' 
retrieveDomains <- function(UniProtID){
  library(XML)
  
  #retrive data from pfam database
  xmlout <- xmlParse(gettextf("http://pfam.xfam.org/protein/%s?output=xml", UniProtID))
  xmldata <- xmlToList(xmlout)
  
  if(length(xmldata) > 1){
    domains <- xmldata$entry$matches
    domainsN <- length(domains)
    
    #extract relevant data
    start <- rep(NA, domainsN)
    end <- rep(NA, domainsN)
    pfamID <- rep(NA, domainsN)
    pfamName <- rep(NA, domainsN)
    for(idx in 1:length(domains)){
      start[idx] <- domains[idx]$match$location['ali_start']
      end[idx] <- domains[idx]$match$location['ali_end']
      pfamID[idx] <- domains[idx]$match$.attrs['accession']
      pfamName[idx] <- domains[idx]$match$.attrs['id']
    }
    
    domainData <- data.frame(start=start,
                             end=end,
                             pfamID=pfamID,
                             pfamName=pfamName,
                             stringsAsFactors=FALSE)
    return(domainData)
  } else {
    return(NULL)
  }

}
