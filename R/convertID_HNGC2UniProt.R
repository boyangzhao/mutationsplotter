#' Convert HGNC gene symbol to UniProt Accession IDs
#'
#' Convert HGNC gene symbol to UniProt Accession IDs using BioMart
#' @param HGNCSymbol HUGO gene symbol
#' @return UniProt Accession ID
#' @export
#' 
convertID_HGNC2UniProt <- function(HGNCSymbol){
  library(biomaRt)
  
  #listMarts()
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  #listFilters(mart)
  #listAttributes(mart)
  results <- getBM(attributes = c("uniprot_swissprot_accession"), filters="hgnc_symbol", values=HGNCSymbol, mart=mart)
  uniprotID <- unique(results$uniprot_swissprot_accession)
  
  if(length(uniprotID) == 1){
    return(uniprotID)
  } else {
    warning(gettextf("Multiple IDs found for given gene %s", HGNCSymbol))
    return(NULL)
  }
}
