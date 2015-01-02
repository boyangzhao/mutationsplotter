#' Convert HGNC gene symbol to UniProt Accession IDs
#'
#' Convert HGNC gene symbol to UniProt Accession IDs using BioMart
#' @param HGNCSymbol HUGO gene symbol
#' @return UniProt Accession ID
#' @export
#' 
convertID_UniProt2HGNC <- function(UniProtID){
  library(biomaRt)
  
  #listMarts()
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  #listFilters(mart)
  #listAttributes(mart)
  results <- getBM(attributes = c("hgnc_symbol"), filters="uniprot_swissprot_accession", values=UniProtID, mart=mart)
  geneSymbol <- unique(results$hgnc_symbol)
  
  if(length(geneSymbol) == 1){
    return(geneSymbol)
  } else {
    warning(gettextf("Multiple HGNC symbols found for given UniProt ID %s", UniProtID))
    return(NULL)
  }
}