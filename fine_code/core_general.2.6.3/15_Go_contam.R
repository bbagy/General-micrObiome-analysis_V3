#' A Go_contam
#'
#'
#' @param Go_contam
#' @keywords Go_contam
#' @export
#' @examples
#' Go_contam


Go_contam <- function(physeq, contams, samples) {
  contams.n1 = prune_taxa(taxa_sums(contams)>=1, contams) 
  samples.n1 = prune_taxa(taxa_sums(samples)>=1, samples) 
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(contams.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samples.n1))
}