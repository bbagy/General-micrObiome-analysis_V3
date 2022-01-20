#' A Go_bdiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Go_bdiv()

Go_dist <- function(psIN, project, distance_metrics){
  # out dir
  #out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out)) dir.create(out)
  #out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_path)) dir.create(out_path)
  #out_perm <- file.path(sprintf("%s_%s/table/perm ",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  # run distance
  dm <- list()
  for (distance_metric in distance_metrics) {
    dm[[length(dm)+1]] <- phyloseq::distance(psIN, method=distance_metric)
  }
  
  names(dm) <- distance_metrics
  class(dm)
  
  return(dm)
}
  