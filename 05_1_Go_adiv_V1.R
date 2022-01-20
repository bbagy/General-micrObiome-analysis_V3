#' A Go_adiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Alpha diversity plot
#' @export
#' @examples
#' Go_adiv


Go_adiv <- function(psIN, project, alpha_metrics){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_adiv <- file.path(sprintf("%s_%s/table/adiv",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_adiv)) dir.create(out_adiv)
  
  
  # adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics) # se.chao1 stand error
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  cat(sprintf("adiv table is saved in %s.\n",out_path))
  cat("                                                       \n")
  write.csv(adiv, quote = FALSE, col.names = NA, 
            file=sprintf("%s/adiv.%s.%s.csv",out_adiv,project, format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
  return(adiv)
} 
  
