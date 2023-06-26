#' A Go_deseq2_heat
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Heatmap
#' @export
#' @examples
#' Go_deseq2_heat()

Go_mergeTab <- function(pattern, file_path){
  
   # add input files
  path <- file_path
 
  
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  filenames <- list.files(path, pattern=pattern);filenames
  
  
  cat(sprintf("Files location: %s\n",path))
  cat("=======================================================================\n")
  cat("Merged files:\n")
  cat(sprintf("%s\n",filenames))

  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  return(df)
}
