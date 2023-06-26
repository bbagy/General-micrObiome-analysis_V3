#' A Go_qq
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords qq plot
#' @export
#' @examples
#' Go_qq()


Go_qq <- function(psIN, project, alpha_metrics, name, height, width){
    if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
      # logic for out file
  pdf(sprintf("%s/QQplot.%s%s.%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # 1st adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  adiv$ShannonLn <-log(adiv$Shannon)
  # show last column name
  rev(names(adiv))[1]
  
  #----------- QQ plot and histogram -----------#
  par(mfrow = c(3,2))
  mes <- c(alpha_metrics, rev(names(adiv))[1])
  for (am in mes){
    test <- shapiro.test(adiv[,am])
    hist(adiv[,am], freq=F, xlab= am, main=sprintf("Histogram of %s (%s)", project, am ), cex.main=1) 
    lines(density(adiv[,am])) 
    rug(adiv[,am])
    # remove inf
    adiv.inf <- adiv[!is.infinite(adiv[,am]),]
    
    qqnorm(adiv.inf[,am], main=sprintf("Normal Q-Q Plot (%s p=%.2g)", "shapiro", test$p.value), cex.main=1)
    qqline(adiv.inf[,am])
    
    print(sprintf("%s %s shapiro test (p=%.2g)",am, project, test$p.value))
  }
  
  dev.off()
}




