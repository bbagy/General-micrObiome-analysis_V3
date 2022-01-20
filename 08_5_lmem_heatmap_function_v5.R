#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_lmem_heat <- function(project,file_path, alpha, pattern, facet, name, orders, height, width){
  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  
  
  print(path)
  print(sample.names)
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1 & length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  
  df.sel <- df
  resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$Estimate),]
  # resSig$smvar <- factor(resSig$smvar)
  
  
  resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)==1, "up", "down"), "NS")
  print(1)
  for (plot in unique(resSig$metadata)){
    resSig.sel <- subset(resSig, metadata == plot)
    print(2)
    if (length(unique(resSig.sel$coefficient)) >=1 ){
      resSig.sel$comparison <-  gsub(sprintf("%s", unique(resSig.sel$metadata)) ,"" , resSig.sel$coefficient)
      resSig.sel$comparison <- factor(resSig.sel$comparison , levels = orders)
      resSig.sel$stars <- cut(resSig.sel$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
      print(3)
      if (length(facet) == 1) {
        resSig.sel[,facet] <- factor(resSig.sel[,facet] , levels = orders)
      }
      print(4)
      p <- ggplot(resSig.sel, aes(x=reorder(taxa,Estimate), y=comparison, color=comparison)) + #, alpha=padj)) +
        theme_classic()+ coord_flip() + geom_tile(aes(fill = Estimate), colour = "white") + 
        geom_text(aes(label=stars), color = "black") + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
        ggtitle(sprintf("LMEM All comparison group (%s p < %s) ", plot, alpha)) +
        theme(plot.title = element_text(hjust = 0.5))+ #0.5
        theme(legend.position= "right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + labs(y = "comparison group") +labs(x = NULL)
      p = p + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
      
      
      if (length(facet) == 1) {
        print(5)
        ncol <- length(unique(resSig.sel[,facet]))*length(unique(resSig.sel[,"comparison"]))
        p = p + facet_wrap(as.formula(sprintf("~ %s+%s", facet, "comparison")), scales="free_x", ncol = ncol)
      }
      else {
        p = p + facet_wrap(~  comparison, scales="free_x", ncol = 10) 
      }
      
      #plotlist[[length(plotlist)+1]] <- p
      print(p)
    }
  }  
  dev.off()
}
