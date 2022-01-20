#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
#dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
dircolors <- c("#4f86f7", "#e10000", "grey"); names(dircolors) <- c("down", "up", "NS")
Go_lmem_fore <- function(project,file_path, alpha, pattern, name, order, height, width){
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
  print(filenames)

  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/9_lmem.forest.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/9_lmem.forest.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }

  
  
  for (sn in sample.names) {
    file <- list.files(path, pattern = sprintf("%s%s", sn, sprintf("%s", pattern)), full.names =T)
    df <- read.csv(file,row.names=NULL ,check.names=FALSE)
    df.sel <- df
    resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$Estimate),]
    # resSig$smvar <- factor(resSig$smvar)
    print(1)
    if (dim(resSig)[1] == 0)
      next

    print(2)
    resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)==1, "up", "down"), "NS")
    print(3)
    
    # 중복 이름 처리 하기
    headers <- vector(dim(resSig)[1], mode="character")
    
    for (i in 1:dim(resSig)[1]) {
      headers[i] <- paste("ASV", i, sep="_")
    }
    resSig$ASV <- headers
    
    
    for (plot in unique(resSig$metadata)){
      resSig.sel <- subset(resSig, metadata == plot)
      if (length(unique(resSig.sel$coefficient)) ==1 ){
        
        lims <- max(abs(resSig.sel$Estimate) + abs(resSig.sel$SE))*1.0
        p1 <- ggplot(resSig.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
          geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
          geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
          scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel$ASV), labels = resSig.sel$taxa)
        
        p1 <- p1+ ggtitle(sprintf("LMEM-%s (%s p < %s) ",plot, sn, alpha)) +
          labs(y = "Estimate") +labs(x = NULL)
        
        print(p1)
      } else if (length(unique(resSig.sel$coefficient)) >=1){
        for (cof in unique(resSig.sel$coefficient)){
          resSig.sel.sel <- subset(resSig.sel, coefficient == cof)
          lims <- max(abs(resSig.sel.sel$Estimate) + abs(resSig.sel.sel$SE))*1.0
          p1 <- ggplot(resSig.sel.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
            geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
            geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
            scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel.sel$ASV), labels = resSig.sel.sel$taxa)
          
          p1 <- p1+ ggtitle(sprintf("LMEM-%s (%s p < %s) ",cof, sn, alpha)) +
            labs(y = "Estimate") +labs(x = NULL)
          print(p1)
        }
      }
    }
  }
  dev.off()
}





