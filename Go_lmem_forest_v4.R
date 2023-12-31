#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_lmem_fore()
# dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
# dircolors <- c("#4f86f7", "#e10000", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_lmem_fore <- function(project,file_path, mycols=NULL, fdr=0.05,est=1, files, name, order, height, width){
    if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=files);filenames
  
  print(path)
  print(filenames)
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  
  # out file
  pdf(sprintf("%s/lmem.forest.%s.%s(fdr=%s.est=%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr,
              est,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  
  
  
  
  
  for (fn in 1:length(filenames)) {
    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")

    df.sel <- df
    
    df.sel.sel <- as.data.frame(subset(df.sel, abs(df.sel$Estimate) > est))
    
    resSig <- as.data.frame(subset(df.sel.sel, pvalue < fdr )); resSig <- resSig[order(resSig$Estimate),]
    # resSig$smvar <- factor(resSig$smvar)
    print(1)
    if (dim(resSig)[1] == 0)
      next

    print(2)
    resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)== 1, "up", "down"), "NS")
    
  
    
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
        
        baseline <- unique(resSig.sel$baseline)
        compare <- unique(resSig.sel$coefficient)
        compare <- gsub(unique(resSig.sel$metadata), "", compare);compare
        
        
        # color
        resSig.sel$dir <- gsub('down',baseline, gsub('up',compare, resSig.sel$dir))
        resSig.sel$dir <- factor(resSig.sel$dir, levels = c(as.character(baseline), "NS", as.character(compare)))
        
        # color end
        
        if(!is.null(mycols)){
          dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
        }else{
          dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
        }
        
        legend.labs <- 
          c(paste(baseline, " (n=", unique(resSig.sel$bas.count),")",sep=""),
            paste("NS"),
            paste(compare, " (n=", unique(resSig.sel$coef.count), ")",sep=""))
        
        
        
        
        lims <- max(abs(resSig.sel$Estimate) + abs(resSig.sel$SE))*1.0
        p1 <- ggplot(resSig.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
          geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
          geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
          scale_color_manual(values=dircolors,labels=legend.labs)  + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel$ASV), labels = resSig.sel$taxa)+theme(plot.title = element_text(size=8, hjust = .5))
        
        #sn <- gsub("\\..*","",sn);sn
        p1 <- p1+ ggtitle(sprintf("LMEM-%s ( pvalue < %s, cutoff=%s (colored fdr < 0.05)) ", plot, fdr, est)) + labs(y = "Estimate") +labs(x = NULL)
        
        print(p1)
      } else if (length(unique(resSig.sel$coefficient)) >=1){
        for (cof in unique(resSig.sel$coefficient)){
          resSig.sel.sel <- subset(resSig.sel, coefficient == cof)
          
          baseline <- unique(resSig.sel.sel$baseline)
          compare <- unique(resSig.sel.sel$coefficient);compare
          compare <- gsub(unique(resSig.sel.sel$metadata), "", compare);compare
          
          
          # color
          resSig.sel.sel$dir <- gsub('down',baseline, gsub('up',compare, resSig.sel.sel$dir))
          resSig.sel.sel$dir <- factor(resSig.sel.sel$dir, levels = c(as.character(baseline), "NS", as.character(compare)))
          
          # color end
          
          if(!is.null(mycols)){
            dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
          }else{
            dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
          }
          
          legend.labs <- 
            c(paste(baseline, " (n=", unique(resSig.sel.sel$bas.count),")",sep=""),
              paste("NS"),
              paste(compare, " (n=", unique(resSig.sel.sel$coef.count), ")",sep=""))
          
          
          
          lims <- max(abs(resSig.sel.sel$Estimate) + abs(resSig.sel.sel$SE))*1.0
          
          p1 <- ggplot(resSig.sel.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
            geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
            geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
            scale_color_manual(values=dircolors, labels=legend.labs) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel.sel$ASV), labels = resSig.sel.sel$taxa)+theme(plot.title = element_text(size=8, hjust = 1))
          
          #sn <- gsub("\\..*","",sn);sn
          
          p1 <- p1+ ggtitle(sprintf("LMEM  (in %s for %s pvalue < %s, cutoff=%s (colored fdr < 0.05)) ",
                                    unique(resSig.sel.sel$metadata), compare, fdr,est)) +
            labs(y = "Estimate") +labs(x = NULL)
          print(p1)
        }
      }
    }
  }
  dev.off()
}





