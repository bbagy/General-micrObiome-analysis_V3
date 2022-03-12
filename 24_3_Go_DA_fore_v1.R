#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
#dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")


Go_DA_fore <- function(project,file_path, files,type, alpha, beta,font, name, height, width){
    
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
  pdf(sprintf("%s/DA.forest.%s.%s(%s.%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              alpha, 
              beta, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  


  for (fn in 1:length(filenames)) {
    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")

    
    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)
    
    df$deseq2 <- ifelse(df$padj < alpha, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.sel <- df
    resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > beta))
    if (dim(resSig)[1] == 0 | dim(resSig.top)[1] == 0 ){
      next
    }
    
    resSig$smvar <- factor(resSig$smvar)
    lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
    
    # save top for fishtaco
    if(type == "function"){
      top.ko <- resSig.top$KO
      #out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))); if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
      #write.csv(top.ko, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.VS.%s.%s.%s.%s.csv",out_deseq2,mvar, basline, smvar, project,format(Sys.Date(), "%y%m%d"),"Forfishtaco",sep="/"))
    }else{
    }
    
    # colors and names
    resSig.top$deseq2<- gsub('down',basline, gsub('up',smvar, resSig.top$deseq2))
    
    resSig.top$deseq2 <- factor(resSig.top$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
    
    dircolors <- c("#1170aa", "grey","#fc7d0b"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    
    
    resSig.top$ancom[is.na(resSig.top$ancom)] <- FALSE
    ancomshape <- c(18,5); names(ancomshape) <- c(TRUE, FALSE) # 16,1,1
    
    #dircolors <- c("#f7022a", "#4f86f7","grey"); names(dircolors) <- c("down", "up", "NS")
    
    p1 <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=log2FoldChange, color=deseq2)) + 
      geom_hline(yintercept=0) + geom_point(aes(shape=ancom)) + coord_flip() + theme_classic() + 
      scale_color_manual(values=dircolors) + scale_shape_manual(values = ancomshape) + #guides(shape = "none") +
       #theme_classic() +theme_bw() 
      geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2)  + 
      ylim(c(-lims, lims))+ xlab("Taxa") + ylab("log2FoldChange")+
      theme(text = element_text(size=font), plot.title = element_text(size=font, hjust = 1),
            axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font,face = "italic")) #hjust =1
 

    if(type == "taxonomy" | type == "taxanomy"){
      p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s__%s__%s",as.character(resSig$Phylum),as.character(resSig$Family), as.character(resSig$ShortName))) 
    }else if(type == "function"){
      p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s__%s",as.character(resSig$Path.des),as.character(resSig$KO.des))) 
    }else if(type == "bacmet" ){
      p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s",as.character(resSig$ShortName))) 
    }
    
    
    if(!is.null(resSig.top$des)){
      des <- unique(resSig.top$des)
      p1 <- p1 + ggtitle(sprintf("%s-%s, %s vs %s (padj < %s,cutoff=%s) ", mvar,des, basline, smvar,  alpha, beta)) 
      
    }else{
      p1 <- p1 + ggtitle(sprintf("%s, %s vs %s (padj < %s,cutoff=%s) ", mvar, basline, smvar,  alpha, beta)) 
      
    }
        

    print(p1)
  } 
  dev.off()
}

