#' A Go_deseq2_volc
#'


Go_DA_volc <- function(project, file_path,files, type,alpha,beta, name,font, height, width){
    
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
  pdf(sprintf("%s/DA.volcano.%s.%s(%s.%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              alpha, 
              beta, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  


  for (fn in 1:length(filenames)) {

    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")
    # remove NA
    df[df==""] <- "NA"
    df$deseq2 <- ifelse(df$padj < alpha & abs(df$log2FoldChange) > beta, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$dir), ]

    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)


    # colors and names
    df.na$deseq2<- gsub('down',basline, gsub('up',smvar, df.na$deseq2))
    
    df.na$deseq2 <- factor(df.na$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
    
    dircolors <- c("#1170aa", "grey","#fc7d0b"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))

    df.na$ancom[is.na(df.na$ancom)] <- "NS"
    ancomshape <- c(18,5,5); names(ancomshape) <- c(TRUE, FALSE, "NS")# 16,1,1
    
    
    p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue),colour=deseq2)) + theme_bw() +
      scale_color_manual(values=dircolors) + 
      xlab("log2 fold change") + ylab("-log10 (p-value)")+ 
      geom_vline(xintercept = -beta,col = "#1170aa", linetype = "dotted", size = 1) + 
      geom_vline(xintercept = beta,col = "#fc7d0b", linetype = "dotted", size = 1) + 
      theme(text = element_text(size=font+8),plot.title = element_text(size=font+8), legend.text=element_text(size=font+8), 
            legend.position="bottom", legend.title = element_blank())# + theme()

    # ancom
    p1 = p1 + geom_point(aes(shape=ancom), alpha=1, size=font-1.5) + scale_shape_manual(values = ancomshape) + guides(shape = FALSE)
    


    
    if(type == "taxonomy" | type == "taxanomy" |type == "bacmet" ){
      p1 <- p1 +  geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < alpha & abs(df.na$log2FoldChange) > beta, as.character(ShortName),'')), size=font, segment.alpha = 0.25, fontface="italic")
    }else if(type == "function"){
      p1 <- p1 +  geom_text_repel(aes(label=ifelse(KOName != "NA" & df.na$padj < alpha & abs(df.na$log2FoldChange) > beta, as.character(KOName),'')), size=font)
    }
    
    if(!is.null(df.na$des)){
      des <- unique(df.na$des)
      p1 <- p1 + ggtitle(sprintf("%s-%s, %s vs %s (padj < %s,cutoff=%s) ", mvar, des,basline, smvar,  alpha, beta))
    }else{
      p1 <- p1 + ggtitle(sprintf("%s, %s vs %s (padj < %s,cutoff=%s) ", mvar, basline, smvar,  alpha, beta))
    }

    print(p1)
  } 
  dev.off()
}

