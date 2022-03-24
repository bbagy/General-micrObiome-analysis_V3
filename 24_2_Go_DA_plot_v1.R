#' A Go_deseq2_volc
#'


Go_DA_plot <- function(project, file_path,files, type, plot = "volcano", fdr, fc,mycols=NULL, name, overlaps=10, font, height, width){
    
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
  
  pdf(sprintf("%s/DA.%s%s.%s(%s.%s).%s.pdf", out_path, 
              ifelse(is.null(plot), "", paste(plot, ".", sep = "")), 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr, 
              fc, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  


  for (fn in 1:length(filenames)) {

    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")
    # remove NA
    df[df==""] <- "NA"
    df$deseq2 <- ifelse(df$padj < fdr & abs(df$log2FoldChange) > fc, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$dir), ]

    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)


    # colors and names
    df.na$deseq2<- gsub('down',basline, gsub('up',smvar, df.na$deseq2))
    
    df.na$deseq2 <- factor(df.na$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
    

   if(!is.null(mycols)){
    dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }else{
     dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }
    
    
    

    df.na$ancom[is.na(df.na$ancom)] <- FALSE
    ancomshape <- c(18,5); names(ancomshape) <- c(TRUE, FALSE)# 16,1,1
    
    #------------------#
    #   plot style     #
    #------------------#
    if(plot == "volcano"){
      print("Generating Volcano plots.")
      p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue),colour=deseq2)) + 
        xlab("log2 fold change") + ylab("-log10 (p-value)")+ 
        geom_vline(xintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1) 
      
    } else if(plot == "maplot"){
      print("Generating M (log ratio) A (mean average)  plots.")
      p1 <-  ggplot(df.na, aes(x=log2(baseMean +1), y=log2FoldChange, colour=deseq2)) +
        xlab("Log2 mean expression") + ylab("Log2 fold change")+ 
        geom_hline(yintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1)
      
    } else if(plot == "forest"){
      print("Generating forest plots.")
      resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
      resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > beta))
      if (dim(resSig)[1] == 0 | dim(resSig.top)[1] == 0 ){
        next
      }
      
      resSig$smvar <- factor(resSig$smvar)
      lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
      resSig.top$deseq2<- gsub('down',basline, gsub('up',smvar, resSig.top$deseq2))
      resSig.top$deseq2 <- factor(resSig.top$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
      resSig.top$ancom[is.na(resSig.top$ancom)] <- FALSE

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
      } else if(type == "bacmet" ){
        p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s",as.character(resSig$ShortName))) 
      }
    }
    
    
    if(plot == "volcano" |  plot == "maplot"){
      p1 <- p1 + theme_bw() + scale_color_manual(values=dircolors) + theme(text = element_text(size=font+8),plot.title = element_text(size=font+8), legend.text=element_text(size=font+8),  legend.position="bottom",legend.justification = "left",legend.box = "vertical")  +
        geom_point(aes(shape=ancom), size=font-1.5) + scale_shape_manual(values = ancomshape)
      if(type == "taxonomy" | type == "taxanomy" |type == "bacmet" ){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(ShortName),'')), size=font, segment.fdr = 0.25, fontface="italic",max.overlaps = overlaps )
      }else if(type == "function"){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(KOName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(KOName),'')), size=font,max.overlaps = overlaps)
      }
    }

    
    if(!is.null(df.na$des)){
      des <- unique(df.na$des)
      p1 <- p1 + ggtitle(sprintf("%s-%s, %s vs %s (padj < %s,cutoff=%s) ", mvar, des,basline, smvar,  fdr, fc))
    }else{
      p1 <- p1 + ggtitle(sprintf("%s, %s vs %s (padj < %s,cutoff=%s) ", mvar, basline, smvar,  fdr, fc))
    }
    print(p1)
  } 
  dev.off()
}

