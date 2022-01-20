#' A Go_deseq2_volc
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Volcano plot
#' @export
#' @examples
#' Go_deseq2_volc()

Go_deseq2_volc <- function(project, file_path,files, type,alpha,beta, name, height, width){
    
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
  if (!is.null(name)) {
    pdf(sprintf("%s_%s/pdf/6_deseq2.volcano.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }   else {
      pdf(sprintf("%s_%s/pdf/6_deseq2.volcano.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, beta, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }


  for (fn in 1:length(filenames)) {

    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")
    # remove NA
    df[df==""] <- "NA"
    df$dir <- ifelse(df$padj < alpha & abs(df$log2FoldChange) > beta, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$dir), ]

    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)


    # colors and names
    df.na$dir<- gsub('down',basline, gsub('up',smvar, df.na$dir))
    
    df.na$dir <- factor(df.na$dir, levels = c(as.character(basline), "NS", as.character(smvar)))
    
    dircolors <- c("#f7022a", "grey","#4f86f7"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))

    
    p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue), colour=dir)) +theme_bw() +
      geom_point(alpha=1, size=2) + scale_color_manual(values=dircolors) +xlab("log2 fold change") + ylab("-log10 (p-value)")+ geom_vline(xintercept = -beta,col = "blue", linetype = "dotted", size = 1) + geom_vline(xintercept = beta,col = "red", linetype = "dotted", size = 1) + theme(plot.title = element_text(size=8), legend.position="bottom", legend.title = element_blank())# + theme()
    

    
    if(type == "taxonomy" | type == "taxanomy" |type == "bacmet" ){
      p1 <- p1 +  geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < alpha & abs(df.na$log2FoldChange) > beta, as.character(ShortName),'')), size=2.5, segment.alpha = 0.25, fontface="italic")
    }else if(type == "function"){
      p1 <- p1 +  geom_text_repel(aes(label=ifelse(KOName != "NA" & df.na$padj < alpha & abs(df.na$log2FoldChange) > beta, as.character(KOName),'')), size=2.5)
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

