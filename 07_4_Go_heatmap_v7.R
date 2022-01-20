#' A Go_deseq2_heat
#' 
Go_deseq2_heat <- function(df, project, data_type, facet,groupby,font, alpha,beta, orders, name, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  

  
  
  # out file
  if (!is.null(name)) {
    pdf(sprintf("%s_%s/pdf/deseq2.heatmap.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (!is.null(facet)) {
    pdf(sprintf("%s_%s/pdf/deseq2.heatmap.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (!is.null(facet) & !is.null(name)) {
    pdf(sprintf("%s_%s/pdf/deseq2.heatmap.%s.%s.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, name, alpha,beta,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/deseq2.heatmap.%s.(%s.%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha,beta, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  resSig <- as.data.frame(subset(df, padj < alpha)); resSig <- resSig[order(resSig$log2FoldChange),]
  resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > beta))
  #print("c")
  #if (length(unique(resSig$smvar)) >=2 ){
  
  if (dim(resSig)[1] >= 1) {
    # re-order
    if (length(orders) >= 1) {
      if (groupby == "smvar"){
        resSig.top$smvar <- factor(resSig.top$smvar, levels = orders)
        #resSig.top$des <- factor(resSig.top$des, levels = orders)
        if (length(unique(resSig.top$smvar)) <= 1) 
          next
        resSig.top$smvar <- factor(resSig.top$smvar, levels = orders)
        #resSig.top$des <- factor(resSig.top$des, levels = orders)
      } else{
        resSig.top$des <- factor(resSig.top$des, levels = orders)
        if (length(unique(resSig.top$smvar)) > 1){
          resSig.top$smvar <- factor(resSig.top$smvar, levels = orders)
        }
      }
    } else {
      if (groupby == "smvar"){
        if (length(unique(resSig.top$smvar)) <= 1) 
          next
        resSig.top$smvar <- factor(resSig.top$smvar)
      } else{
        if (length(unique(resSig.top$des)) <= 1) 
          next
        resSig.top$des <- factor(resSig.top$des)
      }
    }
    
    print(1)
    if (groupby == "smvar"){
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=smvar, color=smvar)) + theme_classic()+ coord_flip() #x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
 
    }  else {
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=smvar, color=smvar)) + theme_classic()+ coord_flip()#x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
    }
    
    
    p = p + geom_tile(aes(fill = log2FoldChange), colour = "white") + 
      labs(y = "Comparison Group") +labs(x = NULL) +
      scale_fill_gradient2(low = "#1170aa", mid = "white", high = "#fc7d0b")+
      ggtitle(sprintf("%s baseline %s vs %s (padj < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig$basline), "All groups",  alpha,beta))  + 
      theme(plot.title = element_text(hjust = 0.5),legend.position= "right")+ #0.5
      theme(axis.text.x = element_text(angle=0, vjust=0.5, hjust=1, size=8),
             axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=8,face = "italic")) 
    
    
    print(2)
    if (data_type == "dada2" | data_type == "DADA2") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$Phylum, resSig$ShortName)))
    } else if (data_type == "Other" | data_type == "other") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$KOName)))
    }
    
    print(3)
    if (groupby == "smvar"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"smvar"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "smvar", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  smvar, scales="free_x", ncol = 10)
      }
    }else if (groupby == "des"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"des"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "des", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  des, scales="free_x", ncol = 10)
      }
    }
    #print(4)
    #plotlist[[length(plotlist)+1]] <- p
    p3 = p2 + theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + theme(text = element_text(size=font), plot.title = element_text(hjust=1))
   # print(p3)
  }else{
    next
  }

  
  p4 <- ggplotGrob(p3)
  id <- which(p4$layout$name == "title")
  p4$layout[id, c("l","r")] <- c(1, ncol(p4))
  #grid.newpage()
  grid.draw(p4)
  dev.off()
}
