# 200505 ladas genefamily 로제작
# table to deseq2

Go_deseq2_tab <- function(project, tab, map, metadata, name,alpha, height, width){
  
  # package
  if(!'colorspace' %in% installed.packages()){
    install.packages('colorspace')
  }else{
    library('colorspace')
  }
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # map 정리 2
  sel <- intersect(rownames(meta), colnames(map)); head(sel, "3")
  meta.sel <- meta[sel,, drop=F];head(meta.sel)
  map.sel <- map[rownames(map), sel, drop=F];head(map.sel)
  
  # 중요 match map and gene table
  # remove NA row and 반올림 하여 정수 만들기 (integer)
  tab <- tab*10000
  tab.ceiling <- ceiling(tab[-c(99),])
  gene3 <- tab.ceiling+1
  
  
  
  sel <- intersect(rownames(map.sel), colnames(gene3)); head(sel)
  gene4 <- gene3[,sel, drop=F];head(gene4)
  map.sel.sel <- map.sel[sel,, drop=F];head(map.sel.sel)
  
  print(sprintf("%s %s","table", dim(gene4)))
  print(sprintf("%s %s","map", dim(map.sel.sel)))
  
  res <- {}
 
  
  dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
  for (mvar in rownames(subset(meta, Go_deseq2=="yes"))) {

    print(sprintf("Analyzong for %s",mvar))
    
    # map 정리 2
    if (meta[mvar, "type"] == "factor") {
      map[,mvar] <- factor(map[,mvar])
      if (!(is.na(meta[mvar, "baseline"])) && meta[mvar, "baseline"] != "") {
        map[,mvar] <- relevel(map[,mvar], meta[mvar, "baseline"])
      }
    } else if (meta[mvar, "type"] == "numeric") {
      map[,mvar] <- as.numeric(as.character(map.sel[,mvar]))
    } else if (meta[mvar, "type"] == "date") {
      map[,mvar] <- as.Date(sprintf("%06d", map.sel[,mvar]), format="%m%d%y")
      map[,mvar] <- factor(as.character(map[,mvar]), levels=as.character(unique(sort(map.sel[,mvar]))))
    }
    
    # run deseq2
    dds <- DESeqDataSetFromMatrix(countData = gene4,
                                  colData = map.sel.sel,
                                  design= as.formula(sprintf("~ %s", mvar)))
    dds <- DESeq(dds)
    tmp <- results(dds)
    
    print(1)
    for (smvar in levels(map[,mvar])) {
      if(smvar == meta.sel[mvar, "baseline"] | smvar == "" )
        next
      basline <- meta.sel[mvar, "baseline"]
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))
      tmp$dir <- ifelse(tmp$padj < 0.05, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      tmp$mvar <- mvar
      res <- rbind(res, tmp)
      resSig <- as.data.frame(subset(tmp, padj<0.05)); resSig <- resSig[order(resSig$log2FoldChange),]
      # 다음이 왜 있는지 모르겠네. 멈주니까 잠시 지우자 
      # resSig$taxa <- factor(resSig$taxa, levels=resSig$taxa)
      
      # p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange)) + geom_bar(stat="identity", fill="#aaaaaa") + geom_text(aes(label=taxa), y=0, size=2, hjust=0.5) + coord_flip() + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + theme(axis.text.y=element_blank())
      #print(p)
      # forest plot of significant results (FDR adjusted for this variable only)
      
      print(2)
      
      lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
      p <- ggplot(resSig, aes(x=taxa, y=log2FoldChange, color=dir)) + geom_point() + geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("DESeq2 hits (%s)", mvar)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims))
      
      
      if(length(name) == 1){
        pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.%s.pdf",out_path, project, mvar,name, format(Sys.Date(), "%y%m%d")), height = height, width=width)
        print(p)
      } else{
        pdf(sprintf("%s/7_%s.forest.deseq2.ntd.%s.%s.pdf",out_path, project, mvar,format(Sys.Date(), "%y%m%d")), height = height, width=width)
        print(p)
      }
    }
    dev.off()
    
    # heatmap
    print(3)
    ntd <- normTransform(dds)
    assay.ntd <- assay(ntd)
    
    sub.res <- subset(res, padj < alpha)
    sub.res.uniq <- sub.res[!duplicated(sub.res$taxa),]
    # 이름 정리
    sub.res.uniq$taxa <- gsub(".*:_", "", sub.res.uniq$taxa);head(sub.res.uniq$taxa)
    
    assay.ntd.sel <- assay.ntd[rownames(sub.res.uniq),]
    rownames(assay.ntd.sel) <- sub.res.uniq$taxa
    
    dim(assay.ntd.sel)
    
    
    
    ####  my_colours 아 진짜 내가 이걸 해내는 구만.. mvar별로 다른 색 입히기.
    #display.brewer.pal(6, "Set2")
    #display.brewer.pal(8, "Set3")
    mvar <- rownames(subset(meta, Go_deseq2=="yes"))
    my_colours <- list()
    for(i in 1:length(mvar)){
      if (i==1){
        cols <- colorRampPalette(brewer.pal(2, "Set1"))
      } else if (i==2){
        cols <- colorRampPalette(brewer.pal(5, "Paired"))
      }else if (i==3){
        cols <- colorRampPalette(brewer.pal(8, "Set3"))
      }
      
      colours <- cols(length(unique(map[,mvar[i]])))
      names(colours) <- unique(map[,mvar[i]])
      my_colour <- list(x = colours)
      names(my_colour) <- mvar[i]
      my_colours <- append(my_colours,my_colour)
    }
    
    
    df <- as.data.frame(colData(dds)[,rownames(subset(meta, Go_deseq2=="yes"))])
    
    if (length(rownames(subset(meta, Go_deseq2=="yes"))) ==1){
      colnames(df) <- mvar
    }
    
    rownames(df) <- colnames(assay.ntd.sel)
    
    # 한개 이상 return 하기
  
  # multiple list 만들기

  }
  
  functionReturningTwoValues <- function() {
    results <- list()
    results$ntd[[mvar]] <- assay.ntd.sel
    results$cols[[mvar]] <- my_colours
    results$df[[mvar]] <- df
    return(results) 
  }
  cat("\n")
  print("$ntd, $df and $cols are returned.")
  functionReturningTwoValues()
}
