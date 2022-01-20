# 200505 ladas genefamily 로제작
# table to deseq2

Go_deseq2_ntd <- function(psIN, project,metaData, taxRanks, adjust, name,  order, data_type,alpha){
  
  # package
  if(!'colorspace' %in% installed.packages()){
    install.packages('colorspace')
  }else{
    library('colorspace')
  }
  
  ranks <- taxRanks
  taxaname <- ranks
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_ntd <- file.path(sprintf("%s_%s/table/deseq2_ntd",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_ntd)) dir.create(out_ntd)
  
  # map 정리 2
  map <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(map)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  map.sel <- map[rownames(map), sel, drop=F];head(map.sel)
  

  
  for(t in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    else if (data_type == "other" | data_type == "Other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    
    # continue
    otu.filt[,taxaname[t]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks, level=taxaname[t])
    
    print(1)
    
    #-- for give taxa name part  "print("pass3")" --#
    funTab <- as.matrix(tax_table(psIN))
    rownames(funTab) <- funTab[,taxaname[t]]
    
    
    if (dim(otu.filt)[2] == 2){
      next
    }
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[t])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[t]]
    agg <- agg[,-1]
    rownames(agg) <- genera
    
    
    #saving table
    write.csv(agg, quote = FALSE, col.names = NA, file=sprintf("%s/%s.deseq2.ntd.%s.%s.csv", out_ntd, project,taxaname[t], format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
    
    
    
    # agg <- normalizeByCols(agg)
    # rownames(agg) <- genera
    #dim(agg)
    #ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    #agg <- agg[intersect(ftk,ftk),]
    #agg <- agg*10000
    # control data set after filter
    if (dim(agg)[1] == 0)
      next
    
    #agg[,taxaname[t]] <- rownames(agg)
    
    print("complete agg")
    # 중요 match map and gene table
    # remove NA row and 반올림 하여 정수 만들기 (integer)
    tab <- 1000000*agg
    tab.ceiling <- ceiling(tab[-c(99),])
    gene3 <- tab.ceiling # tab.ceiling+1
    
    
    #gene3 <- tab
    sel <- intersect(rownames(map.sel), colnames(gene3)); head(sel)
    gene4 <- gene3[,sel, drop=F];head(gene4)
    map.sel.sel <- map.sel[sel,, drop=F];head(map.sel.sel)
    
    print(sprintf("%s %s","table", dim(gene4)))
    print(sprintf("%s %s","map", dim(map.sel.sel)))
    
    print("pass1")
    res <- {}
    for (mvar in rownames(subset(metadata, Go_deseq2=="yes"))) {
      print(sprintf("Analyzing for %s",mvar))
      
      # map 정리 2
      if (metadata[mvar, "type"] == "factor") {
        map[,mvar] <- factor(map[,mvar])
      } else if (metadata[mvar, "type"] == "numeric") {
        map[,mvar] <- as.numeric(as.character(map.sel[,mvar]))
      } else if (metadata[mvar, "type"] == "date") {
        map[,mvar] <- as.Date(sprintf("%06d", map.sel[,mvar]), format="%m%d%y")
        map[,mvar] <- factor(as.character(map[,mvar]), levels=as.character(unique(sort(map.sel[,mvar]))))
      }
      
      #======== run deseq2 ==========#
      #dds <- DESeqDataSetFromMatrix(countData = gene4,  colData = map.sel.sel, design= as.formula(sprintf("~ %s", mvar)))
      # dds <- DESeq(dds)  기본구조 였음20210305 부터 사용안함
      # tmp <- results(dds) 기본구조 였음20210305 부터 사용안함
      
      #-- DESeq2 --#
      print("run deseq2")
      print("pass2")
      if (length(adjust) >= 1) {
        form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
        print(form)
        dds <- DESeqDataSetFromMatrix(countData = gene4,
                                      colData = map.sel.sel,
                                      design= form)

      } else {
        form = as.formula(sprintf("~ %s", mvar))
        dds <- DESeqDataSetFromMatrix(countData = gene4,
                                      colData = map.sel.sel,
                                      design= form)
        print(form)
      }
      
      # option 1
       gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
       }
       geoMeans = apply(counts(dds), 1, gm_mean)
       dds = estimateSizeFactors(dds, geoMeans = geoMeans)
       dds = estimateDispersions(dds)
       vst = getVarianceStabilizedData(dds)
       dds = DESeq(dds, fitType="mean") # local
       resultsNames(dds)

      
      # option 2
      # keep <- rowSums(counts(dds)) >= 10
      # dds <- dds[keep,]
      # dds = DESeq(dds)
      
      # paired analysis
      # make a comnination for stat
      map.sel.sel[,mvar] <- factor(map.sel.sel[,mvar], levels = orders)
      
      # normal transformation for heatmap
      print("ntd")
      print("pass7")
      ntd <- normTransform(dds)
      assay.ntd <- assay(ntd);assay.ntd
      res.ntd <- results(dds);res.ntd
      res.ntd$padj <- p.adjust(res.ntd$pvalue, method="fdr")
      
      dim(res.ntd)
      
      sub.res <- subset(res.ntd, pvalue < alpha);dim(sub.res)
      sub.res.uniq <- sub.res[!duplicated(rownames(sub.res)),];dim(sub.res.uniq)
      if (dim(sub.res.uniq)[1] == 1){
        next
      }
      # 이름 정리
      #sub.res.uniq$taxa <- gsub(".*:_", "", sub.res.uniq$taxa);head(sub.res.uniq$taxa)
      
      assay.ntd.sel <- assay.ntd[rownames(sub.res.uniq),]
    
      rownames(assay.ntd.sel) <- rownames(sub.res.uniq)
      
      
      
      
      ####  my_colours 아 진짜 내가 이걸 해내는 구만.. mvar별로 다른 색 입히기.
      #display.brewer.pal(6, "Set2")
      #display.brewer.pal(8, "Set3")
      #mvar <- rownames(subset(metadata, Go_deseq2=="yes"))
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
      
      
      df <- as.data.frame(colData(dds)[,mvar])

      if (length(rownames(subset(metadata, Go_deseq2=="yes"))) ==1){
        colnames(df) <- mvar
      }
      rownames(df) <- colnames(assay.ntd.sel)
      # 한개 이상 return 하기
      # multiple list 만들기
      functionReturningTwoValues <- function() {
        results <- list()
        results$ntd <- assay.ntd.sel
        results$cols <- my_colours
        results$df <- df
        return(results) 
      }
      cat("\n")
      print("$ntd, $df and $cols are returned.")
      result <- functionReturningTwoValues()
      
      # saveRDS for list
      
      if (!is.null(name)) {
        saveRDS(result, sprintf("%s/deseq2_ntd.%s.%s.%s.%s.%s.rds", out_ntd, project, mvar,taxaname[t],name,format(Sys.Date(), "%y%m%d")))
      }else{
        saveRDS(result, sprintf("%s/deseq2_ntd.%s.%s.%s.%s.rds", out_ntd, project, mvar,taxaname[t],format(Sys.Date(), "%y%m%d")))
      }

    }
  }
}

