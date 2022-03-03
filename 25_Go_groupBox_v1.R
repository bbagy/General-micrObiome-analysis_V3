#' A Go_groupBox
#'

Go_groupBox <- function(psIN, mainGroup, project, orders=NULL, top=NULL, name =NULL, rank, cutoff, color=NULL, ylim=NULL,flip,height, width){
  
  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  
  
  
  
  if(!is.null(top)){
    Top = names(sort(taxa_sums(psIN), TRUE)[1:top])
    ps.top = prune_taxa(Top, psIN);ps.top
  }else{
    ps.top = psIN
  }
  
  
  ### log transformation
  ps.top.rel <- transform_sample_counts(ps.top, function(x) x / log2(x)) # log(1+x) 를 하면 NaN가 많이 나온다.
  
  tab = data.frame(otu_table(ps.top.rel))
  
  nsamps_threshold <- 0.01 # fraction of relabund to call a sample positive
  filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
  nperm <- 100000
  
  ### aggregation by rank
  otu.filt <- as.data.frame((otu_table(ps.top.rel))) # for dada2  t(otu_table(ps.relative)
  otu.filt$func <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.top.rel), taxRanks =colnames(tax_table(psIN)), level= rank)
  agg <- aggregate(. ~ func, otu.filt, sum, na.action=na.pass);dim(agg)
  
  
  
  funcNames <- agg$func
  agg <- agg[,-1]
  rownames(agg) <- funcNames
  ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
  agg <- agg[intersect(ftk,ftk),]
  agg$Taxa <- rownames(agg)
  # agg_t <- t(agg)
  
  
  
  aggDf<-as.data.frame(agg, row.names = agg$Taxa)
  aggDf <- aggDf[,-1]
  agg_t <- t(aggDf)
  
  colnames(agg_t)
  rownames(agg_t)
  
  # Add grouping information
  map <- sample_data(ps.top.rel);dim(map)
  df <- data.frame(agg_t, Group = map[,mainGroup]) #, name = map$StudyID,  NoOfFMT= map$NoOfFMT );head(df)
  
  df[,mainGroup] <- as.character(df[,mainGroup]);df[,mainGroup]
  df[,mainGroup][df[,mainGroup]==""] <- "NA";df[,mainGroup]
  df.na <- subset(df, df[,mainGroup] != "NA");df.na[,mainGroup]  # subset 를 사용한 NA 삭제
  df.na[,mainGroup] <- as.factor(df.na[,mainGroup]);df.na[,mainGroup]  
  
  

  group_1 <- as.factor(df.na[,mainGroup]); group_1
  
  # N <- 20
  # taxaname<-colnames(df)[1:N];taxaname
  
  
  #df.na1 <- subset(df.na, select = c(length(colnames(agg_t))) )
  #df.na <- df.na[1:N];df.na
  
  
  kruskal.wallis.table <- data.frame()
  for (i in 1:dim(df.na)[2]) {
    ks.test <- kruskal.test(df.na[,i], g=group_1)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                  data.frame(id=names(df.na)[i],
                                             p.value=ks.test$p.value
                                  ))
    # Report number of values tested
    cat(paste("Kruskal-Wallis test for ",names(df.na)[i]," ", i, "/", 
              dim(df.na)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
  }
  
  
  
  kw <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing = FALSE), ] # decreasing, increasing
  
  kw.sig <- kw[which(kw$p.value < cutoff),];dim(kw.sig)[1]
  
  cat(paste(dim(kw.sig)[1]," was p < ", cutoff,".","\n", sep=""))
  
  kw.mat <- as.matrix(kw.sig);dim(kw.mat)
  
  funcNames.sig <- kw.mat[,1];length(funcNames.sig)
  
  df.sel <- df.na[funcNames.sig]
  df.sel <- data.frame(df.sel, Group = map[,mainGroup]) 
  

  
  
  df.sel.melt <- melt(df.sel, id.vars = mainGroup, measure.vars = funcNames.sig)
  df.sel.melt$value <- as.numeric(df.sel.melt$value)
  df.sel.melt.clean <- subset(df.sel.melt, variable != "Group" &  value > 0)
  
  
  
  
  if (!is.null(orders)) {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup], levels = rev(orders))
  } else {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup])
  }
  
  df.sel.melt.clean$variable <- as.character(df.sel.melt.clean$variable)
  
  df.sel.melt.clean <- df.sel.melt.clean[order(df.sel.melt.clean$variable ,  decreasing = F), ]
  
  
  p <- ggplot(df.sel.melt.clean, aes_string(x="variable", y="value", fill=mainGroup)) +  geom_boxplot(outlier.shape = NA,lwd=0.3) + 
    theme_bw() + theme(strip.background = element_blank()) + 
    labs(y="Relative abundance (log2)", x= NULL) + ggtitle(sprintf("kruskal wallis p < %s",cutoff))
  
  # + stat_compare_means(aes_string(group = mainGroup),label = "p.format") + 
    
  #+ scale_x_discrete(limits = rev)
  
  
  if(!is.null(color)){
    p <- p + scale_fill_manual(values = color)
  }else{
    p <- p
  }

  if(flip == T){
    p <- p+ coord_flip()
  }else{
    p <- p + theme(text=element_text(size=9), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  }
  

  
  # tt <- try(mycols, T)
  # if(class(tt) == "try-error"){
  #  p <- p
  # }else{
  #   p <- p + scale_fill_manual(values = mycols)
  # }

  if(!is.null(ylim)){
    p = p + ylim(ylim[1] , ylim[2])
  }else{
    p=p
  }
  
  #=== image size ===#
  #height <- 0.4*length(unique(df.sel.melt.clean[,mainGroup])) + 0.4*dim(kw.sig)[1];height
  #width <- log((max(nchar(funcNames.sig)))*max(nchar(as.character(unique(df.sel.melt.clean[,mainGroup])))));width
  print(p)
  
  pdf(sprintf("%s/groupBox.%s.%s.%s%s%s.pdf", out_path, 
              project, 
              mainGroup,
              ifelse(is.null(rank), "", paste(rank, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  print(p)
  
  dev.off()
}








