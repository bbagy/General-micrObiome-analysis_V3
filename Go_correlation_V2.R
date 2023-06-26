#' A Go_correlation
#'
#' This function allows you to do correlation anaylysis
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords filter
#' @export
#' @examples
#' Go_filter


Go_correlation <- function(psIN, project,
                           rank,
                           KW = 0.05,
                           maingroup=NULL,
                           top = NULL,
                           Table,
                           facet="TAB",
                           X="GROUP",
                           Y="PS",
                           method = "spearman",
                           xanlgle = 90,
                           ncols = 23,
                           name = NULL,
                           orders=NULL,
                           height,width
                           ){

  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  #=======================
  #      For ps table 
  #======================
  if(!is.null(top)){
    Top = names(sort(taxa_sums(psIN), TRUE)[1:top])
    ps.top = prune_taxa(Top, psIN);ps.top
  }else{
    ps.top = psIN
  }
  
  ### decide for  log transformation
  
  if(max(data.frame(otu_table(ps.top))) < 1){
    psIN.top <- ps.top
  }else{
    #ps.top.rel <- transform_sample_counts(ps.top, function(x) x / log2(x)) # log(1+x) 를 하면 NaN가 많이 나온다. 
    psIN.top <- transform_sample_counts(ps.top, function(x) x / sum(x)) # percent
  }
  
  
  
  otu.filt <- as.data.frame(t(otu_table(psIN.top))) 
  tt <- try(otu.filt[,rank]  <- getTaxonomy(otus=rownames(otu.filt), taxRanks = colnames(tax_table(psIN.top)), tax_tab=tax_table(psIN.top), level=rank),T)
  
  if(class(tt) == "try-error"){
    print("other table")
    otu.filt <- as.data.frame(otu_table(psIN.top)) 
    otu.filt[,rank]  <- getTaxonomy(otus=rownames(otu.filt), taxRanks = colnames(tax_table(psIN.top)), tax_tab=tax_table(psIN.top), level=rank)
  }else{
    otu.filt <- as.data.frame(t(otu_table(psIN.top))) 
    print("DADA2 table")
    otu.filt[,rank]  <- getTaxonomy(otus=rownames(otu.filt), taxRanks = colnames(tax_table(psIN.top)), tax_tab=tax_table(psIN.top), level=rank)
  }
  
  agg.ps <- aggregate(as.formula(sprintf(". ~ %s" , rank)), otu.filt, sum, na.action=na.pass)
  
  rankName <- agg.ps[,rank]
  rownames(agg.ps) <- rankName
  agg.ps[,rank] <- NULL
  
  
  ps_table <-t(agg.ps);dim(agg.ps)
  
  # Add grouping information
  map <- sample_data(psIN.top);dim(map)
  df <- data.frame(ps_table, Group = map[,maingroup]) #, name = map$StudyID,  NoOfFMT= map$NoOfFMT );head(df)
  
  #=========== KW test for function screening =========#
  
  df$Group <- NULL
  df[,maingroup] <- as.character(df[,maingroup]);df[,maingroup]
  df[,maingroup][df[,maingroup]==""] <- "NA";df[,maingroup]
  df.na <- subset(df, df[,maingroup] != "NA");df.na[,maingroup]  # subset 를 사용한 NA 삭제
  df.na[,maingroup] <- as.factor(df.na[,maingroup]);df.na[,maingroup]  
  
  
  set.seed(151) 
  group_1 <- as.factor(df.na[,maingroup]); group_1
  
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
  
  kw.sig <- kw[which(kw$p.value < KW),];dim(kw.sig)[1]
  
  cat(paste(dim(kw.sig)[1]," was p < ", KW,".","\n", sep=""))
  
  kw.mat <- as.matrix(kw.sig);dim(kw.mat)
  
  funcNames.sig <- kw.mat[,1];length(funcNames.sig)
  
  df.sel <- df.na[funcNames.sig]
  #df.sel <- data.frame(df.sel, Group = map[,maingroup]) 
  
  
  
  # genus <- df.sel[,maingroup]
  #rownames(df.sel) <- genus
  df.sel[,maingroup] <- NULL
  df.sel$Treatment.1 <- NULL
  df.sel_t <- data.frame(t(df.sel))
  df.sel_t$Func <- NULL
  
  class(df.sel_t)
  
  df.sel_t[] <- lapply(df.sel_t, function(x) as.numeric(as.character(x)))
  sapply(df.sel_t, class)
  #========== Cut number of taxa =======#
  
  if (dim(kw.sig)[1] > 20 ){
    # N <- dim(FuncTab)[1]
    
    psTab <- cbind(df.sel_t, total = rowSums(df.sel_t))
    psTab <- psTab[order(psTab$total ,  decreasing = TRUE), ];psTab$total
    
    psTab <- psTab[1:20,];psTab$total
    
    ps_table.kw <-t(psTab)
    print(1)
  } else if (dim(kw.sig)[1] == 1){
    next
    print(2)
  } else{
    ps_table.kw <- t(df.sel_t)
    print(3)
  }
  
  
  #============================
  #      For metabolite table
  #============================
  # try table type
  table <- Table
  #Table$Group <- NULL
  #agg.metab <- metabolite;dim(agg.metab)
  
  
  
  #=======================
  # map 정리 of function table and  metabolite 
  #======================
  rownames(ps_table.kw) <- gsub("X", "", rownames(ps_table.kw));head(rownames(ps_table.kw));dim(ps_table.kw)
  sel <- intersect(rownames(ps_table.kw), rownames(Table)); head(sel, "3")
  
  ps_table.sel <- ps_table.kw[sel,, drop=F];head(ps_table.sel);dim(ps_table.sel)

  table <- table[rownames(ps_table.sel),] 
  
  ps_table.sel <- ps_table.sel[rownames(ps_table.sel),] # 중요! 이것을 안하면 order에 따라 결과가 달라진다.
  
  
  sel_env <- colnames(ps_table.sel);sel_env 
  
  #NA 제거
  for (f in sel_env){
    ps_table.sel[,f] <- as.character(ps_table.sel[,f]);ps_table.sel[,f]
    ps_table.sel[,f][ps_table.sel[,f]==""] <- "NA";ps_table.sel[,f]
    # clin.sel[,maingroup]<- as.integer(clin.sel[,maingroup]);clin.sel[,maingroup]
    ps_table.sel[,f]<- as.numeric(ps_table.sel[,f]);ps_table.sel[,f]
  }
  # rownames NA 제거 
  ps_table.sel <- ps_table.sel[grepl("^NA", rownames(ps_table.sel))==F,]
  
  sel_env_label<-t(as.data.frame(sel_env))
  sel_env_label<-as.data.frame(sel_env_label)
  colnames(sel_env_label)<-c("Trans")
  sel_env_label$Trans<-as.character(sel_env_label$Trans)
  
  #Now get a filtered table based on sel_env
  ps_table_filtered <- ps_table.sel[,sel_env];head(ps_table_filtered) # meta_table to func_table
  table_filtered <- as.data.frame(table[rownames(table),]);head(table_filtered) # abund_table to  meta_table
  
  
  
  
  #Apply normalisation (either use relative or log-relative transformation)
  #x<-abund_table_filtered/rowSums(abund_table_filtered)
  
  #---------------------------#
  #---- group informaiton ----# kendall, spearman, or pearson 
  #---------------------------#

  # x, xy 는 data frame이어야 한다.
  # x <- log((abundTab+1)/(rowSums(abundTab)+dim(abundTab)[2]))
  # change class to as.numeric
  # sapply(ps_table_filtered, class)
  # sapply(table_filtered, class)
  

  
  #ps_table_filtered[] <- lapply(ps_table_filtered, function(x) as.numeric(as.character(x)))
  table_filtered[] <- lapply(table_filtered, function(x) as.numeric(as.character(x)))
  
  x <- table_filtered[,order(colSums(table_filtered),decreasing=TRUE)];dim(table_filtered)
  
  # y <- meta_table_filtered
  y <- ps_table_filtered;dim(y)
  mapping.sel <- mapping
  x <- data.frame(x)
  y <- data.frame(y)
  
  y[] <- lapply(y, function(x) as.numeric(as.character(x)))
  
  sapply(y, class)
  
  groups <- data.frame(map)[,maingroup]
  
  pdf(sprintf("%s/Correlation_%s.%s.%s.(KW=%s).%s%s.pdf",
              out_path,
              rank,
              method,
              maingroup,
              KW,
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width=width)
  
  for (des in c(maingroup)){
    #groups <- map[,des] # group 정보는 초기 clin에서 가져오자 
    # Now calculate the correlation between individual Taxa and the environmental data
    
    df<-NULL
    for(i in colnames(x)){
      for(j in colnames(y)){
        for(k in unique(groups)){
          a <- x[unique(groups)==k,i,drop=F]
          b <- y[unique(groups)==k,j,drop=F]
          tmp <- c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),
                   cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value, k)
          if(is.null(df)){
            df<-tmp  
          }
          else{
            df<-rbind(df,tmp)
          }    
        }
      }
    }
    
    df<-data.frame(row.names=NULL,df)
    
    
    colnames(df)<-c("TAB","PS","Correlation","Pvalue","GROUP") # c("Taxa","Env","Correlation","Pvalue","Type")
    df$Pvalue<-as.numeric(as.character(df$Pvalue))
    df$AdjPvalue<-rep(0,dim(df)[1])
    df$Correlation<-as.numeric(as.character(df$Correlation))
    
    unique(df$Correlation)
    
    
    #You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
    # 1 -> donot adjust
    # 2 -> adjust Env + Type (column on the correlation plot)
    # 3 -> adjust Taxa + Type (row on the correlation plot for each type)
    # 4 -> adjust Taxa (row on the correlation plot)
    # 5 -> adjust Env (panel on the correlation plot)
    adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
    adjustment<-5
    
    if(adjustment==1){
      df$AdjPvalue<-df$Pvalue
    } else if (adjustment==2){
      for(i in unique(df$PS)){
        for(j in unique(df$GROUP)){
          sel<-df$PS==i & df$GROUP==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==3){
      for(i in unique(df$TAB)){
        for(j in unique(df$GROUP)){
          sel<-df$TAB==i & df$GROUP==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==4){
      for(i in unique(df$TAB)){
        sel<-df$TAB==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    } else if (adjustment==5){
      for(i in unique(df$PS)){
        sel<-df$PS==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
    
    #Now we generate the labels for signifant values
    
    # df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("**", "*", "", ""))
    df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("**", "*", "", ""))
    
    #We ignore NAs
    df<-df[complete.cases(df),]
    
    #We want to reorganize the Env data based on they appear
    #df$Env<-factor(df$Env,as.character(df$Env))
    
    #We use the function to change the labels for facet_grid in ggplot2
    Env_labeller <- function(variable,value){
      return(sel_env_label[as.character(value),"Trans"])
    }
    
    # order
    if (length(orders) >= 1) {
      df$GROUP <- factor(df$GROUP , levels = orders)
    }
    else {
      df$GROUP  <- factor(df$GROUP)
    }

  
    
    
    
    df <- df[order(df$PS ,  decreasing = F), ]
    df$PS <- as.factor(df$PS)
    class(df$PS)
    

    p <- ggplot(aes_string(x=X, y=Y, fill="Correlation"), data=df) + theme_bw() + theme(strip.background = element_blank())
    p <- p + geom_tile() + scale_fill_gradient2(low="#0C5BB0FF", mid="white", high="#EE0011FF") 
    p <- p + geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
    p <- p + ggtitle(sprintf("%s", maingroup)) 
    
    p <- p + theme(axis.text.x = element_text(angle = xanlgle, size=8,hjust = 1, vjust=0.5,face = "italic"),
                   axis.text.y = element_text(angle=0, vjust=0.5, hjust=1))
    
  
    
    p <- p+ facet_wrap(as.formula(sprintf("~ %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", ncol = ncols) 
    
    # p<-p+facet_wrap(~  TAB, scales="free_x", ncol = ncols)
    
    print(p)
    dev.off()
    
    print("For method you can use kendall, spearman, or pearson")
    print("GROUP, PS, and TAB are only available for rearranging the plot xy-axis.")
  }
}

