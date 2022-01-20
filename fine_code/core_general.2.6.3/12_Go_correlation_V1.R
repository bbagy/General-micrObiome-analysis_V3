
#' A Go_correlation
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_correlation()
#' May 22 2020



Go_correlation <- function(project, metaData, metabolicTab, abundTab, map, method, xanlgle, ncol,name, orders, height, width){
  # output
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  # map 정리 2
  sel.meta <- intersect(rownames(metabolicTab), rownames(map)); head(sel, "3")
  metabolicTab.sel <- metabolicTab[sel.meta,, drop=F];dim(metabolicTab.sel)
  
  sel.abun <- intersect(rownames(abundTab), rownames(map)); head(sel, "3")
  abundTab.sel <- abundTab[sel.abun,, drop=F];dim(abundTab.sel)
  
  
  
  x <- log((abundTab.sel+1)/(rowSums(abundTab.sel)+dim(abundTab.sel)[2]))
  x <- x[,order(colSums(x),decreasing=TRUE)]
  y <- metabolicTab.sel
  

  
  if(length(name) == 1){
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.%s.pdf",out_path,project, method,name,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }else{
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.pdf",out_path,project, method,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }
  
  
  for (des in rownames(subset(metadata, Go_correlation=="yes"))){
    groups<-map[,des]
    #Now calculate the correlation between individual Taxa and the environmental data
    df<-NULL
    for(i in colnames(x)){
      for(j in colnames(y)){
        for(k in unique(groups)){
          a <- x[groups==k,i,drop=F]
          b <- y[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
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
    colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
    df$Pvalue<-as.numeric(as.character(df$Pvalue))
    df$AdjPvalue<-rep(0,dim(df)[1])
    df$Correlation<-as.numeric(as.character(df$Correlation))
    
    adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
    adjustment<-5
    
    if(adjustment==1){
      df$AdjPvalue<-df$Pvalue
    } else if (adjustment==2){
      for(i in unique(df$Env)){
        for(j in unique(df$Type)){
          sel<-df$Env==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==3){
      for(i in unique(df$Taxa)){
        for(j in unique(df$Type)){
          sel<-df$Taxa==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==4){
      for(i in unique(df$Taxa)){
        sel<-df$Taxa==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    } else if (adjustment==5){
      for(i in unique(df$Env)){
        sel<-df$Env==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
    
    #Now we generate the labels for signifant values
    
    df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("***", "**", "*", ""))
    
    #We ignore NAs
    df<-df[complete.cases(df),]
    
    #We want to reorganize the Env data based on they appear
    #df$Env<-factor(df$Env,as.character(df$Env))
    
    #We use the function to change the labels for facet_grid in ggplot2
    Env_labeller <- function(variable,value){
      return(sel_env_label[as.character(value),"Trans"])
    }
    
    df$Env <- factor(df$Env)
    df$Type <- factor(df$Type, levels = orders)
    p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
    p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
    p<-p+theme(axis.text.x = element_text(angle = xanlgle, hjust = 1, vjust=0.5))
    p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
    p<-p+facet_wrap (~  Env, scales="free_x", ncol = ncol)
    #p<-p+facet_wrap (~  Env, ncol = 11)
    p<- p+ ggtitle(sprintf("%s %s", des, "Ladas")) 
    print(p)
  }
  dev.off()
}

