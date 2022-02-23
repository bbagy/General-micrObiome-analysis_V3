#' A Go_box_plot
#'

Go_dualYplot <- function(df, TaxTab, metaData, project, orders=NULL, Box, Line1, Line2=NULL,
                       title= NULL, name= NULL, 
                       xanlgle=90,  height, width){
  
  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  
  # out file

  if (!is.null(name)) {
    pdf(sprintf("%s/4_dualYplot.%s.%s.%s.pdf",out_path,project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } 
  else {
    pdf(sprintf("%s/4_dualYplot.%s.%s.pdf",out_path,project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  ## fix factor  and  numeric
  df$etc <- NULL
  df2 <- read.csv(sprintf("%s",TaxTab),header=T,row.names=1,check.names=FALSE);head(df2)
  
  rownames(df2)<- df2$Species
  df2$Species <-NULL
  rownames(df2) <- gsub(" ","_",rownames(df2));rownames(df2)
  
  
  df2 <- as.data.frame(t(df2))
  

  
  
  for (var in rownames(subset(metadata, Go_box =="yes"))) {
    print(var)
    if (metadata[var, "type"] == "factor") {
      df[,var] <- factor(df[,var])
    } else if (metadata[var, "type"] == "numeric") {
      df[,var] <- as.numeric(as.character(df[,var]))
    }
  }
  

  
  # plot
  for (mvar in rownames(subset(metadata, Go_box =="yes"))) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    
    # merge adiv and taxa table
    adiv <- merge(df, df2, by="row.names");head(adiv)
    
    # NA remove
    adiv[,mvar] <- as.character(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    adiv.na[,mvar] <- as.factor(adiv.na[,mvar]);adiv.na[,mvar]  
    
    
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(adiv.na)[1], dim(adiv)[1]))
    
    if (length(unique(adiv.na[,mvar])) ==1) {
      next
    }
    
    summary.adiv.na <- summary(adiv.na[,mvar])
    
    # re-order
    if (length(orders) >= 1) {
      adiv.na[,mvar] <- factor(adiv.na[,mvar], levels = orders)
    } else {
      adiv.na[,mvar] <- factor(adiv.na[,mvar])
    }
    
    #===============================#
    # Visualization for Dual Y axis #
    #===============================#
  
    # for Line1
    mean.line1 <- aggregate(adiv.na[,Line1], list(adiv.na[,mvar]), FUN=mean)
    colnames(mean.line1) <- c(mvar, Line1);mean.line1
    mean.line1[,Line1] <- mean.line1[,Line1]*10

    
    
    p <- ggplot() + theme_bw() + theme(strip.background = element_blank()) + #theme_ipsum() +
      geom_boxplot(data=adiv.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), outlier.shape = NA, show.legend = FALSE) +
      theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5)) +
      # theme(legend.position="none") +
      scale_color_manual(NULL, values = mycols) 
    
    p1 <- p + geom_line(data = mean.line1, 
                        mapping = aes(x = !!sym(mvar), y = !!sym(Line1), group=1, linetype=""), 
                        inherit.aes = FALSE,color="#FF9DA7", size=1) + 
      scale_linetype_manual(NULL,labels = Line1, values = 1) 
    

      
    
    
    # for Line2
    if (!is.null(Line2)){
      mean.line1 <- aggregate(adiv.na[,Line1], list(adiv.na[,mvar]), FUN=mean)
      colnames(mean.line1) <- c(mvar, Line1);mean.line1
      mean.line1[,Line1] <- mean.line1[,Line1]*10
      
      
      mean.line2 <- aggregate(adiv.na[,Line2], list(adiv.na[,mvar]), FUN=mean)
      colnames(mean.line2) <- c(mvar, Line2);mean.line1
      mean.line2[,Line2] <- mean.line2[,Line2]*10
      
      mean.line <- merge(mean.line1, mean.line2, by=mvar);head(mean.line)
      
      mean.line.melt <- melt(mean.line)

      
      

      p1 <- p + geom_line(data = mean.line.melt, 
                          mapping = aes(x = !!sym(mvar), y = value, group=variable,color=variable), 
                          inherit.aes = FALSE, size=1) + 
        scale_linetype_manual(NULL, values = 1) 
    }
    
    
    
    p1 <- p1 + scale_y_continuous(sec.axis = sec_axis(~.*10, name="Relative abundance (%)")) 
    
    if (!is.null(title)) {
      p1 <- p1 + ggtitle(title)
    } else{
      p1 <- p1 + ggtitle(sprintf("%s", mvar))
    }
    print(p1)
  }
  dev.off()
}

