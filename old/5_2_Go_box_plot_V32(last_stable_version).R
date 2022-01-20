#' A Go_box_plot
#'
#' This function allows you to express your love of cats.
#' @paroc love Do you love cats? Defaults to TRUE.
#' @keywords basic statistic and simple plots
#' @export
#' @exocples
#' Go_box_plot


Go_box_plot <- function(df, metaData, project, orders=NULL, outcomes, Chao1_Normality="no", Shannon_Normality="no", star="no",title= NULL, facet= NULL, paired=NULL, name= NULL, xanlgle=90,  height, width, plotCols, plotRows){

  if(!is.null(dev.list())) dev.off()
  # plot color
  colorset = "Dark2" # Dark2 Set1 Paired
  

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
  if (length(facet) >= 1) {
    if (!is.null(name)) {
      pdf(sprintf("%s/4_box.%s.%s.%s.%s.pdf",out_path, project,facet,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s/4_box.%s.%s.%s.pdf",out_path,project,facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (!is.null(name)) {
      pdf(sprintf("%s/4_box.%s.%s.%s.pdf",out_path,project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s/4_box.%s.%s.pdf",out_path,project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }

  ## fix factor  and  numeric
  df$etc <- NULL
  for (var in rownames(subset(metadata, Go_box =="yes"))) {
    if (metadata[var, "type"] == "factor") {
      df[,var] <- factor(df[,var])
    } else if (metadata[var, "type"] == "numeric") {
      df[,var] <- as.numeric(as.character(df[,var]))
    }
  }
  
  
  # plot
  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_box =="yes"))) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}
    
    # remove Na
    adiv <- data.frame(df)
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
    
    # make a comnination for stat
    cbn <- combn(x = levels(adiv.na[,mvar]), m = 2)

    my_comparisons <- {}
    for(i in 1:ncol(cbn)){
      x <- cbn[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    # check normality for Chao1 and Shannon
    for(oc in outcomes){
      if (oc == "Chao1"){
        if ( Chao1_Normality  == "no"  & length(unique(adiv.na[,mvar])) > 2 ) {
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          testmethod <- "t.test"
        }

        if (Chao1_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          testmethod <- "wilcox.test"
        } else if (Chao1_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          testmethod <- "t.test"
        }       
      } else if (oc == "Shannon") {
        if ( Shannon_Normality  == "no" & length(unique(adiv.na[,mvar])) > 2 ) {
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) > 2 ) {
          testmethod <- "t.test"
        }
        
        if (Shannon_Normality  == "no" & length(unique(adiv.na[,mvar])) < 3 ) {
          testmethod <- "wilcox.test"
        } else if (Shannon_Normality  == "yes" & length(unique(adiv.na[,mvar])) < 3 ) {
          testmethod <- "t.test"
        }
      } else{
        testmethod <- "wilcox.test"
      }
      
      # re-order
      if (length(orders) >= 1) {
        adiv.na[,mvar] <- factor(adiv.na[,mvar], levels = orders)
      } else {
        adiv.na[,mvar] <- factor(adiv.na[,mvar])
      }
      
      # remove NA for facet
      if (length(facet) >= 1) {
        for (fc in facet){
          adiv.na[adiv.na[,fc] == ""] <- "NA"
          adiv.na.sel <- adiv.na[!is.na(adiv.na[,fc]), ]
          adiv.na <- adiv.na.sel 
          # facet or not
          adiv.na[,fc] <- factor(adiv.na[,fc], levels = orders)
        }
      }

      
      p1 <- ggplot(adiv.na, aes_string(x=mvar, y=oc, colour=mvar)) + theme_bw() + labs(y=oc, x=NULL) +
        theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5)) +
        scale_color_brewer(palette=colorset)
      
      # Close an image
        if (!is.null(title)) {
          p1 <- p1 + ggtitle(title)
        } else{
          p1 <- p1 + ggtitle(sprintf("%s", mvar))
        }
      
      if (star == "no") {  
      p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
      }  else if (star == "yes") {
        p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
      }
      
      #paired plot type
      
       if (!is.null(paired)) {
        #p1 = p1 + geom_point(size = 1) 
        p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3) 
        p1 = p1 + geom_point(aes_string(shape=paired)) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
        p1 = p1 + guides(color = FALSE, size = FALSE) + theme(legend.title = element_blank(), 
                                                              legend.position="bottom", 
                                                              legend.justification="left",
                                                              legend.box.margin = ggplot2::margin(0,0,0,-1,"cm")) 
      }  else{
        p1 = p1 + geom_boxplot(outlier.shape = NA)  + theme(legend.position="none") 
        p1 = p1 + geom_jitter(shape=16, alpha = 3, size = 1.5, position=position_jitter(0.2)) # alpha=0.3
      } 
      

      

      
      # facet
      if (length(facet) >= 1) {
        facetCol <- length(unique(adiv[,facet]))
        p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
        
      } else {
        p1 = p1 
      }
      
      plotlist[[length(plotlist)+1]] <- p1 
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}

