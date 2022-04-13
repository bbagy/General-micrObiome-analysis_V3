#' A Go_box_plot
#'

Go_boxplot <- function(df, cate.vars, project, outcomes,
                       orders=NULL, 
                       mycols=NULL, 
                       combination=NULL,
                       ylim =NULL,
                       title= NULL, 
                       facet= NULL, 
                       paired=NULL, 
                       name= NULL, 
                       statistics = "yes", 
                       parametric= "no", 
                       star="no",
                       xanlgle=90,  
                       height, width, plotCols, plotRows){
  
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  pdf(sprintf("%s/box.%s.%s%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # plot
  plotlist <- list()
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}
    
    # remove Na
    df <- data.frame(df)
    
    df[,mvar] <- as.character(df[,mvar]);df[,mvar]
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    df.na[,mvar] <- as.factor(df.na[,mvar]);df.na[,mvar]  
    
    df.na[,mvar] <- factor(df.na[,mvar], levels = intersect(orders, df.na[,mvar]))
    
    
    
    
    # Add number of samples in the group
    renamed_levels <- as.character(levels(df.na[,mvar]));renamed_levels
    oldNames <- unique(df.na[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(df.na[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(df.na[[mvar]])[levels(df.na[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(df.na)[1], dim(df)[1]))
    
    if (length(unique(df.na[,mvar])) ==1) {
      next
    }
    
    summary.df.na <- summary(df.na[,mvar])
    
    #------------------------------#
    # for group combination or not #
    #------------------------------#
    
    
    if (!is.null(combination)){
      group.cbn <- combn(x = levels(df.na[,mvar]), m = combination)
      
      #print(count(group.cbn))
      
      group_comparisons <- {}
      for(i in 1:ncol(group.cbn)){
        x <- group.cbn[,i]
        group_comparisons[[i]] <- x
      };group_comparisons
      
      print(1)
      for(i in 1:length(group_comparisons)){
        print(group_comparisons[i])
        group.combination <- unlist(group_comparisons[i]);group.combination
        
        if(combination ==2){
          basline <- group.combination[1]
          smvar <- group.combination[2]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar)) 
        } else if(combination ==3){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2)) 
        }else if(combination ==4){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3)) 
        }else if(combination ==5){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4)) 
        }else{
          print("combination should be 2, 3, 4 and 5 only.")
          break
        }
        
        unique(df.cbn[,mvar])
        
        
        # make a comnination for stat
        df.cbn[,mvar] <- factor(df.cbn[,mvar])
        cbn <- combn(x = levels(df.cbn[,mvar]), m = 2)
        
        
        my_comparisons <- {}
        for(i in 1:ncol(cbn)){
          x <- cbn[,i]
          my_comparisons[[i]] <- x
        };my_comparisons
        
        if(combination != 2){
          combination.N <- combination - 1
          my_comparisons <- my_comparisons[1:combination.N]
        }
        
        
        
        # check statistics method
        for(oc in outcomes){
          if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
            if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
              testmethod <- "t.test"
            }else{
              testmethod <- "wilcox.test"
            }
          } 
          
          
          
          # remove NA for facet
          if (length(facet) >= 1) {
            for (fc in facet){
              df.cbn[,fc] <- as.character(df.cbn[,fc]);df.cbn[,fc]
              df.cbn[,fc][df.cbn[,fc] == ""] <- "NA"
              df.cbn.sel <- df.cbn[!is.na(df.cbn[,fc]), ]
              df.cbn <- df.cbn.sel 
              # facet or not
              df.cbn[,fc] <- factor(df.cbn[,fc], levels = orders)
            }
          }
          
  
          p1 <- ggplot(df.cbn, aes_string(x=mvar, y=oc, colour=mvar))  + labs(y=oc, x=NULL) + 
            theme_bw() + theme(strip.background = element_blank()) +
            theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
                  plot.title=element_text(size=9)) # ,face="bold"

          
          if(!is.null(mycols)){
            p1 <- p1 + scale_color_manual(values = mycols)
          }else{
            p1 <- p1
          }
          
        
          if (!is.null(title)) {
            p1 <- p1 + ggtitle(title)
          } else{
            p1 <- p1 + ggtitle(sprintf("%s", mvar))
          }
          
          if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
            if (star == "no") {  
              p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
            }  else if (star == "yes") {
              p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
            }
          }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
            p1 <- p1 
          }
          
          # plot design
          if (height*width <= 6){
            dot.size = 0.7
            box.tickness = 0.3
          }else if (height*width > 6 & height*width < 10){
            dot.size = 1
            box.tickness = 0.4
          }else{
            dot.size = 1.5
            box.tickness = 0.5
          }
          
          if(oc != "Chao1"){
            if(!is.null(ylim)){
              p1 = p1 + ylim(ylim[1] , ylim[2])
            }else(
              p1=p1
            )
          }
          # paired plot type
          if (!is.null(paired)) {
            #p1 = p1 + geom_point(size = 1) 
            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
            p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3),show.legend = F)  #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
            p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3)) 
            p1 = p1  + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm")) 
            
          }  else{
            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
            
            
            # count or table for number of variable
            if (max(table(df.cbn[,mvar])) > 250 & max(table(df.cbn[,mvar])) < 500){
              dot.size <- dot.size/2
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            } else  if (max(table(df.cbn[,mvar])) < 250 ){
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            }else if(max(table(df.cbn[,mvar])) > 500) {
              dot.size <- dot.size/3
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))
            }
            
            
          } 
          
          # facet
          if (length(facet) >= 1) {
            facetCol <- length(unique(df[,facet]))
            p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          } else {
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          }
          
          plotlist[[length(plotlist)+1]] <- p1 
        }
      }
    }else{
      # make a comnination for stat
      cbn <- combn(x = levels(df.na[,mvar]), m = 2)
      
      my_comparisons <- {}
      for(i in 1:ncol(cbn)){
        x <- cbn[,i]
        my_comparisons[[i]] <- x
      };my_comparisons
      # check statistics method
      for(oc in outcomes){
        if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
          if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
            testmethod <- "t.test"
          }else{
            testmethod <- "wilcox.test"
          }
        } 
        
        
        
        # remove NA for facet
        if (!is.null(facet)) {
          for (fc in facet){
            df.na[,fc] <- as.character(df.na[,fc]);df.na[,fc]
            df.na[,fc][df.na[,fc] == ""] <- "NA"
            df.na.sel <- df.na[!is.na(df.na[,fc]), ]
            df.na <- df.na.sel 
            # facet or not
            df.na[,fc] <- factor(df.na[,fc], levels = orders)
          }
        }
        
        
        p1 <- ggplot(df.na, aes_string(x=mvar, y=oc, colour=mvar))  + labs(y=oc, x=NULL) + 
          theme_bw() + theme(strip.background = element_blank()) +
          theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
                plot.title=element_text(size=9))  #,face="bold"


        # scale_color_brewer(palette=colorset)
        
        if(!is.null(mycols)){
          p1 <- p1 + scale_color_manual(values = mycols)
        }else{
          p1 <- p1
        }
        
        
        # Close an image
        if (!is.null(title)) {
          p1 <- p1 + ggtitle(title)
        } else{
          p1 <- p1 + ggtitle(sprintf("%s", mvar))
        }
        
        if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
          if (star == "no") {  
            p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
          }  else if (star == "yes") {
            p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
          }
        }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
          p1 <- p1 
        }
        # plot design
        if (height*width <= 6){
          dot.size = 0.7
          box.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 1
          box.tickness = 0.4
        }else{
          dot.size = 1.5
          box.tickness = 0.5
        }
        
        # y axis limit      
        if(oc == "Shannon"){
          if(!is.null(ylim)){
            p1 = p1 + ylim(ylim[1] , ylim[2])
          }else{
            p1=p1
          }
        }
        
        # paired plot type
        if (!is.null(paired)) {
          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)   #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
          p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3)) 
          p1 = p1  + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm")) 
          
        } else{
          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          
            # count or table for number of variable
            if (max(table(df[,mvar])) > 250 & max(table(df[,mvar])) < 500){
              dot.size <- dot.size/2
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            } else  if (max(table(df[,mvar])) < 250 ){
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            }else if(max(table(df[,mvar])) > 500) {
              dot.size <- dot.size/3
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))
            }
        } 
        # facet
        if (length(facet) >= 1) {
          facetCol <- length(unique(df[,facet]))
          p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        } else {
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        }
        
        plotlist[[length(plotlist)+1]] <- p1 
      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
