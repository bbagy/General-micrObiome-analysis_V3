#' A Go_bdiv
#'


Go_bdiv <- function(psIN, category.vars, project, orders, mycols=NULL,combination=NULL,
distance_metrics, plot="PCoA", shapes = NULL, ID = NULL, ellipse="yes", facet=NULL, name=NULL, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

 # out file
   # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  pdf(sprintf("%s/ordi.%s.%s%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  
  

  
  plotlist <- list()
  for (mvar in category.vars) {
    mapping <- data.frame(sample_data(psIN))

    mapping[,mvar] <- factor(mapping[,mvar])

    sample_data(psIN) <- mapping
    
    
    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}
    
    if (length(shapes) >= 1){
      if (shapes == mvar){
        next
      }
    } else {}
    
    #------------------------------#
    # for group combination or not #
    #------------------------------#
    
    if (!is.null(combination)){
      mapping[,mvar] <- factor(mapping[,mvar], levels = intersect(orders, mapping[,mvar]))
      
      
      group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
      
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
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar)) 
        } else if(combination ==3){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar1, smvar2)) 
        }else if(combination ==4){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar1, smvar2,smvar3)) 
        }else{
          print("combination should be 2, 3, and 4 only.")
          break
        }
        
        psIN.cbn <- psIN
        sample_data(psIN.cbn) <- mapping.cbn
        for(distance_metric in distance_metrics){
          # remove na
          mapping.sel <- data.frame(sample_data(psIN.cbn))
          mapping.sel[mapping.sel==""] <- "NA"
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          na.count <- length(mapping.sel.na)
          psIN.cbn.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.cbn.na ))
          
          
          if (!is.null(facet)) {
            print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                          facet,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          } else{
            print(sprintf("##-- %s (total without NA: %s/%s) --##",
                          mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          }
          
          
          
          ## fix factor  and  numeric
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
          
          
          ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
          plist = llply(as.list(ord_meths), function(i, psIN.cbn.na, distance_metric){
            ordi = ordinate(psIN.cbn.na, method=i, distance=distance_metric)
            plot_ordination(psIN.cbn.na, ordi, type = "samples", color= mvar)
          }, psIN.cbn.na, distance_metric)
          
          
          
          names(plist) <- ord_meths
          
          pdataframe = ldply(plist, function(x){
            df = x$data[, 1:2]
            colnames(df) = c("Axis_1", "Axis_2")
            return(cbind(df, x$data))
          })
          names(pdataframe)[1] = "method"
          
          pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)
          
          pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)
          
          
          # Plots
          p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
          
          
          if (!is.null(shapes)) {
            
            pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
            p = p +  geom_point(aes_string(shape=shapes), size=1.5, alpha = 3) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
            
          }else{
            p = p + geom_point(size=1.5, alpha = 3)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
          }
          
          p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
          p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
          p = p + theme(legend.position = "bottom", 
                        legend.title = element_blank(),
                        legend.justification="left", 
                        legend.box = "vertical",
                        legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))
          
          
          if(!is.null(mycols)){
            p <- p + scale_color_manual(values = mycols)
          }else{
            p <- p
          }
          
          
          
          # ID variation
          if (!is.null(ID)) {
            p = p + geom_text_repel(aes_string(label = ID), size = 2)
          } else {
            p = p 
          }
          
          # ellipse variation
          if (ellipse == "yes" | ellipse == "Yes" ) {
            p = p + stat_ellipse(type = "norm", linetype = 2) 
          } else if (ellipse == "no" | ellipse == "No" ){
            p = p 
          }
          
          if (!is.null(facet)) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }
          else {
            p = p
          }
          
          #plotlist[[length(plotlist)+1]] <- p
          
          print(p)
          
        }
      }
    }    else{
      for(distance_metric in distance_metrics){
        # remove na
        mapping.sel <- data.frame(sample_data(psIN))
        mapping.sel[mapping.sel==""] <- "NA"
        mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
        na.count <- length(mapping.sel.na)
        psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
        mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
        
        
        if (!is.null(facet)) {
          print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                        facet,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
        } else{
          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
        }
        
        
        
        ## fix factor  and  numeric
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
        
        
        
        ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
        plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
        }, psIN.na, distance_metric)
        
        names(plist) <- ord_meths
        
        pdataframe = ldply(plist, function(x){
          df = x$data[, 1:2]
          colnames(df) = c("Axis_1", "Axis_2")
          return(cbind(df, x$data))
        })
        names(pdataframe)[1] = "method"
        
        pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)
        
        pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)
        
        
        # Plots
        p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
        
        
        if (!is.null(shapes)) {
          
          pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
          p = p +  geom_point(aes_string(shape=shapes), size=1.5, alpha = 3) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
          
        }else{
          p = p + geom_point(size=1.5, alpha = 3)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
        }
        
        p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
        p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
        p = p + theme(legend.position = "bottom", 
                      legend.title = element_blank(),
                      legend.justification="left", 
                      legend.box = "vertical",
                      legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                      plot.title=element_text(size=9,face="bold"))
        
        if(!is.null(mycols)){
          p <- p + scale_color_manual(values = mycols)
        }else{
          p <- p
        }
        
        # ID variation
        if (!is.null(ID)) {
          p = p + geom_text_repel(aes_string(label = ID), size = 2)
        } else {
          p = p 
        }
        
        # ellipse variation
        if (ellipse == "yes" | ellipse == "Yes" ) {
          p = p + stat_ellipse(type = "norm", linetype = 2) 
        } else if (ellipse == "no" | ellipse == "No" ){
          p = p 
        }
        
        
        
        
        
        if (!is.null(facet)) {
          ncol <- length(unique(mapping.sel.na.rem[,facet]))
          p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
        }
        else {
          p = p
        }
        
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  dev.off()
}
