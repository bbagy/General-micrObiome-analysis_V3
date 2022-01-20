#' A Go_bdiv
#'


Go_bdivcbn <- function(psIN, metaData, project, orders, distance_metrics,
                       plot="PCoA", shapes = NULL, ID = NULL, ellipse="yes", facet=NULL, name=NULL,
                       combination, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  #colorset = "Dark2" # Dark1 Set1 Paired
  Tableau10 = c("#1170aa", "#fc7d0b",  "#76B7B2", "#E15759","#59A14F","#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BABOAC") 
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))



  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
    mapping <- data.frame(sample_data(psIN))
    
    if (metadata[mvar, "type"] == "factor") {
      mapping[,mvar] <- factor(mapping[,mvar])
    } else if (metadata[mvar, "type"] == "numeric") {
      mapping[,mvar] <- as.numeric(as.character(mapping[,mvar]))
    }
    
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
    
    #-----------------------#
    # for group combination #
    #-----------------------#
    mapping[,mvar] <- factor(mapping[,mvar], levels = orders)
    mapping[,mvar] <- factor(mapping[,mvar])
    
    group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
    
    print(count(group.cbn))
    
    
    
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
        if (metadata[mvar, "type"] == "factor") {
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
          sample_data(psIN.cbn.na) <- mapping.sel.na.rem
        } else if (metadata[mvar, "type"] == "numeric") {
          next
        }
        
        
        
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
        p = p + scale_color_manual(values = Tableau10)
        p = p + theme(legend.position = "bottom", 
                      legend.title = element_blank(),
                      legend.justification="left", 
                      legend.box = "vertical",
                      legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))
        
        
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
        
        plotlist[[length(plotlist)+1]] <- p
        
      }
    }
    
    # out file
  }
  if (!is.null(name)) {
    if (!is.null(facet)) {
      pdf(sprintf("%s_%s/pdf/ordi.%s.%s.%s.(cbn=%s).%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, name, combination,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      
    } else{
      pdf(sprintf("%s_%s/pdf/ordi.%s.%s.(cbn=%s).%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,name,combination,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }else{
    if (!is.null(facet)) {
      pdf(sprintf("%s_%s/pdf/ordi.%s.%s.(cbn=%s).%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, combination,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }else{
      pdf(sprintf("%s_%s/pdf/ordi.%s.(cbn=%s).%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,combination,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  multiplot(plotlist=plotlist)
  dev.off()
}
