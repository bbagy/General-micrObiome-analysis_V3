#' A Go_bdiv
#'


Go_bdivPlots <- function(psIN, metaData, project, orders, ordination, shapes = NULL, ID = NULL, ellipse="yes", facet=NULL, name=NULL, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

  ordi <- readRDS(ordination)
 # out file

  
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
    if (metadata[mvar, "type"] == "factor") {
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    } else if (metadata[mvar, "type"] == "numeric") {
      next
    }
    
    
    # extract a certain word fto finding plot_method and distance_metric
    words <- gsub(project,"", gsub("ordi","",ordination));words
    plot_method <- word(words,3,sep = fixed("."));print(plot_method)
    distance_metric <- word(words,4,sep = fixed("."));print(distance_metric)
    
    # plot
    # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
    plist = llply(as.list(plot_method), function(i, psIN.na, distance_metric){
      plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
    }, psIN.na, distance_metric)
    
    names(plist) <- plot_method
    
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
    p = p + scale_color_manual(values = mycols)
    p = p + theme(legend.position = "bottom", 
                  legend.title = element_blank(),
                  legend.justification="left", 
                  legend.box = "vertical",
                  legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                  plot.title=element_text(size=9,face="bold"))
    
    
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
    
    # facet
    if (!is.null(facet)) {
      ncol <- length(unique(mapping.sel.na.rem[,facet]))
      p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
    } else {
      p = p
    }
    
    #plotlist[[length(plotlist)+1]] <- p
    


    pdf(sprintf("%s/ordiPlot.%s.%s.%s.%s%s%s.pdf", out_path, 
                project, 
                plot_method,
                distance_metric,
                ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
                ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                format(Sys.Date(), "%y%m%d")), height = height, width = width)
    
    print(p)
    
    
  }
  dev.off()
}
