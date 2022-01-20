#' A Go_biplot
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity ordination plot
#' @export
#' @examples
#' Go_biplot()

Go_biplot <- function(psIN, metaData, project, orders, distance_metrics, data_type, biplot, shapes, TaxaLevel, ID, facet, des, name, height, width){
  
  if(!is.null(dev.list())) dev.off()
  
  colorset = "Dark2" # Dark1 Set1 Paired
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  cat("#--  Getting Taxa  --#")
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN))) # for dada2 
  }else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame((otu_table(psIN))) # not for dada2 
  }
  

  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN),taxRanks=ranks, level= TaxaLevel)
  agg <- aggregate(. ~ Genus, otu.filt, sum, na.action=na.pass) #add na.action=na.pass if have error "no rows to aggregate"
  genera <- agg$Genus
  agg <- agg[,-1]
  rownames(agg) <- genera
  
  agg_t <-t(agg)
  
  dim(agg)
  
  
  
  # out file
  if (length(des) == 1) {
    if (length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else if (length(facet) == 1 & length(name) == 1) {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,  des, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, des, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }


  # out file
    mapping.sel <- data.frame(sample_data(psIN))

    for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
        for(distance_metric in distance_metrics){
          #na remove
          mapping.sel <- data.frame(sample_data(psIN))
          mapping.sel[mapping.sel==""] <- "NA"
          mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
          na.count <- length(mapping.sel.na)
          psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
          mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))

          if (length(des) == 1) {
            print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                          des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          } else{
            print(sprintf("##-- %s (total without NA: %s/%s) --##",
                          mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
          }


          if (class(mapping.sel.na.rem[,mvar]) == "integer"){
            mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
            sample_data(psIN.na) <- mapping.sel.na.rem
          }

          
          
          ordi <- ordinate(psIN , method = "NMDS", distance = distance_metric)
          #Get site information
          df <- scores(ordi,display=c("sites"));head(df)
          
          # Add grouping information
          df <- merge(df, mapping.sel.na.rem, by="row.names")
          
          #Get the vectors for bioenv.fit
          bio.fit <- envfit(ordi, agg_t[,biplot,drop=F],  perm = 999); head(bio.fit)
          df_biofit <- scores(bio.fit, display=c("vectors"))
          df_biofit <- df_biofit*vegan:::ordiArrowMul(df_biofit)
          df_biofit <- as.data.frame(df_biofit);df_biofit
          

          if (length(shapes) == 1) {
            p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar, shape=shapes), size=2) + 
              stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                           inherit.aes = TRUE) 
          }
          else{
            p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar), size=2) +
              stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                           inherit.aes = TRUE) 
          }
          
          p = p + theme_bw() +  ggtitle(sprintf("%s - %s - %s - %s", project, mvar, "NMDS", distance_metric)) #+ geom_polygon()
          # open(1), cross(10), closed(2)
          p = p + scale_fill_brewer(type="qual", palette=colorset)#scale_colour_manual(values=colours) 
          p = p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                         arrow = arrow(length = unit(0.2, "cm")), color="#808080",alpha=0.7)+
            geom_text(data=as.data.frame(df_biofit*1.1), aes(NMDS1, NMDS2, label = rownames(df_biofit)), 
                      color="#808080",alpha=0.7, size = 3)


          if (length(ID) == 1) {
            p = p +  geom_text_repel(aes_string(label = ID), size = 2)
          }

          else {
            p = p 
          }

          if (length(orders) >= 1) {
            p = p + scale_colour_brewer(type="qual", palette=colorset, breaks=orders) + 
              scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) 
          }
          else {
            p = p + scale_colour_brewer(type="qual", palette=colorset) + 
              scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16)) 
          }
          
          if (length(facet) == 1) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }
          else {
            p = p
          }
          print(p)
        }
      }
    dev.off()
  } else {
      if (length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else if (length(facet) == 1 & length(name) == 1) {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project, facet, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      else {
        pdf(sprintf("%s_%s/pdf/12_biplot.%s.%s.pdf",project,format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
      }
      
      for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
          for(distance_metric in distance_metrics){
            #na remove
            mapping.sel <- data.frame(sample_data(psIN))
            mapping.sel[mapping.sel==""] <- "NA"
            mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
            na.count <- length(mapping.sel.na)
            psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
            mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))


            if (length(des) == 1) {
              print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                            des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            } else{
              print(sprintf("##-- %s (total without NA: %s/%s) --##",
                            mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
            }


            if (class(mapping.sel.na.rem[,mvar]) == "integer"){
              mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
              sample_data(psIN.na) <- mapping.sel.na.rem
            }

            ordi <- ordinate(psIN , method = "NMDS", distance = distance_metric)
            #Get site information
            df <- scores(ordi,display=c("sites"));head(df)
            
            # Add grouping information
            df <- merge(df, mapping.sel.na.rem, by="row.names")
            #Get the vectors for bioenv.fit
            bio.fit <- envfit(ordi, agg_t[,biplot,drop=F],  perm = 999); head(bio.fit)
            df_biofit <- scores(bio.fit, display=c("vectors"))
            df_biofit <- df_biofit*vegan:::ordiArrowMul(df_biofit)
            df_biofit <- as.data.frame(df_biofit);df_biofit
            
            
            if (length(shapes) == 1) {
              p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar, shape=shapes), size=2) + 
                stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                             inherit.aes = TRUE) 
            }
            else{
              p = ggplot()+ geom_point(data=df,aes_string("NMDS1","NMDS2", colour=mvar), size=2) +
                stat_ellipse(data=df, aes_string("NMDS1","NMDS2", colour=mvar), level=0.95, type = "t", linetype = 3, 
                             inherit.aes = TRUE) 
            }
            
            p = p + theme_bw() +  ggtitle(sprintf("%s - %s - %s - %s", project, mvar, "NMDS", distance_metric)) #+ geom_polygon()
            # open(1), cross(10), closed(2)
            p = p + scale_fill_brewer(type="qual", palette=colorset)#scale_colour_manual(values=colours) 
            p = p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                               arrow = arrow(length = unit(0.2, "cm")), color="#808080",alpha=0.7)+
              geom_text(data=as.data.frame(df_biofit*1.1), aes(NMDS1, NMDS2, label = rownames(df_biofit)), 
                        color="#808080",alpha=0.7, size = 3)
            
            
            if (length(ID) == 1) {
              p = p +  geom_text_repel(aes_string(label = ID), size = 2)
            }
            
            else {
              p = p 
            }
            
            if (length(orders) >= 1) {
              p = p + scale_colour_brewer(type="qual", palette=colorset, breaks=orders) + 
                scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) 
            }
            else {
              p = p + scale_colour_brewer(type="qual", palette=colorset) + 
                scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16)) 
            }
            
            if (length(facet) == 1) {
              ncol <- length(unique(mapping.sel.na.rem[,facet]))
              p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
            }
            else {
              p = p
            }
            print(p)
          }
      }
      dev.off()
  }
}

