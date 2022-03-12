#' A Go_bdiv
#'


Go_bdivOrdi <- function(psIN, project, metaData, plot, distance_metrics, name=NULL){
    

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_ordi <- file.path(sprintf("%s_%s/table/odination",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_ordi)) dir.create(out_ordi)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  # create ordinate list
  for (mvar in rownames(subset(metadata, Go_bdiv =="yes"))) {
    for(distance_metric in distance_metrics){
      # remove na
      mapping.sel <- data.frame(sample_data(psIN))
      mapping.sel[mapping.sel==""] <- "NA"
      mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
      na.count <- length(mapping.sel.na)
      psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
      mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
      
      
      
      ## fix factor  and  numeric
      if (metadata[mvar, "type"] == "factor") {
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
        sample_data(psIN.na) <- mapping.sel.na.rem
      } else if (metadata[mvar, "type"] == "numeric") {
        next
      }
      
      
      
      ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
      plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
      ordi = ordinate(psIN.na, method=i, distance=distance_metric)
      saveRDS(ordi,sprintf("%s/ordi.%s.%s.%s.%s.%s%s.rds",out_ordi,
                           project,
                           i, 
                           distance_metric,
                           mvar, 
                           ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                           format(Sys.Date(), "%y%m%d")))
      }, psIN.na, distance_metric)
    }
  }
}
