#' A Go_overview
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_overview()


Go_overview <- function(psIN, metaData, ylabn = "", facet, Color, orders, name,xanlgle, height, width){
    
  if(!is.null(dev.list())) dev.off()
    
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

  
  
  if (length(name) == 1) {
    pdf(sprintf("%s/3_taxa_overview.%s.%s.%s.%s.pdf", out_path, facet, project, name, format(Sys.Date(), "%y%m%d")), height = height, width=width)
  } else{
    pdf(sprintf("%s/3_taxa_overview.%s.%s.%s.pdf", out_path, facet, project ,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }
  

  
  
  print("Calculating relative abundance .....")
  ps3ra = transform_sample_counts(psIN, function(x){x / sum(x)})
  mphyseq <- psmelt(ps3ra)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  for (maingroup in rownames(subset(metadata, Go_overview =="yes"))) {
    mphyseq[,maingroup] <- as.character(mphyseq[,maingroup]);mphyseq[,maingroup]
    mphyseq[,maingroup][mphyseq[,maingroup]==""] <- "NA";mphyseq[,maingroup]
    mphyseq[,maingroup]<- as.factor(mphyseq[,maingroup]);mphyseq[,maingroup]
    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
    mphyseq.na <- subset(mphyseq, mphyseq[,maingroup] != "NA");mphyseq.na[,maingroup] 
    
    if (facet == "Genus") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "g__");mphyseq.na.na[,facet] 
    } else if (facet == "Family") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "f__");mphyseq.na.na[,facet] 
    } else if (facet == "Order") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "o__");mphyseq.na.na[,facet] 
    } else if (facet == "Class") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "c__");mphyseq.na.na[,facet] 
    } else if (facet == "Phylum") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "p__");mphyseq.na.na[,facet] 
    }else if (facet == "Species") {
      mphyseq.na.na <- subset(mphyseq.na, mphyseq.na[,facet] != "s__");mphyseq.na.na[,facet] 
    }
    
    
    if (length(orders) >= 1) {
      mphyseq.na.na[,maingroup] <- factor(mphyseq.na.na[,maingroup], levels = orders)
    }       else {
      mphyseq.na.na[,maingroup] <- factor(mphyseq.na.na[,maingroup])
    }
    
    p<- ggplot(data = mphyseq.na.na,  mapping = aes_string(x = maingroup, y = "Abundance",color = Color, fill = Color)) +
      geom_violin(fill = NA) + theme_bw()  +#scale_colour_brewer(type="qual", palette="Set4") + #+ 
      geom_point(size = 1, alpha = 0.3, position = position_jitter(width = 0.3)) +
      theme(title=element_text(size=8), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5)) +
      facet_wrap(facets = facet) + ylab(ylabn) + scale_y_log10()
    print(p)
  }
  
  dev.off()
}

