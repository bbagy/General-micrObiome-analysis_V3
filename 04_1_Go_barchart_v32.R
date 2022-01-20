#' A Go_barchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_barchart()


Go_barchart <- function(psIN, metaData, project, taxanames, data_type, simple = "no",  
                        x_label="SampleIDfactor", facet=NULL, legend="bottom", orders,
                        cutoff=0.005, name=NULL, ncol=11, height, width,plotCols,  plotRows){
    
  if(!is.null(dev.list())) dev.off()
  
  colorset = "Set1" # Dark1 Set1 Paired
  #taxRanks <- c("Phylum","Class","Order","Family", "Genus","Species")
  
  taxRanks <- taxanames
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_taxa <- file.path(sprintf("%s_%s/table/taxa",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_taxa)) dir.create(out_taxa)

  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

  
  
  
  # logic for out file
  if (!is.null(facet)) {
    if (!is.null(name)) {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet,name, cutoff, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
    else {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project, facet, cutoff, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (!is.null(name)) {
      pdf(sprintf("%s_%s/pdf/3_barchart_simple.%s.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,cutoff, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }      else {
      pdf(sprintf("%s_%s/pdf/3_barchart.%s.(%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,cutoff, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  

  
  # order by bdiv

  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  for(i in 1:length(taxanames)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxanames[i])
    
    if (dim(otu.filt)[2] == 2){
      next
    }
    
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames[i])), otu.filt, sum, na.action=na.pass)
    
    if (taxanames[i] == "Species"){
      agg <- agg[grepl("NA NA", agg$Species)==F,]
     
    }
    
    genera <- agg[,taxanames[i]]
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    inds_to_grey <- which(rowMeans(agg) < cutoff)
    genera[inds_to_grey] <- "[1_#Other]"
    agg[,taxanames[i]] <- genera
    #saving table
    agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
    write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_abundance.(%s).%s.%s.csv", out_taxa, project, cutoff,taxanames[i], format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
    
    df <- melt(agg, variable="SampleID")


    # add StduyID
    
    
    df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , taxanames[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)

    #mapping.sel[df2$SampleID, "StudyID"]
   
    # add groups
    for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
      df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, mvar])
      df2[,mvar] <- mapping.sel[df2$SampleID, mvar]

      # order
      if (length(orders) >= 1) {
        df2[,mvar] <- factor(df2[,mvar], levels = orders)
      }
      else {
        df2[,mvar] <- factor(df2[,mvar])
      }
    }
    
    # adding facet to groups
    if (!is.null(facet)) {
      for (fa in facet){
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
        df2[,fa] <- factor(df2[,fa], levels = orders)
      }
    }
   
   
    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    } 
    


    print(1)
    # color
    colourCount = length(unique(df2[,taxanames[i]]));colourCount
    getPalette = colorRampPalette(brewer.pal(9, colorset))

    # pdf size height = 5, width=9
   
    if (legend == "bottom"){
      if (colourCount < 30) {
        coln <- 4
      }else if (colourCount > 30) {
        coln <- 5
      }
    } else if (legend == "right") {
      if (colourCount < 18) {
        coln <- 1
      } else if (colourCount > 19 & colourCount  < 35) {
        coln <- 2
      } else if (colourCount > 36) {
        coln <- 3
      }
    }

    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    print(2)

    p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxanames[i], order=taxanames[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill=NULL)+
      theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + 
      guides(fill=guide_legend(ncol= coln))  + #guides(col = guide_legend(ncol = coln)) + 
      ylim(c(-.1, 1.01)) + scale_fill_manual(values = getPalette(colourCount)) + labs(y = "Relative abundance")
    

    if (!is.null(facet)) {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
        }
        
        if (facet == mvar) {
          next
        }

        df2[,facet] <- factor(df2[,facet], levels = orders)
        
        print(sprintf("Facet by %s-%s",mvar, facet))
        p <- p+ facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", ncol = ncol) 

        if (!is.null(name)) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }

        print(p)
      }

    }     else if (is.null(facet) & simple == "no") {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }
        print("B")
        print(sprintf("Facet by %s",mvar))
        p <- p+  facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales = "free_x", ncol = ncol)  
        
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i], mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    } else if (is.null(facet) & simple == "yes") {
      for (mvar in rownames(subset(metadata, Go_barchart =="yes"))) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }
        print("C")
        print("Simpe plot")
        
        p = p
        
        if (!is.null(name)) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i],mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  dev.off()
}
