

Go_pheatmap <- function(psIN,project, title, 
                        group1=NULL, group2=NULL, group3=NULL, group4=NULL,
                        Ntax=NULL, 
                        name=NULL, 
                        col_orders=NULL,
                        show_rownames = T,show_colnames = F,
                        cutree_rows = NA, cutree_cols = NA,
                        cluster_rows = T,cluster_cols = T, 
                        showPhylum = T,
                        width){
  # BiocManager::install("ComplexHeatmap")
  # install.packages("Cairo")
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_pdf <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_pdf)) dir.create(out_pdf)
  out_tab <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_tab)) dir.create(out_tab)
  out_pheatmapTab <- file.path(sprintf("%s_%s/table/pheatmapTab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_pheatmapTab)) dir.create(out_pheatmapTab)
  

  #----- normalization relative abundant ---#
   if(max(data.frame(otu_table(psIN))) < 1){
     ps.rel <- psIN
     print("Percentage abundant table data.")
   }else{
     print("Counts or cpm table data. Data normalized by perentage.")
     psIN.prune = prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
     ps.rel <- transform_sample_counts(psIN.prune, function(x) x/sum(x)*100);ps.rel# log(1+x) 를 하면 NaN가 많이 나온다. 
   }
   
print("Check the psIN")
   
  if(is.null(Ntax)){
    Ntax <- dim(tax_table(ps.rel))[1]
    print(sprintf("number of taxa is %s",Ntax))
    ps.rel.sel <- prune_taxa(names(sort(taxa_sums(ps.rel),TRUE)[1:Ntax]), ps.rel);ps.rel.sel
  }else{
    ps.rel.sel <- prune_taxa(names(sort(taxa_sums(ps.rel),TRUE)[1:Ntax]), ps.rel);ps.rel.sel 
  }

   
  # for height
   
   if ( dim(tax_table(ps.rel.sel))[1] > 30){
     h=8
   }else if(dim(tax_table(ps.rel.sel))[1] <= 30 & dim(tax_table(ps.rel.sel))[1] > 25){
     h=7.4
   }else if(dim(tax_table(ps.rel.sel))[1] <= 25 & dim(tax_table(ps.rel.sel))[1] > 20){
     h=6.3
   }else if(dim(tax_table(ps.rel.sel))[1] <= 20 & dim(tax_table(ps.rel.sel))[1] > 15){
     h=5
   }else if(dim(tax_table(ps.rel.sel))[1] <= 15 & dim(tax_table(ps.rel.sel))[1] > 10){
     h=5
   }else if(dim(tax_table(ps.rel.sel))[1] <= 10){
     h=4.5
   }
   
  
  matrix <- data.frame(t(otu_table(ps.rel.sel)))
  

  
  # normalization for log2
  is.na(matrix)<-sapply(matrix, is.infinite)
  matrix[is.na(matrix)]<-0
  
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  colnames(matrix) <- gsub("X","",colnames(matrix))
  
  

  
 # get taxa names
  print("Check the data type")
  taxtab.col <- colnames(data.frame((tax_table(ps.rel.sel))))
  
  if (any(grepl("Species", taxtab.col))){
    taxaTab <- data.frame(tax_table(ps.rel.sel)[,"Species"])
    type <- "taxonomy"
    print(type)
  }else if(any(grepl("KO", taxtab.col))){
    taxaTab <- data.frame(tax_table(ps.rel.sel)[,"KO"])
    type <- "kegg"
    print(type)
  }else if(any(grepl("pathway", taxtab.col))){
    taxaTab <- data.frame(tax_table(ps.rel.sel)[,"pathway"])
    type <- "pathway"
    print(type)
    }
  
  colnames(taxaTab) <- "Rank"
  taxaTab.print <- data.frame(tax_table(ps.rel.sel))
  
  write.csv(taxaTab.print,file=sprintf("%s/pheatmap.%s.tap(%s).%s%s%s.csv",out_pheatmapTab,
                                       project,
                                       Ntax, 
                                       ifelse(is.null(title), "", paste(title, ".", sep = "")), 
                                       ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                       format(Sys.Date(), "%y%m%d"),
                                       sep="/"), quote = FALSE, col.names = NA)


  # map 정리
  mapping <- data.frame(sample_data(ps.rel.sel));dim(mapping)
  
  tt <- try(sel <- intersect(rownames(mapping), colnames(matrix)), T)
  
  if (length(tt) == 0){
    sel <- intersect(rownames(mapping), rownames(matrix)); head(sel, "3")
  }else{
    sel <- intersect(rownames(mapping), colnames(matrix)); head(sel, "3")
  }

  mapping.sel <- mapping[sel,, drop=F];dim(mapping.sel)
  
  
  
  # phylum annotation
  print("Get_annotation")
  annotation_row = data.frame(
    if(type == "taxonomy" | type == "taxanomy" ){
      Phylum = as.factor(tax_table(ps.rel.sel)[, "Phylum"])
    }else if(type == "kegg"){
      Path = as.factor(tax_table(ps.rel.sel)[, "KO"])
    }else if(type == "pathway"){
      Path = as.factor(tax_table(ps.rel.sel)[, "pathway"])
    }
  )
  
  
  if(type == "taxonomy" | type == "taxanomy" ){
    colnames(annotation_row) <- c("Phylum")
    
    getPalette = colorRampPalette(brewer.pal(8, "Paired"))
    phylum_col <- getPalette(length(unique(annotation_row$Phylum)))
    
    #phylum_col <- head(brewer.pal(8, "Dark2"),length(unique(annotation_row$Phylum)))
    names(phylum_col) = levels(annotation_row$Phylum)
    
  } else if(type == "kegg" | type == "pathway"){
    colnames(annotation_row) <- c("Path")
    
    getPalette = colorRampPalette(brewer.pal(8, "Paired"))
    Path_col <- getPalette(length(unique(annotation_row$Path)))
    names(Path_col) = levels(annotation_row$Path)
  }
  


  tt <- try(rownames(annotation_row) <- rownames(matrix), T)
  
  if (class(tt) == "try-error"){
    rownames(annotation_row) <- colnames(matrix)
  }else{
    rownames(annotation_row) <- rownames(matrix)
  }
  
  
  
  # phylum colors
  # phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Dark2")

  
  # add group(s) and color list
  if (!is.null(group2) & is.null(group3)){
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]), 
      group2 = as.factor(mapping.sel[,group2])
    ) 
    
    colnames(annotation_col) <-c(group1, group2)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(1)
    } else{
      group1.col <- head(brewer.pal(9, "Set1"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(2)
    }
    
    
    group2.col <- head(brewer.pal(12, "Paired"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    
    # color list
    if(type == "taxonomy" | type == "taxanomy" ){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        Phylum = phylum_col
      )
      names(ann_colors) <-c(group1, group2, "Phylum")
    }else if(type == "kegg" | type == "pathway"){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, group2,"Path")
    }
    
  } else if (!is.null(group2) & !is.null(group3) & is.null(group4)){
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]), 
      group2 = as.factor(mapping.sel[,group2]), 
      group3 = as.factor(mapping.sel[,group3])
    )
    
    colnames(annotation_col) <-c(group1, group2, group3)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(1)
    } else{
      group1.col <- head(brewer.pal(9, "Set1"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(2)
    }
    
    group2.col <- head(brewer.pal(12, "Paired"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    
    col3 <- c("#B15928", "#FFFF99", "#6A3D9A", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3", "#1170aa", "#fc7d0b", "#76B7B2", "#B07AA1", "#E15759", "#59A14F", "#EDC948", "#FF9DA7", "#9C755F","#BAB0AC","#C84248")
    group3.col <- head(col3,length(unique(mapping.sel[,group3])))
    names(group3.col) <- unique(mapping.sel[,group3])
    
    # color list
    if(type == "taxonomy" | type == "taxanomy" ){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        group3 = group3.col,
        Phylum = phylum_col
      )
      names(ann_colors) <-c(group1, group2,group3, "Phylum")
    }else if(type == "kegg" | type == "pathway"){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        group3 = group3.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, group2,group3,"Path")
    }
  } else if (!is.null(group2) & !is.null(group3) & !is.null(group4)){
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]), 
      group2 = as.factor(mapping.sel[,group2]), 
      group3 = as.factor(mapping.sel[,group3]),
      group4 = as.factor(mapping.sel[,group4])
    )
    
    colnames(annotation_col) <-c(group1, group2, group3,group4)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(1)
    } else{
      group1.col <- head(brewer.pal(9, "Set1"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(2)
    }
    
    group2.col <- head(brewer.pal(12, "Paired"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    
    col3 <- c("#B15928", "#FFFF99", "#6A3D9A", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3", "#1170aa", "#fc7d0b", "#76B7B2", "#B07AA1", "#E15759", "#59A14F", "#EDC948", "#FF9DA7", "#9C755F","#BAB0AC","#C84248")
    group3.col <- head(col3,length(unique(mapping.sel[,group3])))
    names(group3.col) <- unique(mapping.sel[,group3])
    
    group4.col <- head(rev(brewer.pal(8, "Dark2")),length(unique(mapping.sel[,group4])))
    names(group4.col) <- unique(mapping.sel[,group4])
    
    # color list
    if(type == "taxonomy" | type == "taxanomy" ){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        group3 = group3.col,
        group4 = group4.col,
        Phylum = phylum_col
      )
      names(ann_colors) <-c(group1, group2,group3,group4, "Phylum")
    }else if(type == "kegg" | type == "pathway"){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        group3 = group3.col,
        group4 = group4.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, group2,group3,group4,"Path")
    }
  } else{
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]),
      check.names = FALSE
    )
    colnames(annotation_col) <- c(group1)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(3)
    } else{
      group1.col <- head(brewer.pal(9, "Set1"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(4)
    }
    
    # color list
    
    if(type == "taxonomy" | type == "taxanomy" ){
      ann_colors = list(
        group1 = group1.col,
        Phylum = phylum_col
      )
      names(ann_colors) <-c(group1, "Phylum")
    }else if(type == "kegg" | type == "pathway"){
      ann_colors = list(
        group1 = group1.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, "Path")
    }
  }
 
  
  
  rownames(annotation_col) = rownames(mapping.sel)
  
  if(type == "taxonomy" | type == "taxanomy" ){
    matrix = as.matrix(matrix)
  }else if(type == "function"){
    matrix <- t(matrix)
  }
  
  colSums(matrix)
  
  bk <- c(0,0.5,1)
  print("p0")
  
  tt<-try(ComplexHeatmap::pheatmap(matrix, annotation_col = annotation_col),T)
  if (class(tt) == "try-error"){
    matrix <- t(matrix)
  }else{
    matrix <- matrix
  }

  if (showPhylum ==TRUE){
    print("with annotation_row")
    p <- ComplexHeatmap::pheatmap(matrix,  fontsize =8,main = title,
                                  #scale= "row",
                                  annotation_col = annotation_col, 
                                  annotation_row = annotation_row, 
                                  show_rownames = show_rownames,
                                  show_colnames = show_colnames,
                                  cluster_rows = cluster_rows,
                                  cluster_cols = cluster_cols,
                                  labels_row=taxaTab$Rank,
                                  cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                  #color=c("seashell1", "seashell2", "seashell3"),
                                  #breaks= bk,
                                  #legend_breaks= bk,
                                  annotation_colors = ann_colors)

  } else{
    print("without annotation_row")

    if(!is.null(col_orders)){
      order <- col_orders[col_orders %in% colnames(matrix)] # match order to matrix column names
      matrix.orderd <- matrix[, order]
    }else{
      matrix.orderd <- matrix
    }
    
    p <- ComplexHeatmap::pheatmap(matrix.orderd,  fontsize =8, main = title,
                                  annotation_col = annotation_col, 
                                  show_rownames = show_rownames,
                                  show_colnames = show_colnames,
                                  cluster_rows = cluster_rows,
                                  cluster_cols = cluster_cols,
                                  labels_row=taxaTab$Rank,
                                  cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                  annotation_colors = ann_colors)
    
    

      print("p3")
  }


  # logic for out file
  pdf(sprintf("%s/pheatmap.%s.%s%s%spdf", out_pdf, 
              project, 
              ifelse(is.null(col_orders), "", paste("ordered", ".", sep = "")), 
              ifelse(is.null(title), "", paste(title, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = h, width = width)
  
  print("p4")
  print(p)
  
  dev.off()
}

