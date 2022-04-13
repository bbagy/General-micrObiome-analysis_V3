

Go_pheatmap <- function(psIN,project, title, group1=NULL, group2=NULL,group3=NULL,Ntax=NULL, name=NULL,
                        show_rownames = T,show_colnames = F,type,showPhylum = T,
                        cutree_rows = NA, cutree_cols = NA,
                        cluster_rows = T,cluster_cols = T, 
                        width){
  # BiocManager::install("ComplexHeatmap")
  # install.packages("Cairo")
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  Tableau10 = c("#B07AA1","#FF9DA7","#76B7B2", "#59A14F","#EDC948", "#9C755F", "#BABOAC","#1170aa", "#fc7d0b","#E15759") 
  
  psIN.prune = prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
  
  
  #----- normalization relative abundant ---#
  # option 1
   ps.rel <- transform_sample_counts(psIN.prune, function(x) x/sum(x)*100);ps.rel
  # option 2
   # ps.rel <- transform_sample_counts(psIN.prune, function(x) log2(x));ps.rel

   
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
  if(type == "taxonomy" | type == "taxanomy" ){
    taxaTab <- data.frame(tax_table(ps.rel.sel)[,"Species"])
  }else if(type == "function"){
    taxaTab <- data.frame(tax_table(ps.rel.sel)[,"Path.des"])
  }
  


  
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
  annotation_row = data.frame(
    if(type == "taxonomy" | type == "taxanomy" ){
      Phylum = as.factor(tax_table(ps.rel.sel)[, "Phylum"])
    }else if(type == "function"){
      Path = as.factor(tax_table(ps.rel.sel)[, "Path"])
    }
  )
  
  
  if(type == "taxonomy" | type == "taxanomy" ){
    colnames(annotation_row) <- c("Phylum")
    
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    phylum_col <- getPalette(length(unique(annotation_row$Phylum)))
    
    #phylum_col <- head(brewer.pal(8, "Dark2"),length(unique(annotation_row$Phylum)))
    names(phylum_col) = levels(annotation_row$Phylum)
    
  } else if(type == "function"){
    colnames(annotation_row) <- c("Path")
    
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
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
      group2 = as.factor(mapping.sel[,group2]), 
    ) 
    
    colnames(annotation_col) <-c(group1, group2)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(8, "Set2"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(1)
    } else{
      group1.col <- head(brewer.pal(8, "Set2"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(2)
    }
    
    
    
    
    group2.col <- head(brewer.pal(8, "Accent"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    
    # color list
    if(type == "taxonomy" | type == "taxanomy" ){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        Phylum = phylum_col
      )
      names(ann_colors) <-c(group1, group2, "Phylum")
    }else if(type == "function"){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, group2,"Path")
    }
    
  } else if (!is.null(group2) & !is.null(group3)){
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]), 
      group2 = as.factor(mapping.sel[,group2]), 
      group3 = as.factor(mapping.sel[,group3])
    )
    
    colnames(annotation_col) <-c(group1, group2,group3)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(8, "Set2"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(1)
    } else{
      group1.col <- head(brewer.pal(8, "Set2"),length(unique(mapping.sel[,group1])))
      names(group1.col) = levels(as.factor(mapping.sel[,group1]))
      print(2)
    }
    
    group2.col <- head(brewer.pal(8, "Accent"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    
    group3.col <- head(brewer.pal(8, "Dark2"),length(unique(mapping.sel[,group3])))
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
    }else if(type == "function"){
      ann_colors = list(
        group1 = group1.col,
        group2 = group2.col,
        group3 = group3.col,
        Path = Path_col
      )
      names(ann_colors) <-c(group1, group2,group3,"Path")
    }
  } else{
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]),
      check.names = FALSE
    )
    colnames(annotation_col) <- c(group1)
    
    # group colors
    
    if (length(unique(mapping.sel[,group1])) > 8){
      getPalette = colorRampPalette(brewer.pal(8, "Set2"))
      group1.col <- getPalette(length(unique(mapping.sel[,group1])))
      names(Path_col) = levels(annotation_row$Path)
      print(3)
    } else{
      group1.col <- head(brewer.pal(8, "Set2"),length(unique(mapping.sel[,group1])))
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
    }else if(type == "function"){
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
    matrix<- t(matrix)
  }
  
  colSums(matrix)
  
  bk <- c(0,0.5,1)

  if (showPhylum ==TRUE){
    p <- ComplexHeatmap::pheatmap(matrix,  fontsize =8,main = title,
                                  #scale= "row",
                                  annotation_col = annotation_col, 
                                  annotation_row = annotation_row, 
                                  show_rownames = show_rownames,
                                  show_colnames = show_colnames,
                                  cluster_rows = cluster_rows,
                                  cluster_cols = cluster_cols,
                                  labels_row=taxaTab$Species,
                                  cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                  #color=c("seashell1", "seashell2", "seashell3"),
                                  #breaks= bk,
                                  #legend_breaks= bk,
                                  annotation_colors = ann_colors)
  } else{
    p <- ComplexHeatmap::pheatmap(matrix,  fontsize =8,main = title,
                                  #scale= "row",
                                  annotation_col = annotation_col, 
                                  #annotation_row = annotation_row, 
                                  show_rownames = show_rownames,
                                  show_colnames = show_colnames,
                                  cluster_rows = cluster_rows,
                                  cluster_cols = cluster_cols,
                                  labels_row=taxaTab$Species,
                                  cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                  #color=c("seashell1", "seashell2", "seashell3"),
                                  #breaks= bk,
                                  #legend_breaks= bk,
                                  annotation_colors = ann_colors)
  }

  

  # logic for out file
  pdf(sprintf("%s//pheatmap.%s.(%s).%s%s.pdf", out_path, 
              project, 
              title,
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = h, width = width)
  
  
  print(p)
  
  dev.off()
}

