

Go_pheatmap <- function(psIN,project, title, group1=NULL, group2=NULL,Ntax=NULL, name=NULL,
                        show_rownames = T,show_colnames = F,
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
  # ps.rel <- transform_sample_counts(psIN.prune, function(x) x/sum(x)*100);ps.rel
  # option 2
   ps.rel <- transform_sample_counts(psIN.prune, function(x) log2(x/sum(x)*100));ps.rel

  
  
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
  
  taxaTab <- data.frame(tax_table(ps.rel.sel)[,"Species"])
  
  
  
  # map 정리
  mapping <- data.frame(sample_data(ps.rel.sel));dim(mapping)
  sel <- intersect(rownames(mapping), colnames(matrix)); head(sel, "3")
  mapping.sel <- mapping[sel,, drop=F];dim(mapping.sel)
  
  
  
  # phylum annotation
  annotation_row = data.frame(
    Phylum = as.factor(tax_table(ps.rel.sel)[, "Phylum"])
  )
  rownames(annotation_row) = rownames(matrix)
  
  # phylum colors
  # phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Dark2")
  
  phylum_col <- head(brewer.pal(8, "Dark2"),length(unique(annotation_row$Phylum)))
  names(phylum_col) = levels(annotation_row$Phylum)
  
  # add group(s) and color list
  if (!is.null(group2)){
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]), 
      group2 = as.factor(mapping.sel[,group2]), 
      check.names = FALSE
    )
    colnames(annotation_col) <-c(group1, group2)
    
    # group colors
    group1.col <- head(brewer.pal(8, "Set2"),length(unique(mapping.sel[,group1])))
    names(group1.col) = levels(as.factor(mapping.sel[,group1]))
    
    group2.col <- head(brewer.pal(12, "Paired"),length(unique(mapping.sel[,group2])))
    names(group2.col) <- unique(mapping.sel[,group2])
    # color list
    ann_colors = list(
      group1 = group1.col,
      group2 = group2.col,
      Phylum = phylum_col
    )
    names(ann_colors) <-c(group1, group2, "Phylum")
  }else{
    annotation_col = data.frame(
      group1 = as.factor(mapping.sel[,group1]),
      check.names = FALSE
    )
    colnames(annotation_col) <-c(group1)
    
    # group colors
    group1.col <- head(brewer.pal(8, "Set2"),length(unique(mapping.sel[,group1])))
    names(group1.col) = levels(as.factor(mapping.sel[,group1]))
    
    # color list
    ann_colors = list(
      group1 = group1.col,
      Phylum = phylum_col
    )
    names(ann_colors) <-c(group1, "Phylum")
    
  };ann_colors
  
  rownames(annotation_col) = rownames(mapping.sel)
  
  p <- ComplexHeatmap::pheatmap(matrix, scale= "row", 
                                main = title,
                                annotation_col = annotation_col, 
                                annotation_row = annotation_row, 
                                show_rownames = show_rownames,
                                show_colnames = show_colnames,
                                cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                labels_row=taxaTab$Species,
                                annotation_colors = ann_colors)
  
  
  if (!is.null(name)) {
    pdf(sprintf("%s/pheatmap.%s.(%s).%s.%s.pdf", out_path, project,title, name,format(Sys.Date(), "%y%m%d")), height = h, width = width)
  }   else {
    pdf(sprintf("%s/pheatmap.%s.(%s).%s.pdf", out_path, project,title, format(Sys.Date(), "%y%m%d")), height = h, width = width)
  }
  
  print(p)
  
  dev.off()
}

