

Go_pheatmap <- function(psIN,project, title, 
                        group1=NULL, group2=NULL, group3=NULL, group4=NULL,
                        Ntax=NULL, 
                        name=NULL, 
                        col_orders=NULL,
                        show_rownames = T,show_colnames = F,
                        cutree_rows = NA, cutree_cols = NA,
                        cluster_rows = T, cluster_cols = T, 
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
  # Define a function to assign colors to Phylum/Path
  phylumcolor <- c("#E7695DFF", "#6B8993FF", "#F6F0D4FF", "#95CE8AFF", "#D2D2D2FF", "#94D4D4FF", "#969696FF", "#F1F3E8FF", "#88775FFF")
  
  assign_colors <- function(x, colors) {
    getPalette = colorRampPalette(colors)
    col <- getPalette(length(unique(x)))
    names(col) = levels(x)
    return(col)
  }
  
  # Based on the type, generate annotation row and assign colors
  if(type %in% c("taxonomy", "taxanomy")) {
    annotation_row <- data.frame(Phylum = as.factor(tax_table(ps.rel.sel)[, "Phylum"]))
    phylum_col <- assign_colors(annotation_row$Phylum, phylumcolor)
  } else if(type %in% c("kegg", "pathway")) {
    annotation_row <- data.frame(Path = as.factor(tax_table(ps.rel.sel)[, "KO" | "pathway"]))
    Path_col <- assign_colors(annotation_row$Path, Pathcolor)  # Ensure Pathcolor is defined somewhere above
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
  # Function to generate colors
  generate_colors <- function(group, hardcode_colors=NULL) {
    unique_vals <- unique(mapping.sel[, group])
    if(is.null(hardcode_colors)) {
      stop(paste("Group", group, "not found in hardcoded_colors"))
    }
    color <- head(hardcode_colors, length(unique_vals))
    names(color) <- levels(as.factor(unique_vals))
    return(color)
  }
  
  # Function to create annotation_col data frame
  generate_annotation_col <- function(groups) {
    annotation_col <- as.data.frame(lapply(groups, function(x) as.factor(mapping.sel[,x])))
    colnames(annotation_col) <- groups
    return(annotation_col)
  }
  
  # Larger list of hardcoded colors
  all_hardcoded_colors <- list(
    color_set_1 = c("#EB5291FF", "#FBBB68FF", "#F5BACFFF", "#9DDAF5FF", "#6351A0FF","#ECF1F4FF", "#FEF79EFF", "#1794CEFF","#972C8DFF"),
    
    color_set_2 = c("#B15928", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3", "#1170aa", "#fc7d0b", "#76B7B2", "#B07AA1", "#E15759", "#59A14F", "#EDC948", "#FF9DA7", "#9C755F","#BAB0AC","#C84248"),
    color_set_3 = c("#2366C0FF", "#E9D738FF", "#B91226FF", "#A3DA4BFF", "#FF6435FF"),
    color_set_4 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  )
  
  # Define groups based on the inputs
  groups <- c(group1, group2, group3, group4)  # Ensure these are the correct variable names for your groups
  groups <- groups[!is.null(groups)]
  
  # Define hardcoded_colors for the current groups
  hardcoded_colors <- all_hardcoded_colors[1:length(groups)]
  names(hardcoded_colors) <- groups
  
  # Generate group colors
  group_colors <- lapply(groups, function(x) generate_colors(x, hardcoded_colors[[x]]))
  
  # Generate annotation_col
  annotation_col = generate_annotation_col(groups)
  
  # Generate annotation_colors
  ann_colors <- c(group_colors, if(type %in% c("taxonomy", "taxanomy")) list(Phylum = phylum_col) else list(Path = Path_col))
  names(ann_colors) <- c(groups, if(type %in% c("taxonomy", "taxanomy")) "Phylum" else "Path")
  
  
  
  #===== data process
  
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

