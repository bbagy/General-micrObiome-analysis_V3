##############################################
#----          install package         ------#
##############################################
packages <- c("ape", "car","cluster","CLME","compositions","cowplot","crayon", "caret","colorspace",
           "digest","data.table", "devtools","doParallel","ellipse", "emmeans","e1071",
           "gplots","ggplot2","grid","gridExtra","gplots","ggrepel","doRNG",
           "Hmisc","huge","irlba","igraph","irr","lme4","lmerTest","nnet",
           "Matrix","magrittr","MASS","missForest","nlme","phangorn","plot3D",
           "pheatmap","pkgconfig","plyr","parallel","pscl","plotly","rfUtilities",
           "rlang","randomForest","readxl","RColorBrewer","ROCR","reshape","reshape2","yarrr",
           "stringi","S4Vectors","tidyverse","vegan","VGAM") #"venneuler","ShortRead",
# version 1
#for (pack in packs){install.packages(sprintf("%s",pack))}
# version 2 (better version)
for (package in packages){
  if(!package %in% installed.packages()){
    install.packages(package)
  }else{library(package, character.only = TRUE)}
}

#for (package in packages){library(sprintf("%s",package), character.only = TRUE)}



##############################################
#---- install package  by bioconductor ------#
##############################################
# 210322
# install package and reads library is combined

# version 1
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioconductors <- c("phyloseq","microbiome","ANCOMBC","Rhtslib","genefilter","dada2","DESeq2", "dplyr","ggpubr","ggfortify", "ggpmisc",
                   "illuminaio","msa","rstatix","useful","DECIPHER","ComplexHeatmap")

# 


for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor)
  }else{library(bioconductor, character.only = TRUE)}
}







##############################################
#----            Introduction          ------#
##############################################

cat(blue("#--------------------------------------------------------------# \n"))
cat(blue("#------       General analysis Of microbiome (Go)        ------# \n"))
cat(blue("#------    Quick statistics and visualization tools      ------# \n"))
cat(blue("#--------------------------------------------------------------# \n"))
cat(red("                                      Version: Go_tools.3.5.0 \n"))
cat("                                              Write by Heekuk \n")
cat(yellow("All the required packages were installed.\n"))
cat(yellow("All the required packages were loaded.\n"))
cat(blue("#--------------------------------------------------------------# \n"))

#' A Go_huamnn2ps
#'
#'
#' @param huamnn2ps
#' @keywords huamnn2ps
#' @export
#' @examples
#' Go_huamnn2ps


Go_path <- function(project, pdf, table, path){
  # main dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  

  # main pdf
  if(is.null(pdf)){
    print("No pdf dir.")
  } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
    out_pdf <- file.path(sprintf("%s/pdf",out)) 
    if(!file_test("-d", out_pdf)) dir.create(out_pdf)
    
    print("pdf is in your working dir. Use /dir$pdf/ for save.")
  }

  
  
  # main table
  if(is.null(table)){
    print("No table dir.")
  } else if (table == "yes" | table == "Yes"|table == "YES"){
    out_tab <- file.path(sprintf("%s/table",out)) 
    if(!file_test("-d", out_tab)) dir.create(out_tab)

    print("table is in your working dir.Use /dir$tab/ for save.")
  }
  
  if(is.null(path)){
    print("No another dir.")
  } else if(!is.null(path)){
    out_path <- file.path(sprintf("%s/%s",out,path)) 
    if(!file_test("-d", out_path)) dir.create(out_path)
    print("path is in your working dir. Use /dir$path/ for save.")
  }


  # 한개 이상 return 하기
  
  
  functionReturningTwoValues <- function() {
    dirs <- list()
    if(is.null(pdf)){
    } else if (pdf == "yes" | pdf == "Yes"|pdf == "YES"){
      dirs$pdf <- out_pdf
    }
    if(is.null(table)){
      next
    } else if (table == "yes" | table == "Yes"|table == "YES"){
      dirs$tab <- out_tab
    }
    if(is.null(path)){
    } else if(!is.null(path)){
      dirs$path <- out_path
    }

    return(dirs) 
  }

  functionReturningTwoValues ()
  
}
Go_emptyMap <- function(psIN, project){
  
  # out dir
  map <- file.path("3_map") 
  if(!file_test("-d", map)) dir.create(map)
  
  
  # empty map table
  SampleID <- sample_names(psIN)
  Contamination <- sample_names(psIN)
  StudyID <- sample_names(psIN)
  
  emptyMap <- data.frame(SampleID, Contamination, StudyID)
  
  cat(sprintf("empty map is saved in %s.\n",map))
  cat("                                                       \n")
  write.csv(emptyMap, quote = FALSE, col.names = NA, row.names = F,
            file=sprintf("%s/emptyMap.%s.%s.csv",map, project,format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  # empty metadata table
  column.names <- c("StudyID", "Variation1", "Variation2","etc")
  col.count <- length(column.names)
  
  # 	"Go_overview","Go_ancombc","Go_deseq2","Go_box","Go_bdiv",	"Go_barchart","Go_linear","Go_clme","Go_perm",
  analysis <- c("type",	"baseline",	"Go_reg", "Go_mirkat", "Go_lmem","Confounder")

  row.count <- length(analysis)
  
  emptyMetadata <- data.frame(matrix(ncol = col.count, nrow = row.count))
  colnames(emptyMetadata) <- column.names
  rownames(emptyMetadata) <- analysis


  for(an in analysis){
    if (an == "type"){
      emptyMetadata[c(an), ] <- c("", "factor", "numeric", "factor")
    }else if(an == "baseline"){
      emptyMetadata[c(an), ] <- c("", "control", "before", "male")
    }else{
      emptyMetadata[c(an), ] <- c("no", "no", "yes", "yes")
    }
  }
  
  #cat(sprintf("empty metadata is saved in %s.\n",map))
  #cat("                                                       \n")
  #write.csv(emptyMetadata, quote = FALSE, col.names = NA,  row.names = T,
  #          file=sprintf("%s/emptyControlpanel.%s.%s.csv",map, project,format(Sys.Date(), "%y%m%d"),sep="/"))
} 








#' A Go_filter
#'
#' This function allows you to get filtering.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords filter
#' @export
#' @examples
#' Go_filter


Go_filter <- function(psIN, cutoff){ #project


  # remove 0 ASVs
  
  tt = try(psIN.prune <- prune_samples(sample_sums(psIN) > 1, psIN),T)
  if (class(tt) == "try-error"){
    psIN.prune = prune_samples(sample_sums(psIN) > 0, psIN);psIN.prune
  }else{
    psIN.prune <- prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
  }
  
  
  x <- sample_sums(psIN.prune) > 1
  cat("#--  Removing 0 sum column   --#\n")
  cat(sprintf("#--  %s column(s) are removed   --#\n", length(x[x== F])))
  cat("\n")

  phylo_relabun <- transform_sample_counts(psIN.prune, function(x) x / sum(x))
  phylo_filter <- filter_taxa(phylo_relabun, function(x) mean(x) < cutoff, TRUE) #.00005
  rmtaxa <- taxa_names(phylo_filter)
  alltaxa <- taxa_names(phylo_relabun)
  myTaxa <- alltaxa[!alltaxa %in% rmtaxa]
  phylo_relabun_filtered <- prune_taxa(myTaxa,phylo_relabun)
  ps_filtered <- prune_taxa(myTaxa,psIN.prune)

  cat("#--  Before filter  --#\n")
  print(psIN)
  cat("\n")
  cat("#--  After filter   --#\n")
  print(ps_filtered)

  prune_taxa(myTaxa,psIN)
  
  # out dir
  # out <- file.path("2_rds") 
  # if(!file_test("-d", out)) dir.create(out)
  #saveRDS(ps_filtered, sprintf("%s/ps_filtered.%s.(%s).%s.rds", out, project, cutoff,format(Sys.Date(), "%y%m%d")))
  return(ps_filtered)
  #cat("\n")
  #print(sprintf("ps_filtered is saved as 2_rds/ps_filtered.%s.(%s).%s.rds",  project, cutoff,format(Sys.Date(), "%y%m%d")))

}

#' Make a rarefaction curve using ggplot2
#' @param physeq_object A phyloseq class object, from which abundance data are extracted
#' @param step Step Size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string. The name of the variable to map to text labels on the plot. Similar to color option but for plotting text.
#' @param color Default `NULL`. Character string. The name of the variable to map to the colors in the plot. This can be a sample variables among the set returned by sample_variables(physeq_object) or taxonomic rank, among the set returned by rank_names(physeq_object)
#' @param plot default `TRUE`. Logical. Should the graph be plotted
#' @param parallel default `FALSE`. Logical. Should rarefaction be parallelized
#' @param se default `TRUE`. Logical. Should standard errors be calculated.
#' @examples
#' good_taxon_table <- data.frame(sum.taxonomy = c("a;b;c;d;f;u", "p;q;r;s;t;u"),
#' site_1 = c(0,1), site_2 = c(10, 20))
#' good_maps <- data.frame(site = c("site_1", "site_2"),
#' season = c("wet", "dry"), host = c("oak", "sage"))
#' physeq_object <- convert_anacapa_to_phyloseq(good_taxon_table, good_maps)
#' ggrare(physeq_object, step = 20, se = TRUE)
#' @export


Go_rare <- function(physeq_object, step = 10, label = NULL, color = NULL, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  p <- p + xlim(NA, xlimit) + theme_classic()
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

calRare <- function(psdata, measures, depths, parallel=FALSE) {
  require('plyr') # ldply
  require('reshape2') # melt
  require('doParallel')
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}



Go_rare2 <- function(psIN, project, alpha_metrics, color, group, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  ps.rare <- calRare(psIN, alpha_metrics, rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
  summary(ps.rare)
  ps.summary <- ddply(ps.rare, c('Depth', 'Sample', 'Measure'), summarise, 
                      Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
  
  ps.summary$Sample <-  gsub("X", "", ps.summary$Sample); ps.summary$Sample
  ps.summ.verbose <- data.frame(merge(ps.summary, data.frame(sample_data(psIN)),by.x = 'Sample', by.y = 'row.names'))
  
  sample_variables(psIN)
  ps.summ.verbose$Measure
  
  
  # Plot
  plotlist <- list()
  for (am in alpha_metrics){
    ps.summ.verbose.sel <- subset(ps.summ.verbose, Measure == am)
    p <- ggplot(ps.summ.verbose.sel, mapping = aes_string(x = "Depth", y = "Alpha_diversity_mean", 
                                                   ymin = "Alpha_diversity_mean - Alpha_diversity_sd", 
                                                   ymax = "Alpha_diversity_mean + Alpha_diversity_sd", 
                                                   colour = color, group = group))+
      xlim(NA, xlimit) + theme_light() + geom_line() + #geom_pointrange(size = 0.1) +
      theme(legend.position="right", legend.text=element_text(size=8))+ guides(col = guide_legend(ncol = 2)) +
      ggtitle(sprintf("%s-Rarefaction curve", am )) + labs(x = "Sequence Sample Size", y = am)
    # + facet_wrap(facets = ~ StudyID, scales = 'free_y')
    #print(p)
    #plotlist[[length(plotlist)+1]] <- p 
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
  #multiplot(plotlist=plotlist, cols=cols, rows=rows)

}



#' A Go_qq
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords qq plot
#' @export
#' @examples
#' Go_qq()


Go_qq <- function(psIN, project, alpha_metrics, name, height, width){
    if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
      # logic for out file
  pdf(sprintf("%s/QQplot.%s%s.%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # 1st adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  adiv$ShannonLn <-log(adiv$Shannon)
  # show last column name
  rev(names(adiv))[1]
  
  #----------- QQ plot and histogram -----------#
  par(mfrow = c(3,2))
  mes <- c(alpha_metrics, rev(names(adiv))[1])
  for (am in mes){
    test <- shapiro.test(adiv[,am])
    hist(adiv[,am], freq=F, xlab= am, main=sprintf("Histogram of %s (%s)", project, am ), cex.main=1) 
    lines(density(adiv[,am])) 
    rug(adiv[,am])
    # remove inf
    adiv.inf <- adiv[!is.infinite(adiv[,am]),]
    
    qqnorm(adiv.inf[,am], main=sprintf("Normal Q-Q Plot (%s p=%.2g)", "shapiro", test$p.value), cex.main=1)
    qqline(adiv.inf[,am])
    
    print(sprintf("%s %s shapiro test (p=%.2g)",am, project, test$p.value))
  }
  
  dev.off()
}




#' A Go_barchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_barchart()


Go_barchart <- function(psIN, cate.vars, project, taxanames, orders=NULL,
                        simple = FALSE,  
                        mycols=NULL, 
                        relative = T,
                        x_label=NULL, 
                        facet=NULL, 
                        legend="bottom", 
                        cutoff=0.005, 
                        name=NULL, 
                        ncol=NULL, 
                        height, width){
    
  if(!is.null(dev.list())) dev.off()
  
  
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

if(!is.null(x_label)){
  x_label = x_label
}else{
  x_label="SampleIDfactor"
}



  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  

if(relative == T){
  pdf(sprintf("%s/barchart.relative.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}else{
  pdf(sprintf("%s/barchart.absolute.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}




  
  
  # order by bdiv

  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))

  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  for(i in 1:length(taxanames)){

    # try table type
    otu.filt <- as.data.frame(otu_table(psIN)) 
    tt <- try(otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i]),T)
    
    if(class(tt) == "try-error"){
      print("DADA2 table")
      otu.filt <- as.data.frame(t(otu_table(psIN))) 
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }else{
      otu.filt <- as.data.frame(otu_table(psIN)) 
      print("other table")
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }
    
    

    
    if (dim(otu.filt)[2] == 2){
      next
    }
    
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames[i])), otu.filt, sum, na.action=na.pass)
    
    if (taxanames[i] == "Species"){
      agg <- agg[grepl("NA NA", agg$Species)==F,]
     
    }
    
    
    if (relative == TRUE){
      genera <- agg[,taxanames[i]]
      agg <- agg[,-1]
      agg <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      #saving table
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_relative_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                                           project,cutoff,taxanames[i],
                                                                           ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                           format(Sys.Date(),"%y%m%d"))) #,sep="/"


      df <- melt(agg, variable="SampleID")
    }else if(relative == FALSE){
      genera <- agg[,taxanames[i]]
      agg <- agg[,-1]
      agg.rel <- normalizeByCols(agg)
      inds_to_grey <- which(rowMeans(agg.rel) < cutoff)
      genera[inds_to_grey] <- "[1_#Other]"
      agg[,taxanames[i]] <- genera
      #saving table
      agg_other_out <- subset(agg, agg[,taxanames[i]] != "[1_#Other]")
      write.csv(agg_other_out, quote = FALSE, col.names = NA, file=sprintf("%s/%s.taxa_absolute_abundance.(%s).%s.%s%s.csv", out_taxa,
                                                                           project,cutoff,taxanames[i],
                                                                           ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                                           format(Sys.Date(),"%y%m%d"))) #,sep="/"
      df <- melt(agg, variable="SampleID")
    }
    



    # add StduyID
    
    
    df2 <- aggregate(as.formula(sprintf("value ~ %s + SampleID" , taxanames[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)

    #mapping.sel[df2$SampleID, "StudyID"]
   
    # add groups
    for (mvar in cate.vars) {
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
        rownames(mapping.sel) <- as.character(rownames(mapping.sel))
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
    
    if(!is.null(mycols)){
      getPalette = colorRampPalette(mycols)
    }else{
      p=p
    }
    
    
    
    
    
    
    # pdf size height = 5, width=9
   
    if (legend == "bottom"){
      if (colourCount < 30) {
        coln <- 4
      }else if (colourCount > 30) {
        coln <- 5
      }
    } else if (legend == "right") {
      if (colourCount <= 18) {
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

    p <- ggplot(df2, aes_string(x= x_label, y="value", fill=taxanames[i], order=taxanames[i])) + 
      geom_bar(stat="identity", position="stack") + theme_classic()  + labs(fill=NULL)+
      theme(legend.position=legend, # legend.text=element_text(size=8), 
            legend.text = element_text(face = c(rep("italic", 5), rep("plain", 5))),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + 
      guides(fill=guide_legend(ncol= coln))   #guides(col = guide_legend(ncol = coln)) + 
      
    
    
    if (relative == TRUE){
      p <- p + labs(y = "Relative abundance") + ylim(c(-.1, 1.01))
    }else if(relative == FALSE){
      p <- p + labs(y = "Absolute abundance")
    }
    
    

    if(!is.null(mycols)){
      p=p + scale_fill_manual(values = getPalette(colourCount)) 
    }else{
      p=p
    }
    
    
    
    if (!is.null(facet)) {
      for (mvar in cate.vars) {        
        if (facet == mvar) {
          next
        }

        df2[,facet] <- factor(df2[,facet], levels = orders)
        
        print(sprintf("Facet by %s-%s",mvar, facet))

         if (!is.null(ncol)) {
         p <- p+ facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", ncol = ncol) 
         }else{
         p <- p+ facet_grid(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales = "free_x", space = "free") 
         }





        if (!is.null(name)) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }

        print(p)
      }

    }     else if (is.null(facet) & simple == FALSE) {
      for (mvar in cate.vars) {
        print("B")
        print(sprintf("Facet by %s",mvar))

         if (!is.null(ncol)) {
         p <- p + facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales = "free_x", ncol = ncol)  
         }else{
         p <- p + facet_grid(as.formula(sprintf("~ %s"  ,mvar)), scales = "free_x", space = "free_x") 
         }

        
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("%s barplots overall of %s-%s (cut off < %s)",taxanames[i],mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("%s barplots overall of %s (cut off < %s)",taxanames[i], mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    } else if (is.null(facet) & simple == TRUE) {
      for (mvar in cate.vars) {

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

#' A Go_colbarchart
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_colbarchart()
#' 20200525
#' color for phylum

Go_colbarchart <- function(psIN, cate.vars, project, taxanames, data_type, orders, 
                        simple = FALSE,  
                        mycols=NULL, 
                        relative = T,
                        x_label=NULL, 
                        facet=NULL, 
                        legend="bottom", 
                        cutoff=0.005, 
                        name=NULL, 
                        ncol=NULL,
                        height, width){
    if(!is.null(dev.list())) dev.off()
    
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

if(!is.null(x_label)){
  x_label = x_label
}else{
  x_label="SampleIDfactor"
}

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    mycols <- NULL
  }
  
  # logic for out file
if(relative == T){
  pdf(sprintf("%s/colbarchart.relative.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}else{
  pdf(sprintf("%s/colbarchart.absolute.%s.%s%s(%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              cutoff,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
}

  taxRanks <- taxanames
  
  # order by bdiv
  ordi <- ordinate(psIN , method = "PCoA", distance = "bray")
  ordering.pc1 <- names(sort(ordi$vectors[,"Axis.1"]))
  mapping.sel <- data.frame(sample_data(psIN))

  plotlist <- list()
  for(i in 1:length(taxanames)){


    # try table type
    otu.filt <- as.data.frame(otu_table(psIN)) 
    tt <- try(otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i]),T)
    
    if(class(tt) == "try-error"){
      print("DADA2 table")
      otu.filt <- as.data.frame(t(otu_table(psIN))) 
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }else{
      otu.filt <- as.data.frame(otu_table(psIN)) 
      print("other table")
      otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=colnames(tax_table(psIN)),level=taxanames[i])
    }
    

    # continue
    otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxanames[i])
    otu.filt$PhylumCol <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks, level="Phylum")

    if (dim(otu.filt)[2] == 2){
      next
    }

    agg <- aggregate(as.formula(sprintf(". ~ %s + PhylumCol" , taxanames[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxanames[i]]
    PhylumCol <- agg$PhylumCol
    agg[,taxanames[i]] <- NULL
    agg$PhylumCol <- NULL

    agg <- normalizeByCols(agg)
    inds_to_grey <- which(rowMeans(agg) < cutoff)
    genera[inds_to_grey] <- "[1_#Other]"
    agg[,taxanames[i]] <- genera
    agg$PhylumCol <- PhylumCol 
    
    
    
    if (taxanames[i] == "Phylum"){
      agg$Phylum <- genera
    }
    
    df <- melt(agg, variable.name="SampleID")


    # add StduyID

    df2 <- aggregate(as.formula(sprintf("value ~ %s + PhylumCol + SampleID" , taxanames[i])), df, sum)
    df2$SampleID <- as.character(df2$SampleID)
    df2$SampleIDfactor <- factor(df2$SampleID, levels=ordering.pc1)
    df.SampleIDstr <- unique(df2[,c("SampleID", "SampleIDfactor")]);head(df.SampleIDstr)

    #mapping.sel[df2$SampleID, "StudyID"]

    # add groups
    for (mvar in cate.vars) {
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
    if (length(facet) == 1) {
      for (fa in facet){
        df.SampleIDstr$Group <- as.character(mapping.sel[df.SampleIDstr$SampleID, fa])
        df2[,fa] <- mapping.sel[df2$SampleID, fa]
      }
    }
    
    


    if (x_label == "SampleID"| x_label == "SampleIDfactor"){
      df2 <- df2
    } else if (length(x_label) >= 1){
      df2[,x_label] <- mapping.sel[df2$SampleID, x_label]
      df2[,x_label] <- factor(df2[,x_label], levels = orders)
    }  


    print(1)
    #------------------------#
    # ---  Color table   --- #
    #------------------------#
    agg$PhylumCol <- PhylumCol 
    agg[,taxanames[i]] <- genera
    

    #-------- remove other from taxa table --------#
    TaxaTab <- agg[order(agg[,taxanames[i]] ,  decreasing = TRUE), ]
    cdf <- data.frame(subset(TaxaTab, select=c("PhylumCol", taxanames[i])))
    cdf.sel <- subset(cdf, cdf[,taxanames[i]] != "[1_#Other]");dim(cdf.sel)[1]
    
    # 몇개인지 결정후 Phylum 으로 정리
    N <- dim(cdf.sel)[1]
    cdf.sel <- cdf.sel[order(cdf.sel$PhylumCol ,  decreasing = FALSE), ]
    cdf.sel <- data.frame(as.character(cdf.sel$PhylumCol[1:N]), as.character(cdf.sel[,taxanames[i]][1:N]))
    colnames(cdf.sel) <- c("PhylumCol", taxanames[i])
    #cdf.sel[ ,c("Kingdom","Class", "Order", "Family","Genus")] <- list(NULL)
    
    cdf.sel[,taxanames[i]] <-  gsub("p__", "", gsub("c__", "", gsub("o__", "", gsub("f__", "", gsub("g__", "", gsub("s__", "", cdf.sel[,taxanames[i]]))))))
    
    # save species name
    taxaName <- cdf.sel[,taxanames[i]]
    cdf.sel[,taxanames[i]] <- NULL
    
    # -----  create color table   ---- #
    coltab <- Go_color(cdf=cdf.sel, taxaName=taxaName)
    
    # hsv code
    #print(coltab$color_table)
    #coltab$color_table$Phylum
    # color code
    #print(coltab$coloring)
    
    print(2)
    # pdf size height = 5, width=9
    if (legend == "bottom"){
      if (N < 30) {
        col <- 5
      }
    } else if (legend == "right") {
      if (N < 18) {
        col <- 1
      }
      else if (N > 19 & N  < 35) {
        col <- 2
      }
      else if (N > 36) {
        col <- 3
      }
    }

    # plot
    # df2 <- df2[order(df2$value, decreasing=T),]
    print(3)
    level <- unique(df2[,taxanames[i]])
    #facet <- "SampleType"
    #mvar <- "TreatmentGroup"
    df2[,facet] <- factor(df2[,facet], levels = orders)
    if (length(facet) == 1) {
      for (mvar in cate.vars) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))*length(unique(df2[,facet]))
        }
        if (facet == mvar) {
          next
        }
        
        df2[,facet] <- factor(df2[,facet], levels = orders)
        print(4)
        
        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxanames[i]], levels=level), order=taxanames[i])) + geom_bar(stat="identity", position="stack") + theme_classic()  + theme(legend.position=legend, legend.text=element_text(size=8), axis.title.x = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col))  + guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + facet_wrap(as.formula(sprintf("~ %s + %s", paste(setdiff(facet, "SampleType"), collapse="+"), mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance") + labs(fill = taxanames[i])

        
        if (length(name) == 1) {
          p = p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }

        print(p)
        #plotlist[[length(plotlist)+1]] <- p
      }

    } else if (length(facet) != "NULL") {
      for (mvar in cate.vars) {
        if (class(ncol) == "numeric") {
          ncol <- ncol
        }else if(length(unique(df2[,mvar])) >= 1){
          ncol <- length(unique(df2[,mvar]))
        }

        p <- ggplot(df2, aes_string(x= x_label, y="value", fill=factor(df2[,taxanames[i]], levels=level), order=taxanames[i])) + 
        geom_bar(stat="identity", position="stack") + theme_classic()  + 
        theme(legend.position= legend, legend.text=element_text(size=8), axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + guides(fill=guide_legend(ncol= col)) + 
        guides(col = guide_legend(ncol = col)) + ylim(c(-.1, 1.01)) + scale_fill_manual(values=coltab$coloring) + 
        facet_wrap(as.formula(sprintf("~ %s"  ,mvar)), scales="free_x", ncol = ncol) + labs(y = "Relative abundance")+ labs(fill = taxanames[i])
        if (length(name) == 1) {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s-%s (cut off < %s)",mvar,name, cutoff))
        }
        else {
          p= p+ ggtitle(sprintf("Taxa barplots overall of %s (cut off < %s)",mvar, cutoff))
        }
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  #multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
#' A Go_adiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Alpha diversity plot
#' @export
#' @examples
#' Go_adiv


Go_adiv <- function(psIN, project, alpha_metrics){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_adiv <- file.path(sprintf("%s_%s/table/adiv",project, format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", out_adiv)) dir.create(out_adiv)
  
  
  # adiv table
  mapping.sel <- data.frame(sample_data(psIN))
  adiv <- estimate_richness(psIN, measures=alpha_metrics) # se.chao1 stand error
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping.sel)
  adiv <- merge(adiv, mapping.sel, by="row.names")
  rownames(adiv) <- adiv$SampleID
  cat(sprintf("adiv table is saved in %s.\n",out_path))
  cat("                                                       \n")
  write.csv(adiv, quote = FALSE, col.names = NA, 
            file=sprintf("%s/adiv.%s.%s.csv",out_adiv,project, format(Sys.Date(), "%y%m%d"),project,format(Sys.Date(), "%y%m%d"),sep="/"))
  return(adiv)
} 
  
#' A Go_box_plot
#'

Go_boxplot <- function(df, cate.vars, project, outcomes,
                       orders=NULL, 
                       mycols=NULL, 
                       combination=NULL,
                       ylim =NULL,
                       title= NULL, 
                       facet= NULL, 
                       paired=NULL, 
                       name= NULL, 
                       addnumber=TRUE,
                       statistics = "yes", 
                       parametric= "no", 
                       star="no",
                       xanlgle=90,  
                       height, width, plotCols, plotRows){
  
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  
  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  
  
  pdf(sprintf("%s/box.%s.%s%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # plot
  plotlist <- list()
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    if (length(facet) >= 1){
      if (facet == mvar){
        next
      }
    } else {}
    
    # remove Na
    print("Control NA")
    df <- data.frame(df)
    df[,mvar] <- as.character(df[,mvar]);df[,mvar]
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    df.na[,mvar] <- as.factor(df.na[,mvar]);df.na[,mvar]  
    
    df.na[,mvar] <- factor(df.na[,mvar], levels = intersect(orders, df.na[,mvar]))
    
    
    # Add number of samples in the group
    print("Add sample number informations")
    if(addnumber==TRUE){
    renamed_levels <- as.character(levels(df.na[,mvar]));renamed_levels
    oldNames <- unique(df.na[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(df.na[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(df.na[[mvar]])[levels(df.na[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    }else{
      df.na <- df.na
    }

    
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(df.na)[1], dim(df)[1]))
    
    if (length(unique(df.na[,mvar])) ==1) {
      next
    }
    
    summary.df.na <- summary(df.na[,mvar])
    
    #------------------------------#
    # for group combination or not #
    #------------------------------#
    
    
    if (!is.null(combination)){
      print(sprintf("Combination n=", combination))
      group.cbn <- combn(x = levels(df.na[,mvar]), m = combination)
      
      #print(count(group.cbn))
      
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
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar)) 
        } else if(combination ==3){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2)) 
        }else if(combination ==4){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3)) 
        }else if(combination ==5){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4)) 
        }else if(combination ==6){
          basline <- group.combination[1]
          smvar1 <- group.combination[2]
          smvar2 <- group.combination[3]
          smvar3 <- group.combination[4]
          smvar4 <- group.combination[5]
          smvar5 <- group.combination[6]
          df.cbn <- subset(df.na, df.na[,mvar] %in% c(basline,smvar1, smvar2,smvar3,smvar4,smvar5)) 
        }else{
          print("combination should be 2, 3, 4, 5, and 6 only.")
          break
        }
        
        unique(df.cbn[,mvar])
        
        
        # make a comnination for stat
        df.cbn[,mvar] <- factor(df.cbn[,mvar])
        cbn <- combn(x = levels(df.cbn[,mvar]), m = 2)
        
        
        my_comparisons <- {}
        for(i in 1:ncol(cbn)){
          x <- cbn[,i]
          my_comparisons[[i]] <- x
        };my_comparisons
        
        if(combination != 2){
          combination.N <- combination - 1
          my_comparisons <- my_comparisons[1:combination.N]
        }
        
        
        
       
        for(oc in outcomes){
          # remove NA for facet
          if (length(facet) >= 1) {
            for (fc in facet){
              df.cbn[,fc] <- as.character(df.cbn[,fc]);df.cbn[,fc]
              df.cbn[,fc][df.cbn[,fc] == ""] <- "NA"
              df.cbn.sel <- df.cbn[!is.na(df.cbn[,fc]), ]
              df.cbn <- df.cbn.sel 
              # facet or not
              df.cbn[,fc] <- factor(df.cbn[,fc], levels = orders)
            }
          }

           # check statistics method
          if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
            if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
                test.name <- "ANOVA"
                testmethod <-  "t.test"
              } else {
                testmethod <-  "t.test"
                pval <- NULL
                test.name <- "Pairwise T-Test"
              }
            }else{
              if (nlevels(factor(df.cbn[,mvar])) > 2) {
                test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.cbn)
                pval <- round(test$p.value, 4)
                test.name <- "KW"
                testmethod <- "wilcox.test"
              } else {
                testmethod <- "wilcox.test" 
                pval <- NULL
                test.name <- "Pairwise Wilcoxon"
              }
            } 
          }else{
            test.name<-NULL
            pval <- NULL
          }

  
          p1 <- ggplot(df.cbn, aes_string(x=mvar, y=oc, colour=mvar))  + labs(y=oc, x=NULL) + 
            theme_bw() + theme(strip.background = element_blank()) +
            theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
                  plot.title=element_text(size=9)) # ,face="bold"

          
          if(!is.null(mycols)){
            p1 <- p1 + scale_color_manual(values = mycols)
          }else{
            p1 <- p1
          }
          
        
          if (!is.null(title)) {
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")), 
                                       ifelse(is.null(pval), "", paste("p=", " ", sep = "")), 
                                       ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
          } else{
            p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                       ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")), 
                                       ifelse(is.null(pval), "", paste("p=", " ", sep = "")), 
                                       ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
          }
          
          # control statistic on the plot
          
          if(is.null(test.name)){
            p1 <- p1 
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < 0.05){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {  
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
                }  else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1 
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {  
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
              }  else if (star == "yes") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1 
            }
          }

          
          

          
          # plot design
          if (height*width <= 6){
            dot.size = 0.7
            box.tickness = 0.3
          }else if (height*width > 6 & height*width < 10){
            dot.size = 1
            box.tickness = 0.4
          }else{
            dot.size = 1.5
            box.tickness = 0.5
          }
          

          if(!is.null(ylim)){
            if(oc == "Chao1"){
              p1 = p1
            }else{
              p1 = p1 + ylim(ylim[1] , ylim[2])   
            }
          }


          # paired plot type
          if (!is.null(paired)) {
            #p1 = p1 + geom_point(size = 1) 
            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
            p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3),show.legend = F)  #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
            p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3)) 
            p1 = p1  + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm")) 
            
          }  else{
            p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
      
            # count or table for number of variable
            if (max(table(df.cbn[,mvar])) > 250 & max(table(df.cbn[,mvar])) < 500){
              dot.size <- dot.size/2
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            } else  if (max(table(df.cbn[,mvar])) < 250 ){
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
            }else if(max(table(df.cbn[,mvar])) > 500) {
              dot.size <- dot.size/3
              p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))
            }
          } 
          
          # facet
          if (length(facet) >= 1) {
            facetCol <- length(unique(df[,facet]))
            p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          } else {
            p1 = p1 + guides(color = "none", size = "none", shape= "none")
          }
          plotlist[[length(plotlist)+1]] <- p1 
        }
      }
    }else{
      # make a comnination for stat
      
      print("Check combination for statistics")
      cbn <- combn(x = levels(df.na[,mvar]), m = 2)
      
      my_comparisons <- {}
      for(i in 1:ncol(cbn)){
        x <- cbn[,i]
        my_comparisons[[i]] <- x
      };my_comparisons
      # check statistics method
      for(oc in outcomes){
        # remove NA for facet
        if (!is.null(facet)) {
          for (fc in facet){
            df.na[,fc] <- as.character(df.na[,fc]);df.na[,fc]
            df.na[,fc][df.na[,fc] == ""] <- "NA"
            df.na.sel <- df.na[!is.na(df.na[,fc]), ]
            df.na <- df.na.sel 
            # facet or not
            df.na[,fc] <- factor(df.na[,fc], levels = orders)
          }
        }
        
        
        if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
          if (parametric == "yes"| parametric == "YES"|parametric == "Yes"){
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- aov(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(summary(test)[[1]][["Pr(>F)"]][1],4)
              test.name <- "ANOVA"
              testmethod <-  "t.test"
            } else {
              testmethod <-  "t.test"
              pval <- NULL
              test.name <- "Pairwise T-Test"
            }
          }else{
            if (nlevels(factor(df.na[,mvar])) > 2) {
              test <- kruskal.test(as.formula(sprintf("%s ~ %s", oc, mvar)), df.na)
              pval <- round(test$p.value, 4)
              test.name <- "KW"
              testmethod <- "wilcox.test"
            } else {
              testmethod <- "wilcox.test" 
              pval <- NULL
              test.name <- "Pairwise Wilcoxon"
            }
          } 
        }else{
          test.name<-NULL
          pval <- NULL
        }

        
        p1 <- ggplot(df.na, aes_string(x=mvar, y=oc, colour=mvar))  + labs(y=oc, x=NULL) + 
          theme_bw() + theme(strip.background = element_blank()) +
          theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
                plot.title=element_text(size=9)) #,face="bold"  


        # scale_color_brewer(palette=colorset)
        
        if(!is.null(mycols)){
          p1 <- p1 + scale_color_manual(values = mycols)
        }else{
          p1 <- p1
        }
        
        
        # Close an image
        if (!is.null(title)) {
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", title,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")), 
                                     ifelse(is.null(pval), "", paste("p=", " ", sep = "")), 
                                     ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
        } else{
          p1 <- p1 + ggtitle(sprintf("%s%s%s%s", mvar,
                                     ifelse(is.null(test.name), "", paste("\n",test.name, " ", sep = "")), 
                                     ifelse(is.null(pval), "", paste("p=", " ", sep = "")), 
                                     ifelse(is.null(pval), "", paste(pval, " ", sep = "")), sep=""))
        }
        
        # control statistic on the plot
        
        if(is.null(paired)){
          if(is.null(test.name)){
            p1 <- p1 
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < 0.07){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {  
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
                } else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1 
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {  
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2)
              }  else if (star == "yes") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3)
              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1 
            }
          }
        }else{
          print("paired")
          
          if(is.null(test.name)){
            p1 <- p1 
          } else if(test.name == "KW" | test.name == "ANOVA"){
            if(pval < 0.07){
              if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
                if (star == "no") {  
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2,paired = TRUE)
                } else if (star == "yes") {
                  p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3,paired = TRUE)
                }
              }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
                p1 <- p1 
              }
            }else {
              p1 <- p1
            }
          }else if(testmethod == "wilcox.test" | testmethod == "t.test"){
            if (statistics == "yes"| statistics == "YES"|statistics == "Yes"){
              if (star == "no") {  
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.format", comparisons = my_comparisons, size = 2,paired = TRUE)
              }  else if (star == "yes") {
                p1 <- p1 + stat_compare_means(method= testmethod, label = "p.signif", comparisons = my_comparisons, hide.ns = TRUE, size = 3,paired = TRUE)
              }
            }else if(statistics == "no"| statistics == "NO"|statistics == "No"){
              p1 <- p1 
            }
          }
        }
        
      
        # plot design
        if (height*width <= 6){
          dot.size = 0.7
          box.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 1
          box.tickness = 0.4
        }else{
          dot.size = 1.5
          box.tickness = 0.5
        }
        
        # y axis limit      
        if(!is.null(ylim)){
          if(oc == "Chao1"){
            p1 = p1
          }else{
            p1 = p1 + ylim(ylim[1] , ylim[2])   
          }
        }
        
        # paired plot type
        if (!is.null(paired)) {
          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          p1 = p1 + geom_point(aes_string(fill=mvar,group=paired),alpha = 0.8, size = dot.size, position = position_dodge(0.3), show.legend = F)   #scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
          p1 = p1 + geom_line(aes_string(group=paired), color="grey50", size=0.3,position = position_dodge(0.3)) 
          p1 = p1 + theme(legend.title = element_blank(), legend.position="bottom", legend.justification="left",legend.box.margin = ggplot2::margin(0,0,0,-1,"cm")) 
          
        } else{
          p1 = p1 + geom_boxplot(aes_string(colour=mvar),outlier.shape = NA,lwd=box.tickness)  + theme(legend.position="none")
          
          # count or table for number of variable
          if (max(table(df[,mvar])) > 250 & max(table(df[,mvar])) < 500){
            dot.size <- dot.size/2
            p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
          } else  if (max(table(df[,mvar])) < 250 ){
            p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2)) # alpha=0.3
          }else if(max(table(df[,mvar])) > 500) {
            dot.size <- dot.size/3
            p1 = p1 + geom_jitter(aes_string(colour=mvar),shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2))
          }
        } 
        # facet
        if (length(facet) >= 1) {
          facetCol <- length(unique(df[,facet]))
          p1 = p1 + facet_wrap(as.formula(sprintf("~ %s" , paste(setdiff(facet, "SocpleType"), collapse="+"))), scales="free_x", ncol = facetCol) 
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        } else {
          p1 = p1 + guides(color = "none", size = "none", shape= "none")
        }
        
        plotlist[[length(plotlist)+1]] <- p1 
      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
#' A Go_clme
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_clme()


Go_clme <- function(psIN, cate.vars, project, paired, mycols=NULL, addnumber=TRUE,node, decreasing, height,timepoint,ID, orders,xangle, name, width, plotCols, plotRows){
    
  if(!is.null(dev.list())) dev.off()
    
  alpha_metrics = c("Chao1","Shannon")
  
  # Descriptions 분석 하고자 하는 variation에 subgroup
  # paired 환자나 같은 사람 ID
  # node 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # decreasing 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }

  pdf(sprintf("%s/clme.%s.%s(%s.%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              node,
              decreasing,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  # adiv
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  mapping <-data.frame(sample_data(psIN))
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping)
  adiv <- merge(adiv, mapping, by="row.names"); rownames(adiv) <- adiv$SampleID

  if (length(orders) >= 1) {
    adiv[,timepoint] <- factor(adiv[,timepoint], levels = orders)
  }

  
  # clme
  cons <- list(order = "umbrella" , node=node, decreasing = decreasing) 
  # 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  print(cons)
  
  plotlist <- list()
  for (mvar in cate.vars) {
    print(mvar)
    
    if (length(unique(adiv[,mvar])) < 2){
      next
    }
    
    
    
    # Na 제거

    adiv[,mvar] <- data.frame(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]
    
    
    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    adiv <- adiv.na
    


    # Add number of samples in the group
    if(addnumber==TRUE){
    renamed_levels <- as.character(levels(adiv[,mvar]));renamed_levels
    oldNames <- unique(adiv[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(adiv[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(adiv[[mvar]])[levels(adiv[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    }else{
      adiv <- adiv
    }




    if (mvar == timepoint){
      for (am in alpha_metrics){
        form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))
        
        clme.mod <- clme(form, data = adiv, constraints = cons, seed = 2, nsim = 1000)
        clme.sum <- summary(clme.mod, seed=2)
        clme.globalp <- function(model) { label <- substitute(
          italic(p) == globalp,
          list(globalp <- model$p.value) )
        as.character(as.expression(format(globalp, nsmall=3))) 
        }
        
        clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))
        
        
        # plot design
        if (height*width <= 6){
          dot.size = 0.6
          line.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 0.9
          line.tickness = 0.4
        }else{
          dot.size = 1.3
          line.tickness = 0.5
        }
        
        # plot
        p <- ggplot(adiv, mapping = aes_string(x=timepoint, y=am, color=timepoint, group=paired)) + 
          geom_line(color="grey",size=line.tickness,position = position_dodge(0.3)) + 
          geom_point(alpha = 0.8, size = dot.size,position = position_dodge(0.3)) + ylab(sprintf("%s Index\n", am)) + 
          ggtitle(sprintf("%s \n (%s) ", mvar, clme.globalp))  + 
          theme_bw() + theme(strip.background = element_blank()) + 
          theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" )           
          
          
          if(!is.null(mycols)){
           p <- p + scale_color_manual(values = mycols)
           }else{
           p <- p
           }


        
        if (length(ID) == 1) {
          p= p + geom_text_repel(aes_string(label = ID), size = 2)
        }
        plotlist[[length(plotlist)+1]] <- p
      }
    }else{
      for (des in unique(adiv[,mvar])){
        if(dim(subset(adiv, adiv[,mvar] == des))[1] < 3){
          next
        }
        if(timepoint == mvar){
          next
        }
        print(des)
        for (am in alpha_metrics){
          form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))
          
          clme.mod <- clme(form, data = adiv[adiv[,mvar] == des,], constraints = cons, seed = 2, nsim = 1000)
          clme.sum <- summary(clme.mod, seed=2)
          clme.globalp <- function(model) { label <- substitute(
            italic(p) == globalp,
            list(globalp <- model$p.value) )
          as.character(as.expression(format(globalp, nsmall=3))) 
          }
          
          clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))


                  # plot design
        if (height*width <= 6){
          dot.size = 0.6
          line.tickness = 0.3
        }else if (height*width > 6 & height*width < 10){
          dot.size = 0.9
          line.tickness = 0.4
        }else{
          dot.size = 1.3
          line.tickness = 0.5
        }
          
          # plot
          p <- ggplot(adiv[adiv[,mvar]==des,], mapping = aes_string(x=timepoint, y=am, color=timepoint, group=paired)) + 
            geom_line(color="grey",size=line.tickness,position = position_dodge(0.3)) + geom_point(alpha = 0.8, size = dot.size,position = position_dodge(0.3)) + 
            xlab(timepoint) + ylab(sprintf("%s Index\n", am)) + 
            ggtitle(sprintf("%s-%s \n (%s) ", mvar, des, clme.globalp))   + theme_bw() +
            theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" ) 


           if(!is.null(mycols)){
           p <- p + scale_color_manual(values = mycols)
           }else{
           p <- p
           }


          
          if (length(ID) == 1) {
            p= p + geom_text_repel(aes_string(label = ID), size = 2)
          }
          plotlist[[length(plotlist)+1]] <- p
        }
      }
    }
    

  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}


#' A Go_linear
#'
#' This function help to generation of linear regression plots
#' @param test test
#' @keywords basic statistic and simple plots
#' @export
#' @examples

Go_linear <- function(df, cont.vars, project, outcomes, method="lm",
                      mycols =NULL, maingroup=NULL, orders=NULL, name=NULL, 
                      height, width, plotCols, plotRows){
    
  if(!is.null(dev.list())) dev.off()
    
print("Method option: glm, lm, loess, gam")
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  pdf(sprintf("%s/linear.%s.%s.%s%s%s.pdf", out_path, 
              project, 
              method,
              ifelse(is.null(maingroup), "", paste(maingroup, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  

  my.formula <- y ~ x

  
  # plot
  plotlist <- list()
  for (mvar in cont.vars) {
    if (!is.null(maingroup)){
      if (mvar == maingroup){
        next
      }
    }
    # Na 제거
    df[,mvar] <- as.numeric(as.character(df[[mvar]]))
    #df[,mvar] <- as.numeric(as.character(df$percenttwl))
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df[,mvar]<- as.numeric(df[,mvar]);df[,mvar]
    
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(df.na)[1], dim(df)[1]))
    if (length(unique(df.na[,mvar])) ==1) {
      next
    }
    summary.df.na <- summary(df.na[,mvar])
    
    # na제거 in the maingroup
    if (!is.null(maingroup)) {
      df.na[,maingroup] <- as.character(df.na[,maingroup]);df.na[,maingroup]
      df.na[,maingroup][df.na[,maingroup]==""] <- "NA";df.na[,maingroup]
      df.na[,maingroup]<- as.factor(df.na[,maingroup]);df.na[,maingroup]
      df.na.na <- subset(df.na, df.na[,maingroup] != "NA");df.na.na[,maingroup]
      
      if (!is.null(orders)) {
        df.na.na[,maingroup] <- factor(df.na.na[,maingroup], levels = orders)
      }else{
        df.na.na[,maingroup] <- factor(df.na.na[,maingroup])
      }

    }

    
    
    
    for(i in 1:length(outcomes)){
      
      if (outcomes[i] == mvar | outcomes[i] == "Chao1" & mvar == "Shannon" | outcomes[i] == "Shannon" & mvar == "Chao1") {
        print(sprintf("Stop function bacause out was %s and mvar was %s", outcomes[i], mvar))
        next
      }
      
      df.na[,outcomes[i]] <- as.numeric(df.na[,outcomes[i]])
      
      print(outcomes[i])
      
      if (!is.null(maingroup)) {
        df.na.na[,outcomes[i]] <- as.numeric(df.na.na[,outcomes[i]])
        p<- ggplot(df.na.na, aes_string(x=mvar, y=outcomes[i], group= maingroup, color=maingroup, linetype = maingroup))
          
      }else {
        p <- ggplot(df.na, aes_string(x=mvar, y=outcomes[i]))+theme_classic() + geom_point(size = 0.5)
        
      }
      
     p <- p + theme_classic() + geom_point(size = 0.5)
       # scale_colour_brewer(palette = colorset) + 

       if(method == "glm"){
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5, method.args = list(family = "poisson"))   #method.args = list(family = "poisson")
       }else if(method == "lm"){
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5) + 
           #stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE, size = 3) +
           stat_fit_glance(method.args = list(formula = my.formula), method = method, 
                           #geom = 'text', 공식이 한쪽으로 정리가 되지 않고, 라인에 수치가 붙는다.
                           aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g', 
                                               stat(r.squared), stat(p.value))), 
                           parse = TRUE, size = 3)
       }else{
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5)
       }
     
     p <- p + ggtitle(sprintf("%s with %s", mvar, outcomes[i])) + theme(title=element_text(size=10)) + labs(x = NULL)+
       theme(title=element_text(size=10),
             axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
             axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))
     
      
      if(!is.null(mycols)){
        p <- p + scale_color_manual(values = mycols)
      }else{
        p <- p
      }



      plotlist[[length(plotlist)+1]] <- p
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}

Go_regression <- function(data, project,  
                          outcomes, 
                          uni.vars=NULL, 
                          mul.vars=NULL,
                          interaction=NULL, 
                          orders, pvalue=0.05, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/regression",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  
  
  # data control
  
  # fix column types
  if(!is.null(uni.vars) | !is.null(mul.vars)){
    
    varis <- unique(c(uni.vars,mul.vars))
    
    for (mvar in  varis) {
      if (class(data$mvar) == "character") {
        data[,mvar] <- factor(data[,mvar])
        # setting baseline
        data[,mvar] <- factor(data[,mvar], levels = intersect(orders, data[,mvar]))
        # NA 제거
        data[,mvar] <- as.character(data[[mvar]]);data[,mvar]
        data[,mvar][data[,mvar]==""] <- "NA";data[,mvar]
        #data.na <- subset(data, data[,mvar] != "NA");data.na[,mvar]  # subset 를 사용한 NA 삭제
        
      } else if (class(data$mvar) == "numeric" | min(mvar) < 0) {
        # NA 제거
        data[,mvar] <- as.character(data[[mvar]]);data[,mvar]
        data[,mvar][data[,mvar]==""] <- "NA";data[,mvar]
        #data.na <- subset(data, data[,mvar] != "NA");data.na[,mvar]  # subset 를 사용한 NA 삭제
        data[,mvar] <- as.numeric(as.character(data[[mvar]]))
      } 
    }
  }



  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes){
    if (length(unique(data[,outcome])) > 3){
      
      
      
    }else{
      if (class(outcome) == "character") {
        data[,outcome] <- factor(data[,outcome])
        data[,outcome] <- factor(data[,outcome], levels = intersect(orders, data[,outcome]))
        
        # NA 제거
        data[,outcome] <- as.character(data[[outcome]]);data[,outcome]
        data[,outcome][data[,outcome]==""] <- "NA";data[,outcome]
        #data.na <- subset(data, data[,outcome] != "NA");data.na[,outcome]  # subset 를 사용한 NA 삭제
        # set the baseline for outcome
        
        if(length(unique(data[,outcome])) == 2){
          
          data[,outcome] <- factor(data[,outcome])
          out <- levels(data[,outcome])[1]
          
          
          data[,outcome] <- factor(ifelse(data[,outcome]== levels(data[,outcome])[1],0,1), levels=c(0,1), 
                                   labels = levels(data[,outcome]))
          
          
          print(levels(data[,outcome]))
        }
      } else if (class(outcome)  == "numeric") {
        
        # NA 제거
        data[,outcome] <- as.character(data[[outcome]]);data[,outcome]
        data[,outcome][data[,outcome]==""] <- "NA";data[,outcome]
        # data <- subset(data, data[,mvar] != "NA");data[,outcome]  # subset 를 사용한 NA 삭제
        data[,outcome] <- as.numeric(as.character(data[[outcome]]))
        
        data[,outcome] <- as.numeric(data[,outcome])
      } 
      
      res <- {}
      
      for (mvar in uni.vars) {
        
        if (outcome == mvar) {
          next
        }
        
        # print(sprintf("##-- %s (total without NA: %s/%s) --##", 
        #              mvar, dim(data.na)[1], dim(data)[1]))
        
        if (length(unique(data[,mvar])) ==1) {
          next
        }
        
        # get formula
        if(!is.null(uni.vars) & is.null(mul.vars) & is.null(interaction)){
          form <- as.formula(sprintf("%s ~ %s", outcome, mvar))
          print("Univariate anaysis")
          type <- "uni"
        } else if (!is.null(mul.vars)){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(mul.vars, interaction)
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+")))
            print("Multivariate anaysis")
            type <- "multi-interantion"
            
          }else{
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars, "SampleType"), collapse="+")))
            
            print("Multivariate anaysis with interaction effect")
            type <- "multi"
          }
        } 
        
        print(form)
        
        
        # model
        if (class(data[,outcome]) == "numeric"){
          m <- "Regression (glm-poisson)"
          mod <- glm(form, data=data,  family = poisson(link='log'))
        } else if (length(unique(data[,outcome])) == 2){
          m <- "Logistic regression (glm-binomial)"
          mod <- glm(form, data=data,  family=binomial("logit"))
        } 
        print(m)
        
        
        # out for the model
        coef <- as.data.frame(summary(mod)$coefficients)
        coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
        colnames(coef) <- c("Estimate", "SE", "t", "pval")
        
        if (dim(coef)[1] == 0){
          next
        }
        
        
        
        # out for the confidence interval 
        conf <- data.frame(confint(mod))
        conf <- conf[setdiff(rownames(conf), "(Intercept)"),,drop=F]
        conf.na <- na.omit(conf) 
        
        coef$`2.5 %` <- conf.na$`2.5 %`
        coef$`97.5 %` <- conf.na$`97.5 %`
        
        coef$outcome <- outcome
        coef$mvar <- mvar
        coef$model <- m
        coef$deviance <- pchisq(q=mod$null.deviance-mod$deviance,df=mod$df.null-mod$df.residual, lower.tail = FALSE)
        
        # get formula
        if (!is.null(mul.vars)){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(mul.vars, interaction)
            mul.inter.form <- sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+"))
            coef$formula <- mul.inter.form
          }else{
            mul.form <- sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars, "SampleType"), collapse="+"))
            coef$formula <- mul.form
          }
        } else{
          coef$mvar <- mvar
        }
        
        
        res <- rbind(res, coef)
        
        
        # stop looing for multivariate analysis
        if(!is.null(mul.vars) | !is.null(interaction)){
          break
        }
      }
      
      
      res$padj <- p.adjust(res$pval, method="fdr")
      #res <- res[order(res$time_point),]
      res$comp <- factor(rownames(res), levels=rownames(res))
      res$dir <- ifelse(res$pval < pvalue, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
      
      print(res)
      
      write.csv(res, quote = FALSE, col.names = NA,file=sprintf("%s/regression_%s.%s.%s.%s%s.csv",out_table,
                                                                project,
                                                                outcome,
                                                                type,
                                                                ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                                                                format(Sys.Date(), "%y%m%d"), sep="/"))
      # return model
      if(!is.null(mul.vars) | !is.null(interaction)){
        saveRDS(mod,sprintf("%s/regression_%s.%s.%s.%s%s.rds",out_table,
                            project,
                            outcome,
                            type,
                            ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                            format(Sys.Date(), "%y%m%d"), sep="/")) 
        
      }
    }
  }
}






#' A Go_box_plot
#'

Go_dualYplot <- function(df, TaxTab, cate.vars, project,  Box, Line1, Line2=NULL,
                       title= NULL, 
                       name= NULL, 
                       mycols=NULL, 
                       orders=NULL,
                       xanlgle=90,  height, width){
  
  if(!is.null(dev.list())) dev.off()

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  


  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  # out file
  pdf(sprintf("%s/dualYplot.%s.%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(Line1), "", paste(Line1, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)


  ## fix factor  and  numeric
  df$etc <- NULL
  df2 <- read.csv(sprintf("%s",TaxTab),header=T,row.names=1,check.names=FALSE);head(df2)
  
  rownames(df2)<- df2$Species
  df2$Species <-NULL
  rownames(df2) <- gsub(" ","_",rownames(df2));rownames(df2)
  
  
  df2 <- as.data.frame(t(df2))

  for (var in cate.vars) {
    print(var)
      df[,var] <- factor(df[,var])
  }
  

  
  # plot
  for (mvar in cate.vars) {
    if (length(unique(df[,mvar])) < 2){
      next
    }
    
    
    # merge merged.df and taxa table
    merged.df <- merge(df, df2, by="row.names");head(merged.df)
    
    # NA remove
    merged.df[,mvar] <- as.character(merged.df[,mvar]);merged.df[,mvar]
    merged.df[,mvar][merged.df[,mvar]==""] <- "NA";merged.df[,mvar]
    merged.df.na <- subset(merged.df, merged.df[,mvar] != "NA");merged.df.na[,mvar]  # subset 를 사용한 NA 삭제
    merged.df.na[,mvar] <- as.factor(merged.df.na[,mvar]);merged.df.na[,mvar]  
    
    # re-order
    if (length(orders) >= 1) {
      merged.df.na[,mvar] <- factor(merged.df.na[,mvar], levels = orders)
    } else {
      merged.df.na[,mvar] <- factor(merged.df.na[,mvar])
    }
    
    # Add number of samples in the group
    renamed_levels <- as.character(levels(merged.df.na[,mvar]));renamed_levels
    oldNames <- unique(merged.df.na[,mvar]);oldNames
    if (length(renamed_levels) == 0) {
      renamed_levels <- oldNames
    }
    for (name in oldNames) {
      total <- length(which(merged.df.na[,mvar] == name));total
      new_n <- paste(name, " (n=", total, ")", sep="");new_n
      levels(merged.df.na[[mvar]])[levels(merged.df.na[[mvar]])== name] <- new_n
      renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
    }
    
    
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(merged.df.na)[1], dim(merged.df)[1]))
    
    if (length(unique(merged.df.na[,mvar])) ==1) {
      next
    }
    
    summary.merged.df.na <- summary(merged.df.na[,mvar])
    

    
    #===============================#
    # Visualization for Dual Y axis #
    #===============================#
  
    # for Line1
    mean.line1 <- aggregate(merged.df.na[,Line1], list(merged.df.na[,mvar]), FUN=mean)
    colnames(mean.line1) <- c(mvar, Line1);mean.line1
    mean.line1[,Line1] <- mean.line1[,Line1]*10

    if (height*width <= 6){
      dot.size = 0.7
      box.tickness = 0.3
    }else if (height*width > 6 & height*width < 10){
      dot.size = 1
      box.tickness = 0.4
    }else{
      dot.size = 1.5
      box.tickness = 0.5
    }
    
    p <- ggplot() + theme_bw() + theme(strip.background = element_blank()) + #theme_ipsum() +
      geom_boxplot(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), outlier.shape = NA, show.legend = FALSE) +
      theme(text=element_text(size=9), axis.text.x=element_text(angle=xanlgle,hjust=1,vjust=0.5),
            plot.title=element_text(size=9,face="bold"))
      # theme(legend.position="none") +

     if(!is.null(mycols)){
      p <- p + scale_color_manual(values = mycols)
     }else{
       p <- p
     }
    
    
    # count or table for number of variable
    if (max(table(merged.df.na[,mvar])) > 250 & max(table(merged.df.na[,mvar])) < 500){
      dot.size <- dot.size/2
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE) 
    } else  if (max(table(merged.df.na[,mvar])) < 250 ){
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE)  
    }else if(max(table(merged.df.na[,mvar])) > 500) {
      dot.size <- dot.size/3
      p = p + geom_jitter(data=merged.df.na, mapping=aes(x=!!sym(mvar), y=!!sym(Box), colour=!!sym(mvar)), 
                          shape=16, alpha = 0.8, size = dot.size, position=position_jitter(0.2), show.legend = FALSE) 
    }
    
    mean.line1.melt <- melt(mean.line1)
    # get label for geom_text_repel
    label.line1 <- subset(mean.line1.melt, variable == Line1);
    n <- dim(label.line1)[1]
    label.line1.sel <- label.line1[n,]
    label.line1.sel$variable <- gsub("_"," ", label.line1.sel$variable)
    
    
    
    p1 <- p + geom_line(data = mean.line1.melt, 
                        mapping = aes(x = !!sym(mvar), y = value, group=variable, color= variable), 
                        inherit.aes = FALSE, size=1)  + guides(color = "none") +
      scale_linetype_manual(values=c("solid", "solid")) + theme(legend.position="top")+
      geom_text_repel(data = label.line1.sel, aes(x = !!sym(mvar), y = value,label = variable),
                      size=3, fontface="italic")
    

   
      
    # for Line2
    if (!is.null(Line2)){
      mean.line1 <- aggregate(merged.df.na[,Line1], list(merged.df.na[,mvar]), FUN=mean)
      colnames(mean.line1) <- c(mvar, Line1);mean.line1
      mean.line1[,Line1] <- mean.line1[,Line1]*10
      
      
      mean.line2 <- aggregate(merged.df.na[,Line2], list(merged.df.na[,mvar]), FUN=mean)
      colnames(mean.line2) <- c(mvar, Line2);mean.line1
      mean.line2[,Line2] <- mean.line2[,Line2]*10
      
      mean.line <- merge(mean.line1, mean.line2, by=mvar);head(mean.line)
      
      mean.line.melt <- melt(mean.line)

      # get label for geom_text_repel
      label.line1 <- subset(mean.line.melt, variable == Line1);
      n <- dim(label.line1)[1]
      label.line1.sel <- label.line1[n,]
      label.line1.sel$variable <- gsub("_"," ", label.line1.sel$variable)
      
      label.line2 <- subset(mean.line.melt, variable == Line2);
      n <- dim(label.line2)[1]
      label.line2.sel <- label.line2[n,]
      label.line2.sel$variable <- gsub("_"," ", label.line2.sel$variable)
      
      
      p1 <- p + geom_line(data = mean.line.melt, 
                          mapping = aes(x = !!sym(mvar), y = value, group=variable, color= variable), 
                          inherit.aes = FALSE, size=1)  + guides(color = "none") +
        scale_linetype_manual(values=c("solid", "solid")) + theme(legend.position="top") +
        geom_text_repel(data = label.line1.sel, aes(x = !!sym(mvar), y = value,label = variable),
                        size=3, fontface="italic")  +
        geom_text_repel(data = label.line2.sel, aes(x = !!sym(mvar), y = value,label = variable),
                        size=3, fontface="italic") 
          
      
      

    }
    

    
    p1 <- p1 + scale_y_continuous(sec.axis = sec_axis(~.*10, name="Relative abundance (%)")) 
    
    if (!is.null(title)) {
      p1 <- p1 + ggtitle(title)
    } else{
      p1 <- p1 + ggtitle(sprintf("%s", mvar))
    }
    print(p1)
  }
  dev.off()
}

#' A Go_bdiv
#'


Go_bdiv <- function(psIN, cate.vars, project, orders, distance_metrics,
                    cate.conf=NULL,
                    plot="PCoA",
                    ellipse="yes",
                    mycols=NULL,
                    combination=NULL,
                    shapes = NULL, 
                    ID = NULL,  
                    facet=NULL, 
                    name=NULL, 
                    addnumber=TRUE,
                    height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  
  pdf(sprintf("%s/ordi.%s.%s%s%s%s%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(combination), "", paste("(cbn=",combination, ").", sep = "")), 
              ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)

  
  plotlist <- list()
  for (mvar in cate.vars) {
    mapping <- data.frame(sample_data(psIN))
    mapping[,mvar] <- factor(mapping[,mvar])

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
    
    #------------------------------#
    # for group combination or not #
    #------------------------------#
    
    if (!is.null(combination)){
      mapping[,mvar] <- factor(mapping[,mvar], levels = intersect(orders, mapping[,mvar]))
      group.cbn <- combn(x = levels(mapping[,mvar]), m = combination)
      
      #print(count(group.cbn))
      
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
          mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
          
          
          ord_meths = plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
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


           # Add number of samples in the group
          if(addnumber==TRUE){
              renamed_levels <- as.character(levels(pdataframe[,mvar]));renamed_levels
              oldNames <- unique(pdataframe[,mvar]);oldNames
          if (length(renamed_levels) == 0) {
            renamed_levels <- oldNames
          }
          for (name in oldNames) {
            total <- length(which(pdataframe[,mvar] == name));total
            new_n <- paste(name, " (n=", total, ")", sep="");new_n
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== name] <- new_n
            renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
          }
          }else{
            pdataframe <- pdataframe
          }
    








          # Plots
          p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
          
          
          if (!is.null(shapes)) {
            
            pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
            p = p +  geom_point(aes_string(shape=shapes), size=0.8, alpha = 1) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
            
          }else{
            p = p + geom_point(size=0.8, alpha = 1)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
          }
          
          p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
          p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
          p = p + theme(legend.position = "bottom", 
                        legend.title = element_blank(),
                        legend.justification="left", 
                        legend.box = "vertical",
                        legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"))
          
          
          if(!is.null(mycols)){
            p <- p + scale_color_manual(values = mycols)
          }else{
            p <- p
          }
          
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
          
          
          #===================================#
          # Add permanova for two combination #
          #===================================#
          set.seed(1)
          distance <- Go_dist(psIN = psIN.cbn.na, project = project, distance_metrics = distance_metric)
          
          x <- as.dist(distance[[distance_metric]])
          factors <-  mapping.sel.na.rem[,mvar]
          
          R2 <- c()
          p.value <- c()
          F.Model <- c()
          pairs <- c()
          SumsOfSqs <- c()
          Df <- c()
          
          x1=as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]
          
          # run
          map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))
          
          # count to table
          
          if (!is.null(cate.conf)) {
            for(conf in cate.conf){
              map.pair[,conf] <- factor(map.pair[,conf])
            }
            form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
            print(form)
          }else{
            form <- as.formula(sprintf("x1 ~ %s", mvar))
            print(form)
          }
          
          ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL
          
          Df <- c(Df,ad[1,1])
          SumsOfSqs <- c(SumsOfSqs, ad[1,2])
          R2 <- round(c(R2,ad[1,3]), digits=3)
          F.Model <- c(F.Model,ad[1,4]);
          p.value <- c(p.value,ad[1,5])
          
          pairw.res <- data.frame(Df,SumsOfSqs,R2,F.Model,p.value)
          
          class(pairw.res) <- c("pwadonis", "data.frame")
          # end adonis end
          tmp <- as.data.frame(pairw.res)
          tmp$distance_metric <- distance_metric
          tmp$padj <- p.adjust(tmp$p.value, method="bonferroni")
          
          grob <- grobTree(textGrob(paste(distance_metric, "\nR2=",R2,"\nPERMANOVA p=",tmp$padj,sep=""), x=0.01,  y=0.15, hjust=0,
                                    gp=gpar(fontsize=8))) #, fontface="italic"

          
          if (table(map.pair[,mvar])[1] <=2 | table(map.pair[,mvar])[2] <=2){
            p=p
          }else{
            p = p + annotation_custom(grob)
          }
          
          
          
          #plotlist[[length(plotlist)+1]] <- p
          
          print(p)
          
        }
      }
    }  else{
      for(distance_metric in distance_metrics){
        # remove na
        mapping.sel <- data.frame(sample_data(psIN))
        #mapping.sel[mapping.sel==""] <- "NA"
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
        mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
        
        
        
        ord_meths= plot # c("DCA", "CCA", "RDA", "DPCoA", "NMDS","PCoA")
        plist = llply(as.list(ord_meths), function(i, psIN.na, distance_metric){
          ordi = ordinate(psIN.na, method=i, distance=distance_metric)
          plot_ordination(psIN.na, ordi, type = "samples", color= mvar)
        }, psIN.na, distance_metric)
        
        names(plist) <- ord_meths
        
        pdataframe = ldply(plist, function(x){
          df = x$data[, 1:2]
          colnames(df) = c("Axis_1", "Axis_2")
          return(cbind(df, x$data))
        })
        names(pdataframe)[1] = "method"
        
        pdataframe[,facet] <- factor(pdataframe[,facet], levels = orders)
        
        pdataframe[,mvar] <- factor(pdataframe[,mvar], levels = orders)
        
           # Add number of samples in the group
          if(addnumber==TRUE){
              renamed_levels <- as.character(levels(pdataframe[,mvar]));renamed_levels
              oldNames <- unique(pdataframe[,mvar]);oldNames
          if (length(renamed_levels) == 0) {
            renamed_levels <- oldNames
          }
          for (name in oldNames) {
            total <- length(which(pdataframe[,mvar] == name));total
            new_n <- paste(name, " (n=", total, ")", sep="");new_n
            levels(pdataframe[[mvar]])[levels(pdataframe[[mvar]])== name] <- new_n
            renamed_levels <- replace(renamed_levels, renamed_levels == name, new_n);renamed_levels
          }
          }else{
            pdataframe <- pdataframe
          }
        
        
        # Plots
        p = ggplot(pdataframe, aes_string("Axis_1", "Axis_2", color=mvar))
        
        
        if (!is.null(shapes)) {
          pdataframe[,shapes] <- factor(pdataframe[,shapes], levels = orders)
          p = p +  geom_point(aes_string(shape=shapes), size=0.8, alpha = 1) + scale_shape_manual(values = c(1, 16, 8, 0,15, 2,17,11, 10,12,3,4,5,6,7,8,9,13,14)) 
          
        }else{
          p = p + geom_point(size=0.8, alpha = 1)+ ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
        }
        
        p = p + ggtitle(sprintf("%s (%s)",mvar,distance_metric)) 
        p = p + facet_wrap(~ method, scales="free") + theme_bw() + theme(strip.background = element_blank())# open(1), cross(10), closed(2)
        p = p + theme(legend.position = "bottom", 
                      legend.title = element_blank(),
                      legend.justification="left", 
                      legend.box = "vertical",
                      legend.box.margin = ggplot2::margin(0,0,0,-1,"cm"),
                      plot.title=element_text(size=8,face="bold"))
        
        if(!is.null(mycols)){
          p <- p + scale_color_manual(values = mycols)
        }else{
          p <- p
        }
        
        # ID variation
        if (!is.null(ID)) {
          p <- p + geom_text_repel(aes_string(label = ID), size = 2)
        } else {
          p <- p 
        }
        
        # ellipse variation
        if (ellipse == "yes" | ellipse == "Yes" ) {
          p <- p + stat_ellipse(type = "norm", linetype = 2) 
        } else if (ellipse == "no" | ellipse == "No" ){
          p <- p 
        }
        
        if (!is.null(facet)) {
          ncol <- length(unique(mapping.sel.na.rem[,facet]))
          p <- p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
        }
        else {
          p <- p
        }
        
        #===================================#
        # Add permanova                     #
        #===================================#
        set.seed(1)
        distance <- Go_dist(psIN = psIN.na, project = project, distance_metrics = distance_metric)
        
        x <- as.dist(distance[[distance_metric]])
        factors <-  mapping.sel.na.rem[,mvar]
        
        R2 <- c()
        p.value <- c()
        F.Model <- c()
        pairs <- c()
        SumsOfSqs <- c()
        Df <- c()
        
        x1=as.matrix(x)[factors %in% unique(factors), factors %in% unique(factors)]
        
        # run
        map.pair <- subset(mapping.sel.na.rem, mapping.sel.na.rem[,mvar] %in% unique(factors))
        
        # count to table
        
        if (!is.null(cate.conf)) {
          for(conf in cate.conf){
            map.pair[,conf] <- factor(map.pair[,conf])
          }
          form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
          print(form)
        }else{
          form <- as.formula(sprintf("x1 ~ %s", mvar))
          print(form)
        }
        
        ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL
        
        Df <- c(Df,ad[1,1])
        SumsOfSqs <- c(SumsOfSqs, ad[1,2])
        R2 <- round(c(R2,ad[1,3]), digits=3)
        F.Model <- c(F.Model,ad[1,4]);
        p.value <- c(p.value,ad[1,5])
        
        pairw.res <- data.frame(Df,SumsOfSqs,R2,F.Model,p.value)
        
        class(pairw.res) <- c("pwadonis", "data.frame")
        # end adonis end
        tmp <- as.data.frame(pairw.res)
        tmp$distance_metric <- distance_metric
        tmp$padj <- p.adjust(tmp$p.value, method="bonferroni")
        
        grob <- grobTree(textGrob(paste(distance_metric, "\nR2=",R2,"\nPERMANOVA p=",tmp$padj,sep=""), x=0.01,  y=0.15, hjust=0,
                                  gp=gpar(fontsize=8))) #, fontface="italic"
        
        if (table(map.pair[,mvar])[1] <=2 | table(map.pair[,mvar])[2] <=2){
          p=p
        }else{
          p = p + annotation_custom(grob)
        }
        
        
        #plotlist[[length(plotlist)+1]] <- p
        print(p)
      }
    }
  }
  dev.off()
}
#' A Go_bdiv
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Go_bdiv()

Go_dist <- function(psIN, project, distance_metrics){
  # out dir
  #out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out)) dir.create(out)
  #out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_path)) dir.create(out_path)
  #out_perm <- file.path(sprintf("%s_%s/table/perm ",project, format(Sys.Date(), "%y%m%d"))) 
  #if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  # run distance
  dm <- list()
  for (distance_metric in distance_metrics) {
    dm[[length(dm)+1]] <- phyloseq::distance(psIN, method=distance_metric)
  }
  
  names(dm) <- distance_metrics
  class(dm)
  
  return(dm)
}
  #' A Go_perm
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Mar 07 2020
#' adjsted 기능을 추가 하였다.
#' 분석 할때 마다 수치가 조금 변하는 것을 수정 하였다.set.seed(1)
#' dm를 따로 분리 하여 시간을 단축 하였고, dm를 다른 방법으로 분석 할수 있게 되었다.
#' Go_perm()


Go_perm <- function(psIN, cate.vars, project, distance, distance_metrics, cate.conf=NULL, des, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s/table",out)) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_perm <- file.path(sprintf("%s/perm",out_path)) 
  if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  out_distance <- file.path(sprintf("%s/distance",out_path)) 
  if(!file_test("-d", out_distance)) dir.create(out_distance)
  

  # Run
  if (!is.null(des)) {
    # Uni
    print(sprintf("#--- Running Paired-PERMANOVA (%s) ---#", des))
  }  else {
    print("#--- Running Paired-PERMANOVA  ---#")
  }
  set.seed(1)
  mapping.sel <-data.frame(sample_data(psIN))
  
  res.pair <-{}
  # Run
  for (mvar in cate.vars) {
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    if (length(unique(mapping.sel.na[,mvar])) == 1){
      cat(sprintf("there is no group campare to %s\n",unique(mapping.sel[,mvar])))
      next
    }
    for (distance_metric in distance_metrics) {
      
      psIN.sel <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel.na[,mvar]), ]), psIN)
      
      
      
      ## fix factor  and  numeric
      mapping.sel.na[,mvar] <- factor(mapping.sel.na[,mvar])
      
      distance <- Go_dist(psIN = psIN.sel, project = project, distance_metrics = distance_metric)
      
      
      # pairwise.adonis2
      # pair.ado <- pairwise.adonis2(x=as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar], map=mapping.sel.na, cate.conf=adjust, mvar=mvar)
      
      x <- as.dist(distance[[distance_metric]])
      factors <-  mapping.sel.na[,mvar]
      map <- mapping.sel.na
      
      co <- combn(unique(as.character(map[,mvar])),2)
      R2 <- c()
      p.value <- c()
      F.Model <- c()
      pairs <- c()
      SumsOfSqs <- c()
      Df <- c()
      
      
      
      
      
      for(elem in 1:ncol(co)){
        x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                        factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
        
        # run
        map.pair <- subset(map, map[,mvar] %in% c(co[1,elem],co[2,elem]))
        
# count to table
        if (table(map.pair[,mvar])[co[1,elem]] <=2 | table(map.pair[,mvar])[co[2,elem]] <=2){
          next
        }
        
        if (!is.null(cate.conf)) {
          form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
          print(form)
        }else{
          form <- as.formula(sprintf("x1 ~ %s", mvar))
          print(form)
        }
        
        ad <- adonis2(form, data = map.pair, permutations=999, by="terms")# "terms"  "margin" NULL
  
        pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
        Df <- c(Df,ad[1,1])
        SumsOfSqs <- c(SumsOfSqs, ad[1,2])
        R2 <- c(R2,ad[1,3])
        F.Model <- c(F.Model,ad[1,4]);
        p.value <- c(p.value,ad[1,5])
      }
      
      pairw.res <- data.frame(pairs,Df,SumsOfSqs,R2,F.Model,p.value)
      
      class(pairw.res) <- c("pwadonis", "data.frame")
      # end adonis end
      tmp <- as.data.frame(pairw.res)
      tmp$distance_metric <- distance_metric
      tmp$mvar <- mvar
      tmp$adjusted <- paste(setdiff(cate.conf, "SampleType"), collapse="+")
      res.pair <- rbind(res.pair, tmp)
    }
  }
  
  res.pair$padj <- p.adjust(res.pair$p.value, method="bonferroni")
  
  res.pair <- res.pair[,c("pairs", "Df","SumsOfSqs","R2","F.Model", "p.value", "padj", "distance_metric","mvar", "adjusted")]
  
  
  # output
    write.csv(res.pair, quote = FALSE,col.names = NA, sprintf("%s/pair_permanova.%s.%s%s%s%s.csv", out_perm, 
              project, 
              ifelse(is.null(cate.conf), "", paste(cate.conf, "adjusted.", sep = "")), 
              ifelse(is.null(des), "", paste(des, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")),sep="/")
  

  return(res.pair)
}

Go_mirkat<- function(psIN, project, cate.vars, cate.conf = NULL,  orders,name=NULL){
  # install bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  bioconductors <- c("data.table","phyloseq","dirmult","vegan")
  
  for (bioconductor in bioconductors){
    if(!bioconductor %in% installed.packages()){
      library(BiocManager)
      BiocManager::install(bioconductor)
    }else{library(bioconductor, character.only = TRUE)}
  }
  
  

  # install package from install_github
  github <- c("GLMMMiRKAT")
  if(!github %in% installed.packages()){
    install_github("hk1785/GLMM-MiRKAT", force=T)
  }else{library(package, character.only = TRUE)}
  # install package 
  packages <- c("data.table","CompQuadForm","devtools","ecodist","GUniFrac","GLMMMiRKAT","lme4","MASS","Matrix","MiRKAT","permute") 
  for (package in packages){
    if(!package %in% installed.packages()){
      install.packages(package)
    }else{library(package, character.only = TRUE)}
  }
  
  

  #install.packages("data.table", version = "1.13.0")
  
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/mirkat",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  

  
  # check by variables
  mapping <- data.frame(sample_data(psIN))
  for (mvar in  cate.conf) {
   mapping[,mvar] <- factor(mapping[,mvar])
  }
  
  sample_data(psIN) <- mapping
  
  
  res <- {}
  for (mvar in cate.vars) {
    if (length(unique(mapping[,mvar])) == 1){
      print(sprintf("%s has only 1 variation, which wouldn't be able to compare.",mvar))
      next
    }
    #-----------------------#
    # for group combination #
    #-----------------------#
    # order for each variation
    mapping[,mvar] <- factor(mapping[,mvar], levels = orders)
    mapping[,mvar] <- factor(mapping[,mvar])
    
    # list of the combination
    group.cbn <- combn(x = levels(mapping[,mvar]), m = 2)
    
    #print(count(group.cbn))
    
    group_comparisons <- {}
    for(i in 1:ncol(group.cbn)){
      x <- group.cbn[,i]
      group_comparisons[[i]] <- x
    };group_comparisons
    
    set.seed(1)
    
    # looping for using the combination of list 
    for(i in 1:length(group_comparisons)){
      print(group_comparisons[i])
      group.combination <- unlist(group_comparisons[i]);group.combination
      
      basline <- group.combination[1]
      smvar <- group.combination[2]
      mapping.cbn <- subset(mapping, mapping[,mvar] %in% c(basline,smvar)) 
      
      psIN.cbn <- psIN
      sample_data(psIN.cbn) <- mapping.cbn
      dim(psIN.cbn)
      
      # extract phyloseq elements
      otu.cbn <- data.frame(otu_table(psIN.cbn))
      mapping.cbn <- data.frame(sample_data(psIN.cbn))
      tree.cbn <- phy_tree(psIN.cbn)
      
      
      # create empty df for this combination
      df <- data.frame(matrix(ncol = 1, nrow = dim(mapping.cbn)[1]));df[1] <- list(NULL)
      df.covar <- df
      
      df[,mvar] <- as.numeric(mapping.cbn[,mvar]  == unique(mapping.cbn[,mvar] )[1])
      
      # add corvatiate into the df
      
      if (!is.null(cate.conf)){
      for (covar in cate.conf) {
        df.covar[,covar] <- as.numeric(mapping.cbn[,covar])
        if (mvar == covar){
        next
        }
      }
      }

      
      # Create the UniFrac Distances
      unifracs <- GUniFrac(otu.cbn, tree.cbn, alpha = c(0, 0.5, 1))$unifracs
      D.weighted = unifracs[,,"d_1"]
      D.unweighted = unifracs[,,"d_UW"]
      D.generalized = unifracs[,,"d_0.5"]
      D.BC = as.matrix(vegdist(otu.cbn, method="bray"))
      
      # Convert Distance to Kernel Matrices
      K.weighted = D2K(D.weighted)
      K.unweighted = D2K(D.unweighted)
      K.generalized = D2K(D.generalized)
      K.BC = D2K(D.BC)
      Ks = list(K.weighted = K.weighted, K.unweighted = K.unweighted, K.BC = K.BC)
      
      
      if (length(cate.conf) >=1){
        for (covar in cate.conf) {
          # Cauchy
          # cauchy <- MiRKAT(y = df[,mvar], Ks = Ks, X = df.covar, out_type = "D", method = "davies",  
          # omnibus = "cauchy", returnKRV = FALSE, returnR2 = FALSE)
          # print(cauchy)  tt <- try()
          
          
          tt <- try(permutation <- MiRKAT(y = df[,mvar], Ks = Ks, X = df.covar, out_type = "D", method = "davies"
                                          ,omnibus = "permutation", returnKRV = FALSE, returnR2 = FALSE), T)
          
          if (class(tt) == "try-error"){
            print(sprintf("Fail comparison %s vs %s",basline,smvar))
            
            next
          }
        }
        
      }else if(length(cate.conf) ==0){
        
        permutation <- MiRKAT(y = df[,mvar], Ks = Ks, X = NULL, out_type = "D", method = "davies", 
                              omnibus = "permutation", returnKRV = FALSE, returnR2 = FALSE)
      }
      
      # create table
      per.df <- data.frame(unlist(permutation));per.df
      
      colnames(per.df) <- paste(basline,"vs",smvar);per.df
      per.df.t <- as.data.frame(t(per.df));per.df.t
      
      per.df.t$mvar <- mvar
      class(per.df.t)
      
      
      if (length(cate.conf) >=1){
        covars <- cate.conf[mvar != cate.conf]
        per.df.t$cate.conf <- paste(setdiff(cate.conf, "SampleType"), collapse="+")
        per.df.t$covar <- paste(setdiff(df.covar, "SampleType"), collapse="+")
        
      }else if(length(cate.conf) ==0){
        per.df.t <- per.df.t
      }
      
      
      res <- rbind(res, per.df.t);res
      print(res)
    }
  }
  write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/mirkat.%s.%s%s.csv",out_table,
                                                          project,
                                                          ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                                                          format(Sys.Date(), "%y%m%d"), sep="/")) 
  
}
Go_dist_plot <- function(psIN, project, distance_metrics, distance, group,orders, name,height,width,plot = TRUE) {
    
  if(!is.null(dev.list())) dev.off()
   
  colorset = "Paired" # Dark1 Set1
  

  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  if (plot == TRUE) {
    if (length(name) == 1) {
      pdf(sprintf("%s/5_distplot.%s.%s.%s.pdf", out_path, project, name,format(Sys.Date(), "%y%m%d")),height = height, width=width)
    } else{
      pdf(sprintf("%s/5_distplot.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")),height = height, width=width)
    }
  
  # calc distances
  for (dist in distance_metrics){
    wu <- as.dist(distance[[dist]])
    wu.m <- melt(as.matrix(wu))
    
    # remove self-comparisons
    wu.m <- wu.m %>%
      filter(as.character(Var1) != as.character(Var2)) %>%
      mutate_if(is.factor, as.character)
    
    
    wu.m  = wu.m  %>%
      rowwise() %>%      # for each row
      mutate(Samples = paste(sort(c(Var1, Var2)), collapse = "-")) %>%  # sort the teams alphabetically and then combine them separating with -
      ungroup()
    
    wu.m.sel  = distinct(wu.m, Samples, .keep_all=T)
    
    
    # get sample data (S4 error OK and expected)
    mapping <- data.frame(sample_data(psIN))
    mapping$ID <- as.character(rownames(mapping))
    mapping[,group] <- as.character(mapping[,group])

    
    sd <- mapping %>%
      select("ID", group) %>%
      mutate_if(is.factor,as.character)
    
    # combined distances with sample data
    # sample1
    colnames(sd) <- c("Var1", "Type1")
    wu.m.sel$Var1 <- factor(wu.m.sel$Var1)
    sd$Var1 <- factor(sd$Var1)
    wu.sd <- left_join(wu.m.sel, sd, by = "Var1")
    # sample2
    wu.sd$Var2 <- as.factor(wu.sd$Var2)
    colnames(sd) <- c("Var2", "Type2")
    sd$Var2 <- factor(sd$Var2)
    wu.sd <- left_join(wu.sd, sd, by = "Var2")
    
    wu.sd$Type3 <- ifelse(wu.sd$Type1 == wu.sd$Type2, wu.sd$Type1,"across_group")
    wu.sd$Type3 <- factor(wu.sd$Type3)
    # make a combination for stat
    wu.sd.sel <- subset(wu.sd, Type3 != "across_group")
    wu.sd.sel$Type3 <- factor(wu.sd.sel$Type3)

    # cbn <- combn(x = levels(wu.sd$Type3), m = 2)
    baseline <- "across_group"
    cbn <-{}
    for (x in levels(wu.sd.sel$Type3)){
      cbn <- cbind(cbn, c(baseline, x))
    }
    
    cbn.sel<-cbn[, !duplicated(t(cbn))]
    
    my_comparisons <- {}
    for(i in 1:ncol(cbn.sel)){
      x <- cbn.sel[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    
    # plot
    if (length(orders) >= 1) {
      wu.sd$Type2 <- factor(wu.sd$Type2, levels = orders)
    }       else {
      wu.sd$Type2 <- factor(wu.sd$Type2)
    }
    
    p <- ggplot(wu.sd, aes(x = Type3, y = value, colour=Type3)) + theme_classic()+
      geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, alpha = 0.3,position=position_jitter(0.2)) +
      scale_color_brewer(palette=colorset) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      #facet_wrap(~ Type2, scales = "free_x") +
      ylab(dist) + xlab(NULL) +  theme(legend.position="none")+
      stat_compare_means(method= "wilcox.test", label = "p.format", comparisons = my_comparisons, size = 2.5)
    
    if (length(name) == 1) {
      p<- p+ ggtitle(sprintf("%s_%s",dist,name))
    } else{
      p <- p+ ggtitle(sprintf("%s",dist))
    }
    print(p)
  }
    dev.off()
    
  } else {
    return(wu.sd)
  }
}
#' A Go_deseq2_heat
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 Heatmap
#' @export
#' @examples
#' Go_deseq2_heat()

Go_mergeTab <- function(pattern, file_path){
  
   # add input files
  path <- file_path
 
  
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  filenames <- list.files(path, pattern=pattern);filenames
  
  
  cat(sprintf("Files location: %s\n",path))
  cat("=======================================================================\n")
  cat("Merged files:\n")
  cat(sprintf("%s\n",filenames))

  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  return(df)
}
#' A Go_lmem
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Run LMEM
#' @export
#' @examples
#' Go_lmem()

## thresholds for association/permutation tests
# nsamps_threshold <- 0.01 fraction of relabund to call a sample positive
#filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
nperm <- 100000


Go_lmem <- function(psIN, cate.vars, cate.conf=NULL, StudyID, project, pval=0.05, nsamps_threshold, filt_threshold, taxanames, name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/lmem",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  out_lmem.Tab <- file.path(sprintf("%s_%s/table/lmem/tab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_lmem.Tab)) dir.create(out_lmem.Tab)
  
  out_lmem.ps <- file.path(sprintf("%s_%s/table/lmem/ps",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_lmem.ps)) dir.create(out_lmem.ps)
  
  
  

  #taxRanks <- taxanames

  psIN.agg <- aggregate_taxa(psIN, taxanames);psIN.agg
  
  for (cvar in cate.conf) {
    df[, cvar] <- mapping.sel[df$SampleID, cvar]
  }
  
  for (mvar in cate.vars) {
    
  # combination
  mapping.sel <- data.frame(sample_data(psIN))
  mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = orders)
  
  mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
  cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
  
  my_comparisons <- {}
  for(i in 1:ncol(cbn)){
    x <- cbn[,i]
    my_comparisons[[i]] <- x
  };my_comparisons
  
  
  for(i in 1:length(my_comparisons)){
    print(my_comparisons[i])
    combination <- unlist(my_comparisons[i]);combination
    baseline <-combination[1];baseline
    smvar <- combination[2];smvar
    
    mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(baseline, smvar));dim(mapping.sel.cb) # phyloseq subset은 작동을 안한다.
    psIN.cb <- psIN.agg
    sample_data(psIN.cb) <- mapping.sel.cb
    
    for(i in 1:length(taxanames)){
      # dada2 or nephele
      # try table type
      otu.filt <- as.data.frame(t(otu_table(psIN.cb)))
      tt <- try(otu.filt[,rank]  <- getTaxonomy(otus=rownames(otu.filt), taxRanks = colnames(tax_table(psIN.cb)), tax_tab=tax_table(psIN.cb), level=rank),T)
      
      if(class(tt) == "try-error"){
        print("other table")
        otu.filt <- as.data.frame(otu_table(psIN.cb)) 
        otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN.cb), taxRanks=colnames(tax_table(psIN.cb)),level=taxanames[i])
      }else{
        otu.filt <- as.data.frame(t(otu_table(psIN.cb)))
        print("DADA2 table")
        otu.filt[,taxanames[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN.cb), taxRanks=colnames(tax_table(psIN.cb)),level=taxanames[i])
      }
      
    
      agg <- aggregate(as.formula(sprintf(". ~ %s" , taxanames[i])), otu.filt, sum, na.action=na.pass)
      genera <- agg[,taxanames[i]]
      
      agg <- agg[,-1]
      #agg <- normalizeByCols(agg)
      rownames(agg) <- genera
      dim(agg)
      ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
      agg <- agg[intersect(ftk,ftk),]
      # control data set after filter
      if (dim(agg)[1] == 0)
        next
      
      agg[,taxanames[i]] <- rownames(agg)
      
      
      # metatable에서 useForNB는 여러개의 yes가 가능 하지만, useAsConfounder 는 그렇지 않다.
      ## baseline등을 관리 하려면 다음이 필요하다.

      
      
      print(2)
      #--------------    lmer    -------------#
      res <- {}
      for (f in agg[,taxanames[i]]) {
        # clean bacteria name
        if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
          next
        }
        
        df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)
        df$StudyID <- mapping.sel.cb[df$SampleID, StudyID]
        
        
        
        # add groups
        for (cate in cate.conf) {
          df$Group <- as.character(mapping.sel[df$SampleID, cate])
          df[,cate] <- mapping.sel[df$SampleID, cate]
          
          # order
          if (length(orders) >= 1) {
            df[,mvar] <- factor(df[,cate], levels = orders)
          }
          else {
            df[,mvar] <- factor(df[,cate])
          }
        }
        
        
          # na remove

          mapping.sel.cb[mapping.sel.cb==""] <- "NA"
          mapping.na <- mapping.sel.cb[!is.na(mapping.sel.cb[,mvar]), ]
          na.count <- length(mapping.na)
          if (length(unique(mapping.na[,mvar])) == 1)
            next
          
          
          #------------ fix column types------------#
           mapping.na[,mvar] <- factor(mapping.na[,mvar])
          #if (metadata[mvar, "type"] == "factor") {
          #  mapping.na[,mvar] <- factor(mapping.na[,mvar])
          #  if (length(unique(mapping.na[,mvar])) ==1 ){
          #    next
          #  }
          #  #if (metadata[mvar, "baseline"] != "") {
          #  #  mapping.na[,mvar] <- relevel(mapping.na[,mvar], metadata[mvar, "baseline"])
          #  #}
          #} else if (metadata[mvar, "type"] == "numeric") {
          #  mapping.na[,mvar] <- factor(mapping.na[,mvar])
          #}
         

          print(3)
          
          
          # na count
          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.na)[1], dim(mapping)[1]))
          print(4)
          
          df[,mvar] <- mapping.na[df$SampleID, mvar]
          
          #=====================#
          #  Regression method  #
          #=====================#
          if(!is.null(StudyID)){
            reg <- "LMEM"
            form <- as.formula(sprintf("value ~ (1 | StudyID) + %s  %s", mvar, 
                                       ifelse(is.null(cate.conf), "", paste("+",setdiff(cate.conf, "SampleType"), collapse=""))))
            print(form)
            tt <- try(mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore")),T)
            if(class(tt) == "try-error"){
              next
            }else{
              mod <- lmer(form, data=df, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
              # mod2 <- lmer(form, data=df)
            }
          }else{
            reg <- "GLM"
            form <- as.formula(sprintf("value ~ %s  %s", mvar, 
                                       ifelse(is.null(cate.conf), "", paste("+",setdiff(cate.conf, "SampleType"), collapse=""))));form
            print(form)
            #mod <- glm(form, data=df,  family = binomial(link='logit'))
            #mod <- glm(form, data=df,  family = quasipoisson(link = "log"))
            mod <- glm(form, data=df,  family = poisson(link = "log"))
            #summary(mod)
            # pchisq(summary(mod.over)$dispersion * mod.ori$df.residual, mod.ori$df.residual, lower.tail = FALSE) #과산포 test p < 0.05 
          }
            


          #exp(coef(mod))
          ## lmer에서 control=은 "number of levels of each grouping~~" 오류가 있을때만 사용한다.
          ##

          
          coef <- summary(mod)$coefficients
          coef <- coef[grep(mvar, rownames(coef)),,drop=F]
          res <- rbind(res, cbind(f, mvar, rownames(coef), coef, baseline))
          dim(res)
        }
      }
      
      
      #-- create table --#
      res <- as.data.frame(res)
      #colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue", "baseline")
      print(5)
      
      if(!is.null(StudyID)){
        reg <- "LMEM"
        res$pvalue <- as.numeric(as.character(res$`Pr(>|t|)`))
        res$`Pr(>|t|)` <- NULL
      }else{
        reg <- "GLM"
        res$pvalue <- as.numeric(as.character(res$`Pr(>|z|)`))
        res$`Pr(>|z|)` <- NULL
      }
      
      res$Estimate <- as.numeric(as.character(res$Estimate))
      res$SE <- as.numeric(as.character(res$`Std. Error`))
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$method <- reg
      
      res <- res[order(res$pvalue),]
      res.sel <- res
      res.sel <- as.data.frame(subset(res, pvalue < pval))
      taxa_sig <- res.sel$taxa[1:dim(res.sel)[1]]; summary(taxa_sig)
      
      if(dim(res.sel)[1] == 0){
        next
      }else{
        res.sel$bas.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == baseline)))
        res.sel$coef.count <-  unique(sum(with(mapping.na, mapping.na[,mvar] == smvar)))
      }
      

      
       if(dim(res.sel)[1] == 0){
        ps.taxa.sig <- psIN.cb
      }else{
        tt <- try(ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb),T)
        
        if(class(tt) == "try-error"){
          pathwayTab <- data.frame(otu_table(psIN.cb))
          pathwayRank <- data.frame(tax_table(psIN.cb))
          rownames(pathwayRank) <- pathwayRank[,taxRanks]
          rownames(pathwayTab) <- pathwayRank[,taxRanks]
          pathwayRank <- as.matrix(pathwayRank)
          pathwayTab <- as.matrix(t(pathwayTab))
          psIN.cb <- phyloseq(otu_table(pathwayTab, taxa_are_rows=FALSE), tax_table(pathwayRank));psIN.cb
          ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
          print(ps.taxa.sig)
        }else{
          ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
          print(ps.taxa.sig)
        }
      }
      

      #res.sel <- arrange(res.sel, res.sel$pvalue)



      
      write.csv(res.sel, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).Sig%s.%s.%s.%s%s.%s.csv",out_lmem.Tab,
                                                               baseline, 
                                                               smvar,
                                                               dim(res.sel)[1],
                                                               taxanames[i], 
                                                               mvar,
                                                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                               project, 
                                                               reg, 
                                                               sep="/"))
      saveRDS(ps.taxa.sig,sprintf("%s/(%s.vs.%s).Sig%s.%s.%s.%s%s%s.%s.rds",out_lmem.ps,
                                  baseline, 
                                  smvar,
                                  dim(res.sel)[1],
                                  taxanames[i], 
                                  mvar, 
                                  ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                  ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                  project,
                                  reg))
  }
 }
}

#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_lmem_fore()
# dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")
# dircolors <- c("#4f86f7", "#e10000", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_lmem_fore <- function(project,file_path, mycols=NULL, fdr=0.05,est=1, files, name, order, height, width){
    if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=files);filenames
  
  print(path)
  print(filenames)
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  
  # out file
  pdf(sprintf("%s/lmem.forest.%s.%s(fdr=%s.est=%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr,
              est,
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  
  
  
  
  
  for (fn in 1:length(filenames)) {
    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")

    df.sel <- df
    
    df.sel.sel <- as.data.frame(subset(df.sel, abs(df.sel$Estimate) > est))
    
    resSig <- as.data.frame(subset(df.sel.sel, pvalue < fdr )); resSig <- resSig[order(resSig$Estimate),]
    # resSig$smvar <- factor(resSig$smvar)
    print(1)
    if (dim(resSig)[1] == 0)
      next

    print(2)
    resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)== 1, "up", "down"), "NS")
    
  
    
    print(3)
    
    # 중복 이름 처리 하기
    headers <- vector(dim(resSig)[1], mode="character")
    
    for (i in 1:dim(resSig)[1]) {
      headers[i] <- paste("ASV", i, sep="_")
    }
    resSig$ASV <- headers
    
    
    for (plot in unique(resSig$metadata)){
      resSig.sel <- subset(resSig, metadata == plot)
      if (length(unique(resSig.sel$coefficient)) ==1 ){
        
        baseline <- unique(resSig.sel$baseline)
        compare <- unique(resSig.sel$coefficient)
        compare <- gsub(unique(resSig.sel$metadata), "", compare);compare
        
        
        # color
        resSig.sel$dir <- gsub('down',baseline, gsub('up',compare, resSig.sel$dir))
        resSig.sel$dir <- factor(resSig.sel$dir, levels = c(as.character(baseline), "NS", as.character(compare)))
        
        # color end
        
        if(!is.null(mycols)){
          dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
        }else{
          dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
        }
        
        legend.labs <- 
          c(paste(baseline, " (n=", unique(resSig.sel$bas.count),")",sep=""),
            paste("NS"),
            paste(compare, " (n=", unique(resSig.sel$coef.count), ")",sep=""))
        
        
        
        
        lims <- max(abs(resSig.sel$Estimate) + abs(resSig.sel$SE))*1.0
        p1 <- ggplot(resSig.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
          geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
          geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
          scale_color_manual(values=dircolors,labels=legend.labs)  + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel$ASV), labels = resSig.sel$taxa)+theme(plot.title = element_text(size=8, hjust = .5))
        
        #sn <- gsub("\\..*","",sn);sn
        p1 <- p1+ ggtitle(sprintf("LMEM-%s ( pvalue < %s, cutoff=%s (colored fdr < 0.05)) ", plot, fdr, est)) + labs(y = "Estimate") +labs(x = NULL)
        
        print(p1)
      } else if (length(unique(resSig.sel$coefficient)) >=1){
        for (cof in unique(resSig.sel$coefficient)){
          resSig.sel.sel <- subset(resSig.sel, coefficient == cof)
          
          baseline <- unique(resSig.sel.sel$baseline)
          compare <- unique(resSig.sel.sel$coefficient);compare
          compare <- gsub(unique(resSig.sel.sel$metadata), "", compare);compare
          
          
          # color
          resSig.sel.sel$dir <- gsub('down',baseline, gsub('up',compare, resSig.sel.sel$dir))
          resSig.sel.sel$dir <- factor(resSig.sel.sel$dir, levels = c(as.character(baseline), "NS", as.character(compare)))
          
          # color end
          
          if(!is.null(mycols)){
            dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
          }else{
            dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(baseline), "NS", as.character(compare))
          }
          
          legend.labs <- 
            c(paste(baseline, " (n=", unique(resSig.sel.sel$bas.count),")",sep=""),
              paste("NS"),
              paste(compare, " (n=", unique(resSig.sel.sel$coef.count), ")",sep=""))
          
          
          
          lims <- max(abs(resSig.sel.sel$Estimate) + abs(resSig.sel.sel$SE))*1.0
          
          p1 <- ggplot(resSig.sel.sel, aes(x=reorder(ASV,Estimate), y=Estimate, color=dir)) + geom_point() +
            geom_errorbar(aes(x=ASV, ymin=Estimate-SE, max=Estimate+SE), width=0.2) + 
            geom_hline(yintercept=0) + theme_classic()  + coord_flip() +  
            scale_color_manual(values=dircolors, labels=legend.labs) + ylim(c(-lims, lims)) +scale_x_discrete(breaks = as.character(resSig.sel.sel$ASV), labels = resSig.sel.sel$taxa)+theme(plot.title = element_text(size=8, hjust = 1))
          
          #sn <- gsub("\\..*","",sn);sn
          
          p1 <- p1+ ggtitle(sprintf("LMEM  (in %s for %s pvalue < %s, cutoff=%s (colored fdr < 0.05)) ",
                                    unique(resSig.sel.sel$metadata), compare, fdr,est)) +
            labs(y = "Estimate") +labs(x = NULL)
          print(p1)
        }
      }
    }
  }
  dev.off()
}





#' A Go_deseq2_fore
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Deseq2 forest plot
#' @export
#' @examples
#' Go_deseq2_fore()
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

Go_lmem_heat <- function(project,file_path, alpha, pattern, facet, name, orders, height, width){
  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=pattern);filenames
  sample.names <- sapply(strsplit(filenames, pattern), `[`, 1) ;sample.names
  
  
  print(path)
  print(sample.names)
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else if (length(facet) == 1 & length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.%s.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, facet, name, alpha,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/10_lmem.heatmap.%s.(%s).%s.pdf", project, format(Sys.Date(), "%y%m%d"), project, alpha, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  # add input files
  df<-{}
  for (sn in sample.names) {
    file <- file.path(path, paste0(sn, pattern))
    df1 <- read.csv(file, row.names=NULL ,check.names=FALSE)
    df <- rbind(df, df1)
  }
  
  df.sel <- df
  resSig <- as.data.frame(subset(df.sel, padj < alpha)); resSig <- resSig[order(resSig$Estimate),]
  # resSig$smvar <- factor(resSig$smvar)
  
  
  resSig$dir <- ifelse(resSig$padj < 0.05, ifelse(sign(resSig$Estimate)==1, "up", "down"), "NS")
  print(1)
  for (plot in unique(resSig$metadata)){
    resSig.sel <- subset(resSig, metadata == plot)
    print(2)
    if (length(unique(resSig.sel$coefficient)) >=1 ){
      resSig.sel$comparison <-  gsub(sprintf("%s", unique(resSig.sel$metadata)) ,"" , resSig.sel$coefficient)
      resSig.sel$comparison <- factor(resSig.sel$comparison , levels = orders)
      resSig.sel$stars <- cut(resSig.sel$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
      print(3)
      if (length(facet) == 1) {
        resSig.sel[,facet] <- factor(resSig.sel[,facet] , levels = orders)
      }
      print(4)
      p <- ggplot(resSig.sel, aes(x=reorder(taxa,Estimate), y=comparison, color=comparison)) + #, alpha=padj)) +
        theme_classic()+ coord_flip() + geom_tile(aes(fill = Estimate), colour = "white") + 
        geom_text(aes(label=stars), color = "black") + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
        ggtitle(sprintf("LMEM All comparison group (%s p < %s) ", plot, alpha)) +
        theme(plot.title = element_text(hjust = 0.5))+ #0.5
        theme(legend.position= "right", axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8)) + labs(y = "comparison group") +labs(x = NULL)
      p = p + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
      
      
      if (length(facet) == 1) {
        print(5)
        ncol <- length(unique(resSig.sel[,facet]))*length(unique(resSig.sel[,"comparison"]))
        p = p + facet_wrap(as.formula(sprintf("~ %s+%s", facet, "comparison")), scales="free_x", ncol = ncol)
      }
      else {
        p = p + facet_wrap(~  comparison, scales="free_x", ncol = 10) 
      }
      
      #plotlist[[length(plotlist)+1]] <- p
      print(p)
    }
  }  
  dev.off()
}
#' A Go_pickTaxa
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity ordination plot
#' @export
#' @examples
#' Go_pickTaxa()

Go_pickTaxa <- function(psIN, project, TaxaLevel, data_type, bestPick){
  # ------------- Aggregate
  cat("#--  Getting Taxa  --#\n")
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN))) # for dada2 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame((otu_table(psIN))) # not for dada2 
  }
  
  otu.filt$Genus <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN),taxRanks=ranks, level= TaxaLevel)
  agg <- aggregate(. ~ Genus, otu.filt, sum, na.action=na.pass) #add na.action=na.pass if have error "no rows to aggregate"
  genera <- agg$Genus
  agg <- agg[,-1]
  rownames(agg) <- genera
  agg_t <-t(agg)
  dim(agg)
  
  
  # ------------- biplot Genus 선택 
  

  #-----------------Parameters----------------#
  cmethod<-"spearman" 
  #Correlation method to use: pearson, spearman, kendall
  fmethod<-"bray" 
  #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  vmethod<-"bray" 
  #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
  nmethod<-"bray" 
  #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao jaccard
  
  best5 <- bestPick #100
  
  if (data_type == "dada2" | data_type == "DADA2") {
    res.bv.step.biobio <- bv.step(wisconsin(agg_t), wisconsin(agg_t), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    res.bv.step.biobio <- bv.step(wisconsin(agg), wisconsin(agg), 
                                  fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                                  scale.fix=FALSE, scale.var=FALSE, 
                                  max.rho=0.95, min.delta.rho=0.001,
                                  random.selection=TRUE,
                                  prop.selected.var=0.3,
                                  num.restarts= best5,
                                  output.best=best5,
                                  var.always.include=NULL) 
  }
  
  taxaNames<-colnames(agg_t);taxaNames
  bestTaxaFit <- ""
  for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl))){
    bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
  }
  
  bestTaxaFit <- data.frame(bestTaxaFit);bestTaxaFit 
  colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"
  
  bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")));bio.keep
  
  taxaNames <- colnames(agg_t);taxaNames
  
  cat("\n Use this numbers for biplot;  \n")

  
  print(bio.keep)
  for (taxa in bio.keep){
    cat(sprintf("%s : %s \n", taxa, taxaNames[taxa]))
  }
  return(taxaNames)
}
#' A Go_biplot
#'

Go_biplot <- function(psIN, metaData, project, orders, distance_metrics, data_type, biplot, shapes, TaxaLevel, ID, facet, name, height, width){
  
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
  if (length(facet) >= 1) {
    if (!is.null(name)) {
      pdf(sprintf("%s/12_biplot.%s.%s.%s.%s.pdf",out_path, project,facet,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s/12_biplot.%s.%s.%s.pdf",out_path,project,facet, format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
  }
  else {
    if (!is.null(name)) {
      pdf(sprintf("%s/12_biplot.%s.%s.%s.pdf",out_path,project,name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    } 
    else {
      pdf(sprintf("%s/12_biplot.%s.%s.pdf",out_path,project,format(Sys.Date(), "%y%m%d")), height = height, width = width)
    }
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

          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))


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
          bio.fit <- envfit(ordi, agg_t[,biplot, drop=F],  perm = 999); head(bio.fit)
          df_biofit <- scores(bio.fit, display=c("vectors"))
          df_biofit <- df_biofit*vegan:::ordiArrowMul(df_biofit)
          df_biofit <- as.data.frame(df_biofit);df_biofit
          
          
          df[,facet] <- factor(df[,facet], levels = orders)
          
          df[,mvar] <- factor(df[,mvar], levels = orders)
          

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
          p = p + scale_color_manual(values = Tableau10)# scale_fill_brewer(type="qual", palette=colorset)
          p = p+geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                         arrow = arrow(length = unit(0.2, "cm")), color="#808080",alpha=0.7)+
            geom_text(data=as.data.frame(df_biofit*1.1), aes(NMDS1, NMDS2, label = rownames(df_biofit)), 
                      color="#808080",alpha=0.7, size = 3)


          p = p + 
            scale_shape_manual(values=c(16, 1, 10, 2,17,3,4,5,6,7,8,9,11,12,13,14,15,16), breaks=orders) 
          
          if (length(ID) == 1) {
            p = p +  geom_text_repel(aes_string(label = ID), size = 2)
          }else {
            p = p 
          }

          

          
          if (length(facet) == 1) {
            ncol <- length(unique(mapping.sel.na.rem[,facet]))
            p = p + facet_wrap(as.formula(sprintf("~ %s", facet)), scales="free_x", ncol = ncol)
          }else {
            p = p
          }
          
          print(p)
        }
      }
    dev.off()
}

Go_randomForest <- function(psIN, project,data_type,cutoff, numTree,numtry=NULL,name=NULL, categorical,  numeric){
  #---- install package         ------#
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  RF_out <- file.path(sprintf("%s_%s/table/RF_out",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", RF_out)) dir.create(RF_out)
  
  print(0)
  # input
  mapping <- data.frame(sample_data(psIN))
  rows<-rownames(mapping)
  mapping[,categorical] <- factor(mapping[,categorical])
  mapping <- mapping %>% mutate_all(na_if,"");mapping[,categorical]
  rownames(mapping)<-rows
  mapping.na <- mapping[!is.na(mapping[,categorical]), ] 
  mapping.na[mapping.na=="#N/A"] <- "NA"
  mapping.na[,categorical] <- factor(mapping.na[,categorical])
  
  # re construntion ps object(20201030)
  psIN <- prune_samples(rownames(mapping.na), psIN)
  

  
  
  if (length(numeric) == 1) {
    mapping.na <- mapping[!is.na(mapping[,numeric]), ]
    sample_data(psIN) <- mapping.na
  }

  
  print(1)
  # call otutable
  if (data_type == "dada2" | data_type == "DADA2") {
    otu.filt <- as.data.frame(t(otu_table(psIN)))
  } else if (data_type == "Nephele" | data_type == "nephele" | data_type == "Other" | data_type == "other") {
    otu.filt <- as.data.frame(otu_table(psIN))
  }
  
  # get taxa 
  
  if (colnames(tax_table(psIN))[1] == "KO"){
    otu.filt[,"Path"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("KO","Path"),level="Path")
    otu.filt[,"KO"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("KO","Path"),level="KO")
    otu.filt$taxa <- paste(otu.filt[,"Path"], otu.filt$KO)
  } else {
    otu.filt[,"Phylum"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Phylum")
    otu.filt[,"Genus"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Genus")
    otu.filt[,"Species"] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=c("Phylum","Genus","Species"),level="Species")
    otu.filt$taxa <- paste(otu.filt[,"Phylum"], otu.filt$Genus, otu.filt$Species)
  }
  

  
  print(2)
  otu.filt.sel <- otu.filt

  
  if (colnames(tax_table(psIN))[1] == "KO"){
    otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$KO), ]
    otu.filt.sel$Path  <- NULL
    otu.filt.sel$KO  <- NULL
  } else {
    otu.filt.sel <- otu.filt.sel[!is.na(otu.filt.sel$Genus), ]
    otu.filt.sel$Phylum  <- NULL
    otu.filt.sel$Genus  <- NULL
    otu.filt.sel$Species <- NULL
  }
  print(3)
  
  if (dim(otu.filt)[2] == 2){
    next
  }
  
  agg <- aggregate(as.formula(sprintf(". ~ %s" , "taxa")), otu.filt.sel, sum, na.action=na.pass)
  genera <- agg[,"taxa"]
  agg <- agg[,-1]
  
  rownames(agg) <- genera
  
  print(2)
  
  c=dim(agg)[1]
  d=dim(agg)[2]
  cat("===================\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_nonzero_counts <- apply(agg, 1, function(y) sum(length(which(y > 0))))
  #hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")
  
  # cutoff by percent
  remove_rare <- function( table , cutoff_pro ) {
    row2keep <- c()
    cutoff <- ceiling( cutoff_pro * ncol(table) )  
    for ( i in 1:nrow(table) ) {
      row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
      if ( row_nonzero > cutoff ) {
        row2keep <- c( row2keep , i)
      }
    }
    return( table [ row2keep , , drop=F ])
  }
  
  print(4)
  otu_table_rare_removed <- remove_rare(table=agg, cutoff_pro=cutoff)
  
  print(5)
  c=dim(otu_table_rare_removed)[1]
  d=dim(otu_table_rare_removed)[2]
  cat("===================\n")
  cat("   After removed\n")
  cat(sprintf("Input taxa number=%s\n", c))
  cat(sprintf("sample number=%s\n", d))
  cat("===================\n")
  cat("\n")
  
  otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
  otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  
  otu_table_asinh_mean_centred <- scale(asinh(agg), center=TRUE, scale=FALSE)  
  
  set.seed(151) 
  # -----  Running model ---------#
  # table
  otu_table_scaled_Categorical <- data.frame(t(otu_table_scaled))  
  # mapping
  otu_table_scaled_Categorical[,categorical] <- mapping.na[rownames(otu_table_scaled_Categorical), categorical]  


  
  #Run RF to classify inflamed and control samples:
  #numTree = 501 # default  10001등으로 늘릴수 있다.
  otu_table_scaled_Categorical.na <- otu_table_scaled_Categorical[complete.cases(otu_table_scaled_Categorical), ]

  if (numtry > 0) {

      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
  }


  
  cat("\n")
  cat("===================\n")
  cat("    RF_classify \n")
  cat("===================\n")
  print(RF_classify) 

  #Permutation Test
  RF_classify_sig <- rf.significance( x=RF_classify, xdata=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , nperm=1000 , ntree=numTree )
  
  
  #Accuracy Estimated by Cross-validation
  fit_control <- trainControl( method = "LOOCV" ) 
  
  RF_classify_loocv <- train(otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=25 ) , trControl=fit_control )
  
  print(RF_classify_loocv$results)

  
  if (length(numeric) == 1) {
    # add continuous value
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    otu_table_scaled_Numeric <- data.frame(t(otu_table_scaled))
    otu_table_scaled_Numeric[,numeric] <- as.numeric(mapping.na[rownames(otu_table_scaled_Numeric), numeric])
    otu_table_scaled_Numeric.na <- otu_table_scaled_Numeric[complete.cases(otu_table_scaled_Numeric), ]
    
    
    
    if (numtry > 0) {
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, mtry=numtry, importance=TRUE, proximities=TRUE )
    } else{
      RF_classify <- randomForest( x=otu_table_scaled_Categorical.na[,1:(ncol(otu_table_scaled_Categorical.na)-1)] , y=otu_table_scaled_Categorical.na[, ncol(otu_table_scaled_Categorical.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
      
      plot(RF_classify$err.rate[,1], type="l")
      ## NaN 찾기
      #sapply(otu_table_scaled_Numeric[,1:(ncol(otu_table_scaled_Numeric)-1)], function(x) sum(is.na(x)))
      # 수치중 1이 있으면 missing value
      #which(is.na(otu_table_scaled_Numeric$p__Firmicutes.g__Streptococcus.s__))
      #Run RF to regress OTUs against inflammation score (IS):
      RF_regress <- randomForest( x=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[ , ncol(otu_table_scaled_Numeric.na)] , ntree=numTree, importance=TRUE, proximities=TRUE )
    }
    
    cat("\n")
    cat("===================\n")
    cat("    RF_regress \n")
    cat("===================\n")
    print(RF_regress)
    
    
    RF_regress_sig <- rf.significance( x=RF_regress,  xdata=otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , nperm=1000 , ntree=numTree ) 
    
    RF_regress_loocv <- train( otu_table_scaled_Numeric.na[,1:(ncol(otu_table_scaled_Numeric.na)-1)] , y=otu_table_scaled_Numeric.na[, ncol(otu_table_scaled_Numeric.na)] , method="rf", ntree=numTree , tuneGrid=data.frame( mtry=215 ) , trControl=fit_control )
    
    print(RF_regress_loocv$results)
  }
  
  
  
  
  # save data
  if (length(numeric) == 1) {
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.%s.rda", RF_out, project,categorical,numeric, numTree,name, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_agg_regress_model.%s.%s.%s.rda", RF_out, project,categorical,numeric,numTree,name,format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.rda", RF_out, project,categorical,numeric, numTree, format(Sys.Date(), "%y%m%d")))
      saveRDS(RF_regress, file = sprintf("%s/%s.%s.%s.RF_agg_regress_model.%s.%s.rda", RF_out, project,categorical,numeric,numTree,format(Sys.Date(), "%y%m%d")))
      
      functionReturningTwoValues <- function() { 
        results <- list()
        results$classify <- RF_classify
        results$regress <-RF_regress
        return(results) 
      }
      cat("\n")
      print("RF$classify and RF$regress are returned.")
      
      functionReturningTwoValues()
    }
  }else{
    if (length(name) == 1) {
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.%s.rda", RF_out, project,categorical, numTree,name, format(Sys.Date(), "%y%m%d")))
    }else{
      saveRDS(RF_classify, file = sprintf("%s/%s.%s.%s.RF_agg_classify_model.%s.rda", RF_out, project,categorical, numTree, format(Sys.Date(), "%y%m%d")))
    }
    cat("\n")
    print("RF$classify are returned.")
    
    return(RF_classify)
  }

}


Go_confusion <- function(model, project,legend, name, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  rf_model <- readRDS(model)
  
  # run
  conf_mat <- rf_model$confusion[, -ncol(rf_model$confusion)]
  melt_conf_mat <- reshape2::melt(conf_mat, na.rm = TRUE)
  table <- data.frame(melt_conf_mat);table
  
  plotTable <- table %>%
    mutate(goodbad = ifelse(table$Var1 == table$Var2, "good", "bad")) %>%
    group_by(Var2) %>%
    mutate(prop = value/sum(value))
  
  
  # plot
  p <- ggplot(data = plotTable, mapping = aes(x = Var2, y = Var1, fill = goodbad, alpha = prop)) +
    geom_tile() +
    geom_text(aes(label = value), vjust = .5, fontface  = "bold", alpha = 1) +
    scale_fill_manual(values = c(good = "blue", bad = "red")) +
    theme_bw() + xlab("True class") + ylab("Predicted class") + theme(legend.position=legend, plot.title = element_text(hjust = 0.5))
    #xlim(rev(levels(table$Var1))) 

  p <- p+ ggtitle(sprintf("confusionMatrix(%s) ",name)) 
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_2.confusion.%s.%s.%s.pdf",out_path, project, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_2.confusion.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}






Go_roc <- function(model, project,map, categorical, name, height, width){
    
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  if (length(name) == 1) {
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.%s.pdf", out_path, name,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }else{
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.pdf", out_path,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }
  
  rf_model <- readRDS(model)
  
  
  pred <- predict(rf_model, type="prob")
  
  # map 정리 
  sel <- intersect(rownames(map), rownames(pred))
  map <- map[sel,, drop=F]
  #dim(pred);dim(map); class(pred);class(map);
  
  
  #dim(map.sel.sel)
  par(mfrow = c(1,length(colnames(pred))))
  for (i in 1:length(colnames(pred))){
    print(i)
    pred2 <- prediction(pred[,i], as.matrix(map[,categorical]))# or oedered(map[,categorical])
    perf <- performance(pred2, "tpr", "fpr")
    perf.auc <- performance(pred2, "auc")
    pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true = map[,categorical], stringsAsFactors=F); 
    pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
    confusion_matrix <- table(pred_df[, c("true", "predicted")])
    accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
    vec.pred <- as.numeric(pred_df$predicted)-2; vec.true <- as.numeric(pred_df$true)-2 # it was "-1"
    mccvalue <- mcc(vec.pred, vec.true)
    if (i == 1){
      color = "red"
    }else if(i == 2){
      color = "blue"
    }else if(i == 3){
      color = "green"
    }
    
    if (length(name) == 1) {
      plot(perf, main=sprintf("RF_ROC %s (%s)", colnames(pred)[i], name), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }else{
      plot(perf, main=sprintf("RF_ROC %s", colnames(pred)[i]), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }
  }
  dev.off()
}





Go_importance_plot <- function(psIN, model, project,title,aggregate, MDA, name, bySample, height, width){
    
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  RF_model <- readRDS(model)
  tmp <- data.frame(RF_model$importance)
  
  if (aggregate == "YES" | aggregate == "Yes"| aggregate == "yes") {
    tmp$features <- rownames(tmp)
    tmp.sorted <- arrange(tmp, desc(MeanDecreaseAccuracy))
    RF.sorted.sel <- as.data.frame(subset(tmp.sorted, MeanDecreaseAccuracy > MDA))
    RF.sorted.sel$ShortName <- RF.sorted.sel$features
    
  } else if (aggregate == "NO" | aggregate == "No" | aggregate == "no") {
    RF <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF[,taxa] == "NA"
      RF[,taxa]<- as.character(RF[,taxa])
      RF[,taxa][is.na(RF[,taxa])] <- "__"
      for(i in 1:length(RF[,taxa])){
        if (RF[,taxa][i] == "s__" || RF[,taxa][i] == "g__" || RF[,taxa][i] == "f__" || RF[,taxa][i] == "o__" || RF[,taxa][i] == "c__"|| RF[,taxa][i] == "__"){
          RF[,taxa][i] <- ""
        }
      }
    }
    
    RF$ShortName <- paste(RF$Genus,"",RF$Species)
    for(taxa in c("Family", "Order", "Class","Phylum")){
      for(i in 1:length(RF[,taxa])){
        if (RF$ShortName[i] != "   "){
          next
        }      else if (RF$ShortName[i] == "   " & RF[,taxa][i] != ""){
          RF$ShortName[i] <- paste(RF[,taxa][i])
        }
      }
    }
    RF$features <- rownames(RF)
    RF.sorted <- arrange(RF, desc(MeanDecreaseAccuracy)  )
    RF.sorted.sel <- as.data.frame(subset(RF.sorted, MeanDecreaseAccuracy > MDA))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF.sorted.sel[,taxa] = NULL

    }
  }
  
  RF.sorted.sel.melt <- melt(RF.sorted.sel, id.vars=c("features","MeanDecreaseAccuracy","MeanDecreaseGini","ShortName"))
  

  if (bySample == "NO" |bySample == "No" | bySample == "no"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  } else if (bySample == "Yes" | bySample == "Yes"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value, fill=variable))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  }


  
  p<- p+ ggtitle(sprintf("RF %s (MDA>%s) ", title, MDA)) 
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_1.RF.%s.%s.(%s).%s.pdf",out_path, project, name, MDA,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_1.RF.%s.(%s).%s.pdf", out_path, project, MDA, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}





#' A Go_rpkm
#'
#' getting RPKM table
#' @param love getting RPKM table
#' @keywords rpkm
#' @export
#' @examples
#' Go_rpkm


Go_rpkmToPs <- function(project, speciesTab, genoemSizeTab, taxaTab) {
  # out path
  tem_path <- file.path(sprintf("RPKM_tem_%s", format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", tem_path)) dir.create(tem_path)
  rds_path <- file.path(sprintf("RPKM_rds_%s", format(Sys.Date(), "%y%m%d")))
  if(!file_test("-d", rds_path)) dir.create(rds_path)
  

  bacTab <- speciesTab
  
  #-- get simple species and genus names --#  species 이름이 길다.
  bacTab$Species <- rownames(bacTab)
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(bacTab$Species),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab$Genus <- cleanSpecies$X1
  bacTab$Species <- cleanSpecies$Species
  bacTab$SubSpecies <- rownames(bacTab)
  
  # read genome size table
  gSizeTab <- genoemSizeTab
  colnames(gSizeTab) <- c("AccessionVersion", "TaxId", "Completeness", "GenomeLen", "Title")
  cleanSpecies_gSizeTab <- data.frame(do.call('rbind', strsplit(as.character(gSizeTab$Title),' ',fixed=TRUE)))
  cleanSpecies_gSizeTab$Species <- paste(cleanSpecies_gSizeTab$X1, cleanSpecies_gSizeTab$X2, sep="_")
  cleanSpecies_gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Genus <- cleanSpecies_gSizeTab$X1
  gSizeTab$Species <- cleanSpecies_gSizeTab$Species
  
  
  # create new data frame and exrtact genome size for calculation
  genomeFinal <- data.frame(cbind(bacTab$Species, bacTab$Genus)) # 두번째 col은 이름이 안따라온다. 다음 명령어를 사용해야 한다.
  colnames(genomeFinal) <- c("Species","Genus")
  genomeFinal$Genus <- bacTab$Genus
  genomeFinal$GenomeLen <- gSizeTab$GenomeLen[match(bacTab$Genus, gSizeTab$Genus)] # index match
  
  
  write.csv(genomeFinal, quote = FALSE, #col.names = NA, #row.names = FALSE,
            file=sprintf("%s/%s.genomeFinal2.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  bacTab$Genus <- NULL
  bacTab$Species <- NULL
  bacTab$SubSpecies <- NULL
  ## RPKM
  geneLength <- 1
  bacNor <- data.frame(sapply(bacTab , function(column) 10^9 * column / geneLength / sum(column)))
  rownames(bacNor) <- rownames(bacTab); head(rownames(bacNor))
  
  #brac.rpkmS <- bracNor*100/genomeFinal$SlenS
  #brac.rpkmS <- ceiling(brac.rpkmS[-c(99),]) # 반올림
  #brac.rpkmS[is.na(brac.rpkmS)] <- 0 # remove NA
  
  bacTab.rpkmG <- bacNor*100/genomeFinal$GenomeLen
  bacTab.rpkmG <- ceiling(bacTab.rpkmG[-c(99),]) # 반올림
  bacTab.rpkmG[is.na(bacTab.rpkmG)] <- 0 # remove NA
  
  ## add LCA
  cleanSpecies <- data.frame(do.call('rbind', strsplit(as.character(rownames(bacTab.rpkmG)),'_',fixed=TRUE)))
  cleanSpecies$Species <- paste(cleanSpecies$X1, cleanSpecies$X2, sep="_")
  bacTab.rpkmG$Species <- cleanSpecies$Species
  bacTab.rpkmG$Genus <- cleanSpecies$X1
  
  Taxa <- taxaTab
  
  ranks <- c("Rank1", "Phylum", "Class", "Order", "Family")
  for(taxa in ranks){
    bacTab.rpkmG[,taxa]  <- Taxa[,taxa][match(bacTab.rpkmG$Genus, Taxa$Genus)]
  }
  
  bacTab.rpkmG.sel <- subset(bacTab.rpkmG, Rank1 == "Bacteria")
  colnames(bacTab.rpkmG.sel) <-  gsub("X", "", colnames(bacTab.rpkmG.sel));head(colnames(bacTab.rpkmG.sel))
  
  ##########################
  #---      Make ps     ---#
  ##########################
  #------------- tax table -------------#
  headers <- rownames(bacTab.rpkmG.sel);head(headers)
  
  tax <- data.frame(bacTab.rpkmG.sel$Rank1, bacTab.rpkmG.sel$Phylum, bacTab.rpkmG.sel$Class, bacTab.rpkmG.sel$Order, bacTab.rpkmG.sel$Family, bacTab.rpkmG.sel$Genus, bacTab.rpkmG.sel$Species)
  rownames(tax) <- headers
  colnames(tax) <- c("Kingdom","Phylum","Class", "Order", "Family", "Genus", "Species") ;tax
  
  #------------- otu table -------------#
  for(i in colnames(tax)){
    bacTab.rpkmG.sel$Rank1 <- NULL
    bacTab.rpkmG.sel[,i] <- NULL
  }
  
  
  otu <- bacTab.rpkmG.sel
  rownames(otu) <- headers
  otu <- otu[rowSums(otu[, -1])>0, ]
  
  write.csv(otu, quote = FALSE, col.names = NA, #row.names = FALSE,
            file=sprintf("%s/otu_%s.genome_size.%s.csv", tem_path, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  #--- create phyloseq file ---#
  tax<-as.matrix(tax)
  otu<-as.matrix(otu)
  
  OTU <- otu_table(otu, taxa_are_rows = TRUE);head(OTU)
  TAX <- tax_table(tax);head(TAX)
  ps1 <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), tax_table(TAX))
  saveRDS(ps1, sprintf("%s/ps_rpkmG.%s.%s.rds",rds_path, project, format(Sys.Date(), "%y%m%d")))
  
  return(ps1)
}





#' A Go_correlation
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_correlation()
#' May 22 2020



Go_correlation <- function(project, metaData, metabolicTab, abundTab, map, method, xanlgle, ncol,name, orders, height, width){
  # output
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  # map 정리 2
  sel.meta <- intersect(rownames(metabolicTab), rownames(map)); head(sel, "3")
  metabolicTab.sel <- metabolicTab[sel.meta,, drop=F];dim(metabolicTab.sel)
  
  sel.abun <- intersect(rownames(abundTab), rownames(map)); head(sel, "3")
  abundTab.sel <- abundTab[sel.abun,, drop=F];dim(abundTab.sel)
  
  
  
  x <- log((abundTab.sel+1)/(rowSums(abundTab.sel)+dim(abundTab.sel)[2]))
  x <- x[,order(colSums(x),decreasing=TRUE)]
  y <- metabolicTab.sel
  

  
  if(length(name) == 1){
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.%s.pdf",out_path,project, method,name,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }else{
    pdf(sprintf("%s/12_Correlation.%s.%s.%s.pdf",out_path,project, method,format(Sys.Date(), "%y%m%d")), height = height, width=width)
  }
  
  
  for (des in rownames(subset(metadata, Go_correlation=="yes"))){
    groups<-map[,des]
    #Now calculate the correlation between individual Taxa and the environmental data
    df<-NULL
    for(i in colnames(x)){
      for(j in colnames(y)){
        for(k in unique(groups)){
          a <- x[groups==k,i,drop=F]
          b <- y[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
          if(is.null(df)){
            df<-tmp  
          }
          else{
            df<-rbind(df,tmp)
          }    
        }
      }
    }
    
    df<-data.frame(row.names=NULL,df)
    colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
    df$Pvalue<-as.numeric(as.character(df$Pvalue))
    df$AdjPvalue<-rep(0,dim(df)[1])
    df$Correlation<-as.numeric(as.character(df$Correlation))
    
    adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
    adjustment<-5
    
    if(adjustment==1){
      df$AdjPvalue<-df$Pvalue
    } else if (adjustment==2){
      for(i in unique(df$Env)){
        for(j in unique(df$Type)){
          sel<-df$Env==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==3){
      for(i in unique(df$Taxa)){
        for(j in unique(df$Type)){
          sel<-df$Taxa==i & df$Type==j
          df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
        }
      }
    } else if (adjustment==4){
      for(i in unique(df$Taxa)){
        sel<-df$Taxa==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    } else if (adjustment==5){
      for(i in unique(df$Env)){
        sel<-df$Env==i
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
      }
    }
    
    #Now we generate the labels for signifant values
    
    df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.01, 0.05, 0.08, Inf), label=c("***", "**", "*", ""))
    
    #We ignore NAs
    df<-df[complete.cases(df),]
    
    #We want to reorganize the Env data based on they appear
    #df$Env<-factor(df$Env,as.character(df$Env))
    
    #We use the function to change the labels for facet_grid in ggplot2
    Env_labeller <- function(variable,value){
      return(sel_env_label[as.character(value),"Trans"])
    }
    
    df$Env <- factor(df$Env)
    df$Type <- factor(df$Type, levels = orders)
    p <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
    p <- p + geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") 
    p<-p+theme(axis.text.x = element_text(angle = xanlgle, hjust = 1, vjust=0.5))
    p<-p+geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
    p<-p+facet_wrap (~  Env, scales="free_x", ncol = ncol)
    #p<-p+facet_wrap (~  Env, ncol = 11)
    p<- p+ ggtitle(sprintf("%s %s", des, "Ladas")) 
    print(p)
  }
  dev.off()
}

#' A Go_huamnn2ps
#'
#' 
#' @param huamnn2ps 
#' @keywords huamnn2ps
#' @export
#' @examples
#' Go_huamnn2ps




Go_hmPath2ps<-function(project, tool, pathway, name, alpha){ # alpha, for rowsum
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  print("version1")
  
  # input files kegg-orthology stratified
  pt <- read.csv(sprintf("%s",pathway), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
  
  #ko.sta <- round(1000*ko.sta)
  
  # remove na column
  all_na <- function(x) any(!is.na(x))
  pt.na <- pt %>% select_if(all_na)
  print(sprintf("remove column na %s to %s",dim(pt)[2], dim(pt.na)[2]))
 
  ## rowsum 
  # remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance
  # inds_to_remove <- which(rowSums(ko.sta.na) < alpha) 
  # ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  # print(dim(ko.sta.na))
  
  
  
  ## create abundance tab
  ## 각 수치가 > 2 10개 미만이면 빼기
   n <- round(length(colnames(pt.na))/10);n
   inds_to_remove <- which(rowSums(pt.na > alpha) < n) 
   pt.na <- pt.na[setdiff(1:nrow(pt.na), inds_to_remove),];dim(pt.na)
  
  print(dim(pt.na))
  # remove less important
  for (rev in c("UNINTEGRATED","UNMAPPED")){
    pt.na <- subset(pt.na, !(grepl(rev, rownames(pt.na))))
    print(dim(pt.na))
  }
  
  
  ## create pathway tab
  if (tool == "humann2"){
    # split names of kolist
    head(rownames(pt.na))
    rownames(pt.na) <- gsub(",", ".", rownames(pt.na));head(rownames(pt.na))
    pt.na.list <- data.frame(str_split(rownames(pt.na), ": ", simplify = T));head(pt.na.list)
    

    # add header names
    headers <- vector(dim(pt.na.list)[1], mode="character")
    for (i in 1:dim(pt.na.list)[1]) {
      headers[i] <- paste("path", i, sep="_")
    }
    
    rownames(pt.na.list) <- headers
    colnames(pt.na.list) <- c("Path","Path.des");head(pt.na.list)
    pt.na.list$Path.des <- gsub(" ", "_", pt.na.list$Path.des)
    
    rownames(pt.na) <- headers;head(pt.na)
    
    
  }else if(tool == "picrust2"){
    #koList.sta1 <- data.frame(rownames(ko.sta.na))
    #koList.sta1$KO <- rownames(ko.sta.na);head(koList.sta1)
    #koList.sta1$KO.des <- factor(koTOpath$KO.description[match(koList.sta1$KO, koTOpath$KO)]);head(koList.sta1)
    #koList.sta1$Path <- factor(koTOpath$Path[match(rownames(ko.sta.na), koTOpath$KO)]);head(koList.sta1)
    #koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
    
    # add header names
    #rownames(koList.sta1) <- koList.sta1$rownames.ko.sta.na.
    #koList.sta1$rownames.ko.sta.na. <-NULL
    #colnames(koList.sta1) <- c("KO", "KO.des","Path","Path.des");head(koList.sta1)
  }
  

  
  # save information
  write.csv(pt.na.list, quote = F, col.names = F, row.names=T,
            file=sprintf("%s/pt.na.list.%s.csv", out, format(Sys.Date(), "%y%m%d"),sep="\t"))
  
  
  #--- create phyloseq file ---#
  
  tab <- as.matrix(pt.na)
  list <- as.matrix(pt.na.list)
  
  
  TAB <- otu_table(tab, taxa_are_rows = TRUE);dim(TAB)
  LIST <- tax_table(list);dim(LIST)



  #ps1
  ps.path <- phyloseq(otu_table(TAB, taxa_are_rows=FALSE), tax_table(LIST))
  
  print(ps.path)
  if (length(name) == 1) {
    saveRDS(ps.path, sprintf("%s/ps.pathway.%s.%s.%s.rds", out, project, name, format(Sys.Date(), "%y%m%d")))
  } else{
    saveRDS(ps.path, sprintf("%s/ps.pathway.%s.%s.rds", out, project,format(Sys.Date(), "%y%m%d")))
  }
}#' A Go_huamnn2ps
#'
#' 
#' @param huamnn2ps 
#' @keywords huamnn2ps
#' @export
#' @examples
#' Go_huamnn2ps

Go_huamnn2ps<-function(project, tool, koTOpath, koSta, name, alpha){ # alpha, for rowsum
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  print("version2")
  
  # input file for index match file koTOpath
  koTOpath <- read.csv(sprintf("%s",koTOpath),header=T,as.is=T,row.names=1,check.names=F)
  
  # input files kegg-orthology stratified
  if (tool == "humann2"){
    ko.sta <- read.csv(sprintf("%s",koSta), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
  }else if(tool == "picrust2"){
    ko.sta <- read.csv(sprintf("%s",koSta),header=T,as.is=T,row.names=1,check.names=F)
  }
  

  
  
  
  #ko.sta <- round(1000*ko.sta)
  
  # remove na column
  all_na <- function(x) any(!is.na(x))
  ko.sta.na <- ko.sta %>% select_if(all_na)
  print(sprintf("remove column na %s to %s",dim(ko.sta)[2], dim(ko.sta.na)[2]))
 
  ## rowsum 
  # remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance
  # inds_to_remove <- which(rowSums(ko.sta.na) < alpha) 
  # ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  # print(dim(ko.sta.na))
  
  
  ## 각 수치가 > 2 10개 미만이면 빼기
   n <- round(length(colnames(ko.sta.na))/10);n
   inds_to_remove <- which(rowSums(ko.sta.na > alpha) < n) 
   ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  
  print(dim(ko.sta.na))
  # remove less important
  for (rev in c("NO_NAME", "ribosomal protein","UNGROUPED","unclassified","UNMAPPED")){
    ko.sta.na <- subset(ko.sta.na, !(grepl(rev, rownames(ko.sta.na))))
    print(dim(ko.sta.na))
  }
  
  

  if (tool == "humann2"){
    # split names of kolist
    head(rownames(ko.sta.na))
    rownames(ko.sta.na) <- gsub(",", ".", rownames(ko.sta.na));head(rownames(ko.sta.na))
    koList.sta1 <- data.frame(str_split(rownames(ko.sta.na), ": ", simplify = T));head(koList.sta1)
    koList.sta1$X3 <- NULL; head(koList.sta1)
    
    # koList.sta2 <- data.frame(str_split(koList.sta1$X2, "\\|", simplify = T));head(koList.sta2)
    # koList.sta3 <- data.frame(str_split(koList.sta2$X2, "\\.", simplify = T));head(koList.sta3)
    
    # koList.sta1$X2 <- koList.sta2$X1; head(koList.sta1)
    # koList.sta1$X3 <- as.character(koList.sta2$X2); head(koList.sta1)
    # koList.sta1$X3 <- koList.sta3$X1; head(koList.sta1)
    # koList.sta1$X4 <- koList.sta3$X2; head(koList.sta1)
    
    # Add path information
    koList.sta1$Path <- factor(koTOpath$Path[match(koList.sta1$X1, koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
    
    
    # add header names
    headers <- vector(dim(koList.sta1)[1], mode="character")
    for (i in 1:dim(koList.sta1)[1]) {
      headers[i] <- paste("func", i, sep="_")
    }
    
    rownames(koList.sta1) <- headers
    colnames(koList.sta1) <- c("KO", "KO.des","Path","Path.des");head(koList.sta1)
    koList.sta1$funcNo <- headers;head(koList.sta1)
    rownames(ko.sta.na) <- headers #;head(ko.sta)
    
  }else if(tool == "picrust2"){
    koList.sta1 <- data.frame(rownames(ko.sta.na))
    koList.sta1$KO <- rownames(ko.sta.na);head(koList.sta1)
    koList.sta1$KO.des <- factor(koTOpath$KO.description[match(koList.sta1$KO, koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path <- factor(koTOpath$Path[match(rownames(ko.sta.na), koTOpath$KO)]);head(koList.sta1)
    koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
    
    # add header names
    rownames(koList.sta1) <- koList.sta1$rownames.ko.sta.na.
    koList.sta1$rownames.ko.sta.na. <-NULL
    #colnames(koList.sta1) <- c("KO", "KO.des","Path","Path.des");head(koList.sta1)
  }
  

  
  # save information
  write.csv(koList.sta1, quote = F, col.names = F, row.names=T,
            file=sprintf("%s/koList.sta.%s.csv", out, format(Sys.Date(), "%y%m%d"),sep="\t"))
  
  # column 정리 
  #colnames(ko.sta)

  
  #--- create phyloseq file ---#
  
  otu.sta <- as.matrix(ko.sta.na)
  tax.sta <- as.matrix(koList.sta1)
  
  
  OTU.sta <- otu_table(otu.sta, taxa_are_rows = TRUE);dim(OTU.sta)
  TAX.sta <- tax_table(tax.sta);dim(TAX.sta)



  #ps1
  kps1.sta <- phyloseq(otu_table(OTU.sta, taxa_are_rows=FALSE), tax_table(TAX.sta))
  
  print(kps1.sta)
  if (length(name) == 1) {
    saveRDS(kps1.sta, sprintf("%s/ps.funtion.%s.%s.%s.rds", out, project, name, format(Sys.Date(), "%y%m%d")))
  } else{
    saveRDS(kps1.sta, sprintf("%s/ps.funtion.%s.%s.rds", out, project,format(Sys.Date(), "%y%m%d")))
  }
}#' A Go_huamnn2ps
#'
#' 
#' @param huamnn2ps 
#' @keywords huamnn2ps
#' @export
#' @examples
#' Go_huamnn2ps

Go_huamnn2ps_stratified<-function(project, koTOpath, koSta,alpha, name){
  # out dir
  out <- file.path("2_rds") 
  if(!file_test("-d", out)) dir.create(out)
  
  # input file for index match file koTOpath
  koTOpath <- read.csv(sprintf("%s",koTOpath),header=T,as.is=T,row.names=1,check.names=F)
  
  # input files kegg-orthology stratified
  ko.sta <- read.csv(sprintf("%s",koSta), header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")
  
  #ko.sta <- round(1000*ko.sta)
  
  # remove na column
  all_na <- function(x) any(!is.na(x))
  ko.sta.na <- ko.sta %>% select_if(all_na)
  print(sprintf("remove column na %s to %s",dim(ko.sta)[2], dim(ko.sta.na)[2]))

  # remove pathways with <50 reads, detected in less than 10 samples, scale to relative abundance

  
  inds_to_remove <- which(rowSums(ko.sta.na) < alpha) 
  ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  print(dim(ko.sta.na))
  
  n <- round(length(colnames(ko.sta.na))/10)
  inds_to_remove <- which(rowSums(ko.sta.na > alpha) < n) # 각 수치가 > 2 10개 미만이면 빼기
  ko.sta.na <- ko.sta.na[setdiff(1:nrow(ko.sta.na), inds_to_remove),]
  
  print(dim(ko.sta.na))
  
  # remove less important #,"unclassified"
  for (rev in c("NO_NAME", "ribosomal protein","UNGROUPED")){
    ko.sta.na <- subset(ko.sta.na, !(grepl(rev, rownames(ko.sta.na))))
    print(dim(ko.sta.na))
  }
  
  
#view(head(ko.sta.na))
  
  # split names of kolist
  head(rownames(ko.sta.na))
  rownames(ko.sta.na) <- gsub(",", ".", rownames(ko.sta.na));head(rownames(ko.sta.na))
  koList.sta1 <- data.frame(str_split(rownames(ko.sta.na), ": ", simplify = T));head(koList.sta1)
  koList.sta2 <- data.frame(str_split(koList.sta1$X2, "\\|", simplify = T));head(koList.sta2)
  koList.sta3 <- data.frame(str_split(koList.sta2$X2, "\\.", simplify = T));head(koList.sta3)
  
  koList.sta1$X2 <- koList.sta2$X1; head(koList.sta1)
  koList.sta1$X3 <- as.character(koList.sta2$X2); head(koList.sta1)
  koList.sta1$X3 <- koList.sta3$X1; head(koList.sta1)
  koList.sta1$X4 <- koList.sta3$X2; head(koList.sta1)
  
  # Add path information
  koList.sta1$Path <- factor(koTOpath$Path[match(koList.sta1$X1, koTOpath$KO)]);head(koList.sta1)
  koList.sta1$Path.des <- factor(koTOpath$Path.description[match(koList.sta1$Path, koTOpath$Path)]);head(koList.sta1)
  
  # add header names
  headers <- vector(dim(koList.sta1)[1], mode="character")
  for (i in 1:dim(koList.sta1)[1]) {
    headers[i] <- paste("func", i, sep="_")
  }
  
  rownames(koList.sta1) <- headers
  colnames(koList.sta1) <- c("KO", "KO.des","Genus","Species","Path","Path.des");head(koList.sta1)
  koList.sta1$funcNo <- headers;head(koList.sta1)
  
  # save information
  write.csv(koList.sta1, quote = F, col.names = F, row.names=T,
            file=sprintf("%s/koList.sta.%s.csv", out, format(Sys.Date(), "%y%m%d"),sep="\t"))
  
  # column 정리 
  #colnames(ko.sta)
  rownames(ko.sta.na) <- headers#;head(ko.sta)
  
  #--- create phyloseq file ---#
  tax.sta <- as.matrix(koList.sta1)
  otu.sta <- as.matrix(ko.sta.na)
  
  OTU.sta <- otu_table(otu.sta, taxa_are_rows = TRUE);head(OTU.sta)
  TAX.sta <- tax_table(tax.sta);head(TAX.sta)
  
  #ps1
  kps1.sta <- phyloseq(otu_table(OTU.sta, taxa_are_rows=FALSE), tax_table(TAX.sta))
  print(kps1.sta)
  if (length(name) == 1) {
    saveRDS(kps1.sta, sprintf("%s/ps.sta.funtion.%s.%s.%s.rds", out, project, name, format(Sys.Date(), "%y%m%d")))
  } else{
    saveRDS(kps1.sta, sprintf("%s/ps.sta.funtion.%s.%s.rds", out, project,format(Sys.Date(), "%y%m%d")))
  }
}
#' A Go_contam
#'
#'
#' @param Go_contam
#' @keywords Go_contam
#' @export
#' @examples
#' Go_contam


Go_contam <- function(psIN, project,
                      name=NULL,
                      contaminant=NULL,
                      neg=NULL,
                      taxa=NULL) {


  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)

  
  map <- data.frame(sample_data(psIN)) 

  #====== step 1 define contaminated ps and ucontaminated ps
  if(!is.null(neg)){

    contams <- psIN
    map.contams <-  subset(map, map[,contaminant] == neg)
    sample_data(contams) <- map.contams

    samples <- psIN
    
    map.samples <-  subset(map, map[,contaminant] != neg)
    sample_data(samples) <- map.samples

    
    #contams <- subset_samples(psIN, map[,contaminant] == neg);contams
    #samples <- subset_samples(psIN, map[,contaminant] != neg);samples
    #====== step 2 substrate ucontaminated ps -  contaminated ps
    contams.n1 <- prune_taxa(taxa_sums(contams)>=1, contams) 
    samples.n1 <- prune_taxa(taxa_sums(samples)>=1, samples) 
    allTaxa <- names(sort(taxa_sums(psIN),TRUE))
    negtaxa <- names(sort(taxa_sums(contams.n1),TRUE))
    taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
    ps.decontam <- prune_taxa(taxa.noneg,samples.n1)
    print("step2")
  }else{
    ps.decontam <- psIN
  }



  if(!is.null(taxa)){
    seqtab.nochim <- as.matrix(otu_table(ps.decontam))
    tax <- as.matrix(tax_table(ps.decontam))
    is.a <- tax[,"Species"] %in% taxa 
    seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
    tax.a <- tax[!is.a,]

    ps.decontam <- phyloseq(otu_table(seqtab.a, taxa_are_rows=FALSE), tax_table(tax.a));ps.decontam
    sample_names(ps.decontam)
    sample_names(ps.decontam) <- gsub("X","",sample_names(ps.decontam));sample_names(ps.decontam)
  }else{
    ps.decontam <- ps.decontam
  }


  #====== step 3 get decontam.fna
  seqs <- getSequences(otu_table(ps.decontam))
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  print("step3")
  write(fasta, file=sprintf("%s/%s%s.%s.decontam.fna",out, 
                            ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                            project, 
                            format(Sys.Date(), "%y%m%d"),sep="/"))
  
  #====== step 4 get the table
  otu <- as.data.frame(t(otu_table(ps.decontam)));dim(otu)
  tax <- tax_table(ps.decontam);dim(tax)
  
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s%s.%s.asvTable_decontam.csv",out,
                         ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                         project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  saveRDS(ps.decontam, sprintf("%s/ps_decontam.%s%s.%s.rds", rds, 
                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                               project,format(Sys.Date(), "%y%m%d")))
  
  
  print(psIN)
  print(ps.decontam)  
  
  return(ps.decontam)
}

Go_deseq2fishtaco <- function(psIN, project, file_path){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
  out_fishtaco <-file.path(sprintf("%s_%s/table/deseq2/funTabforFishtaco",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_fishtaco)) dir.create(out_fishtaco)
  
  
  
  
  path <- file_path
  print(path)
  filenames <- list.files(path, pattern="Forfishtaco.csv")
  print(filenames)
  
  for(f in filenames){
    top.deseq2 <- read.csv(sprintf("%s/%s",path,f),row.names=1,check.names=FALSE)
    print(head(top.deseq2))
    ranks <-c("KO", "KO.des","Genus","Species","Path","Path.des")
    # otu table
    otu.filt <- as.data.frame(otu_table(psIN)) 
    
    otu.filt$ko <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=ranks,  level="KO")
    agg <- aggregate(. ~ ko, otu.filt, sum)
    
    
    agg.clean <- agg[agg$ko %in% top.deseq2$x, ]# agg[!agg$ko %in% top.deseq2$x, ]  top.deseq2$x만 제거
    
    write.csv(agg.clean, quote = FALSE,col.names = NA,file=sprintf("%s/%s",out_fishtaco,f,sep="/"))
  }
}







# The ASVs were aligned the sequences with DECIPHER using AlignSeqs or AlignTranslation and then cluster them into OTUs at a specific similarity level with DistanceMatrix and IdClusters.


Go_cluster <- function(psIN, project,db, percent){
  library(tibble)
  library(dplyr)
  # install.packages("remotes")
  # remotes::install_github("mikemc/speedyseq")
  library(speedyseq)
  # Packages that are required but not loaded:
  # library(DECIPHER)
  # library(Biostrings)
  
  # out dir
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  
  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  # ----- Input ------#
  project.name <-sprintf("%s_%s",project,percent);project.name
  
  ps <- psIN
  x <- 1-percent/100;x
  
  seqtab <- otu_table(psIN)
  
  
  nproc <- 4 # set to number of cpus/processors to use for the clustering
  
  asv_sequences <- colnames(seqtab);head(asv_sequences)
  sample_names <- rownames(seqtab);(sample_names)
  dna <- Biostrings::DNAStringSet(asv_sequences)
  
  
  ## Find clusters of ASVs to form the new OTUs
  aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
  d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
  clusters <- DECIPHER::IdClusters(
    d, 
    method = "complete",
    cutoff = x , # use `cutoff = 0.03` for a 97% OTU 
    processors = nproc)
  
  ## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
  # prep by adding sequences to the `clusters` data frame
  cluster <- clusters %>%
    add_column(sequence = asv_sequences)
  
  merged_seqtab <- seqtab %>% 
    t %>%
    rowsum(clusters$cluster) %>%
    t
  
  
  # rebuilt ASVs table 
  clustered <- distinct(cluster, cluster, .keep_all = TRUE)

  merged_seqtab.t <- data.frame(t(merged_seqtab))
  merged_seqtab.t$seqs <- factor(clustered$sequence[match(as.factor(rownames(merged_seqtab.t)), as.factor(clustered$cluster))])
  
  rownames(merged_seqtab.t) <- merged_seqtab.t$seqs
  merged_seqtab.t$seqs <- NULL
  
  seqtab <- as.matrix(t(merged_seqtab.t))
  
  #----- save seqs.fna for tree  -----#
  seqs <- getSequences(seqtab)
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  write(fasta, file=sprintf("%s/ASVs%s.%s.%s.seqs.fna",  out, percent,project, format(Sys.Date(), "%y%m%d"),sep="/"))
  

  #----- re classification  -----#
  tax <- assignTaxonomy(seqtab, db, 
                        taxLevels = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"), 
                        minBoot = 80, verbose = TRUE, multithread = TRUE)
  
  
  #----- merge  -----# 
  ps_clustered <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax));ps_clustered
  
  #----- setting fix species names  -----# 
  
  tax <- data.frame(tax_table(ps_clustered))
  tax$Species.1 <- paste(tax$Genus,tax$Species)
  
  tax$Species <- NULL
  colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_table(ps_clustered) <- as.matrix(tax)
  
  
  
  
  cat("#--  Before clustered  --#\n")
  print(psIN)
  cat("\n")
  cat(sprintf("#--  After clustered by %s   --#\n",percent))
  print(ps_clustered)
  
  saveRDS(ps_clustered, sprintf("%s/ps%s_filtered.%s.%s.rds", rds, percent, project.name,format(Sys.Date(), "%y%m%d")))
  
  otu <- as.data.frame(t(otu_table(ps_clustered)));dim(otu)
  tax <- tax_table(ps_clustered);dim(tax)
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,row.names = T, file=sprintf("%s/ASVs%s_clustered.%s.%s.csv", out, percent, project.name,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  
}



# R source file

#-- to make color table---#
# 190529
# fishtaco에서 색깔을 가져 왔다.
# Arsenic 분석을 예시로 하였고, 작동 된다.
# 200524 function으로 제작 

#-------- give color to phylum --------#
# https://alloyui.com/examples/color-picker/hsv.html


Go_color <- function(cdf, taxaName){
  cPalette <-cdf
  num_of_final_phyla = length(unique(cdf$PhylumCol)); num_of_final_phyla
  cPalette$h = 0
  cPalette$s = 0
  cPalette$v = 0

  col1 <- "Actinobacteriota" # Actinobacteriota/ Actinobacteria
  col2 <- "Bacteroidota" # Bacteroidota/ Bacteroidetes
  col3 <- "Firmicutes" # Firmicutes/ Firmicutes
  col4 <- "Fusobacteriota"# Fusobacteriota/ Fusobacteria
  col5 <- "Proteobacteria" # Proteobacteria/Proteobacteria
  col6 <- "Verrucomicrobia" # Verrucomicrobia
  col7 <- "TM7" # Verrucomicrobia
  col8 <- ""   # Patescibacteria/

  num_col1 = length(grep(col1, cPalette$PhylumCol));num_col1
  num_col2 = length(grep(col2, cPalette$PhylumCol));num_col2
  num_col3 = length(grep(col3, cPalette$PhylumCol));num_col3
  num_col4 = length(grep(col4, cPalette$PhylumCol));num_col4
  num_col5 = length(grep(col5, cPalette$PhylumCol));num_col5
  num_col6 = length(grep(col6, cPalette$PhylumCol));num_col6
  num_col7 = length(grep(col7, cPalette$PhylumCol));num_col7

  
  # Synergistetes
  number_of_other_phyla = num_of_final_phyla - ((num_col1 > 0) + (num_col2 > 0) + (num_col3 > 0) +(num_col4 > 0)+  (num_col5 > 0) + (num_col6 > 0)+ (num_col7 > 0))
  
  #print(number_of_other_phyla)
  # col1 = green_pallete
  cPalette[grep(col1, cPalette$PhylumCol), -1] = expand.grid(h=0.4, s=seq(0.3,1,length.out=num_col1), v=0.9)
  
  # col2 = purple_pallete
  cPalette[grep(col2, cPalette$PhylumCol), -1] = expand.grid(h=0.8, s=seq(0.3,1,length.out=num_col2), v=0.9) 

  # col3 = blue_pallete
  cPalette[grep(col3, cPalette$PhylumCol), -1] = expand.grid(h=0.6, s=seq(0.3,1,length.out=num_col3), v=0.9)

  # col4 = orange_pallete
  cPalette[grep(col4, cPalette$PhylumCol), -1] = expand.grid(h=0.2, s=seq(0.3,1,length.out=num_col4), v=0.9)
  
  # col5 = red_pallete
  cPalette[grep(col5, cPalette$PhylumCol), -1] = expand.grid(h=0, s=seq(0.3,1,length.out=num_col5), v=0.9)
  
  # col6 = brown_pallete
  cPalette[grep(col6, cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_col6), v=1)
  
  # col7 = yellow_pallete
  cPalette[grep(col7, cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_col7), v=1)
  
  #print(cPalette)
  #print(number_of_other_phyla)
  
  
  print(cPalette)
  
  # Add other and species name
  cPalette$PhylumCol <- taxaName
  other<-data.frame("[1_#Other]",0,0,0.75)
  names(other)<-c("PhylumCol", "h","s","v")
  color_table <- rbind(other, cPalette)
  
  
  ## hsv to color code ##
  ## taxa위치 변경을 용이 하게 하려면 hsv to color code를 해야 한다. 
  taxa_vs_color = cbind(color_table, apply(color_table[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)];taxa_vs_color
  colnames(taxa_vs_color) <- c("Taxa", "Color");taxa_vs_color
  class(taxa_vs_color)
  coloring <- as.character(taxa_vs_color$Color) 
  names(coloring) <- taxa_vs_color$Taxa;coloring
  
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$color_table <- color_table
    results$coloring <-coloring
    return(results) 
  }
  cat("\n")
  print("$color_table and $coloring are returned.")
  
  functionReturningTwoValues()
}
removeRows <- function(rowNum, data) {
    newData <- data[-rowNum, , drop = FALSE]
    rownames(newData) <- NULL
    newData
}




getTaxonomy <- function(otus, tax_tab, level, taxRanks,na_str = c( "NA")) {
    ranks <- taxRanks
    sel <- ranks[1:match(level, ranks)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
    retval[inds!=match(level, ranks)] <- paste(na_str[1], retval[inds!=match(level, ranks)], sep="_")
    return(retval)
}


getTaxonomyAll <- function(otus, tax_tab, level, na_str = c( "NA")) {
    rank <- rank
    sel <- rank[1:match(level, rank)]
    inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str))))
    retval <- as.data.frame(tax_tab)[cbind(otus, rank[inds])]
    retval[inds!=match(level, rank)] <- paste(na_str[1], retval[inds!=match(level, rank)], sep="_")
    return(retval)
}



# X is indicator matrix of predictions, Y is indicator matrix of truth
# columns are classes, rows are samples
mcc <- function(preds=NULL, actuals=NULL, x=NULL, y=NULL) {
    # if preds and actuals are provided, x and y will be ignored
    if (!is.null(preds)) {
        nclasses <- length(union(preds, actuals))
        x <- matrix(0, nrow=length(preds), ncol=nclasses)
        y <- matrix(0, nrow=length(actuals), ncol=nclasses)
        x[cbind(1:nrow(x), preds+1)] <- 1
        y[cbind(1:nrow(y), actuals+1)] <- 1
    }
    if (!all(dim(x) == dim(y))) {
        stop("X and Y must have the same dimensions")
    }
    
    cov_biased <- function(x, y) {
        sum(sapply(1:ncol(x), function(k) {
            cov(x[,k], y[,k]) # unbiased estimate with (n-1) denominator as opposed to (n), but cancels out anyways so identical result
        }))
    }
    numerator <- cov_biased(x,y)
    denominator <- sqrt(cov_biased(x,x) * cov_biased(y,y))
    numerator / denominator
}


multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    
    i = 1
    while (i < numPlots) {
        numToPlot <- min(numPlots-i+1, cols*rows)
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
        if (numToPlot==1) {
            print(plots[[i]])
        } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            # Make each plot, in the correct location
            for (j in i:(i+numToPlot-1)) {
                # Get the i,j matrix positions of the regions that contain this subplot
                matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
                print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col))
            }
        }
        i <- i+numToPlot
    }
}

normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
    tab <- table(fvec)
    retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
    #    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
    newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
    retval <- factor(retval, levels=newlevels)
    if (originalLevelsAsNames) {
        names(retval) <- fvec
    }
    return(retval)
}






# ============================================================
# Tutorial on plotting significant taxa and environmental variables on an NMDS plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(vegan)
library(ggplot2)
library(grid)


bv.step <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
max.rho=0.95,
min.delta.rho=0.001,
random.selection=TRUE,
prop.selected.var=0.2,
num.restarts=10,
var.always.include=NULL,
var.exclude=NULL,
output.best=10
){
    
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
    require(vegan)
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    
    fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
    
    #an initial removal phase
    var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
    full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
    var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
    RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
    for(i in 1:dim(var.comb)[2]){
        var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
        RES$rho[i] <- temp$estimate
    }
    delta.rho <- RES$rho - full.cor
    exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
    
    if(random.selection){
        num.restarts=num.restarts
        prop.selected.var=prop.selected.var
        prob<-rep(1,ncol(var.mat))
        if(prop.selected.var< 1){
            prob[exclude]<-0
        }
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    } else {
        num.restarts=1
        prop.selected.var=1
        prob<-rep(1,ncol(var.mat))
        n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
    }
    
    RES_TOT <- c()
    for(i in 1:num.restarts){
        step=1
        RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
        attr(RES$step.dir, "levels") <- c("F","B")
        best.comb <- which.max(RES$rho)
        best.rho <- RES$rho[best.comb]
        delta.rho <- Inf
        selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
        while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
            #forward step
            step.dir="F"
            step=step+1
            var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
            if(RES$n.var[best.comb] == 0){
                var.comb.incl<-1:length(var.comb)
            } else {
                var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                temp <- NA*1:length(var.comb)
                for(j in 1:length(temp)){
                    temp[j] <- all(var.keep %in% var.comb[[j]])
                }
                var.comb.incl <- which(temp==1)
            }
            
            RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
            for(f in 1:length(var.comb.incl)){
                var.incl <- var.comb[[var.comb.incl[f]]]
                var.incl <- var.incl[order(var.incl)]
                var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                RES.f$var.incl[f] <- paste(var.incl, collapse=",")
                RES.f$rho[f] <- temp$estimate
            }
            
            last.F <- max(which(RES$step.dir=="F"))
            RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
            best.comb <- which.max(RES$rho)
            delta.rho <- RES$rho[best.comb] - best.rho
            best.rho <- RES$rho[best.comb]
            
            if(best.comb == step){
                while(best.comb == step & RES$n.var[best.comb] > 1){
                    #backward step
                    step.dir="B"
                    step <- step+1
                    var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
                    var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
                    RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
                    for(b in 1:length(var.comb)){
                        var.incl <- var.comb[[b]]
                        var.incl <- var.incl[order(var.incl)]
                        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
                        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
                        RES.b$var.incl[b] <- paste(var.incl, collapse=",")
                        RES.b$rho[b] <- temp$estimate
                    }
                    RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
                    best.comb <- which.max(RES$rho)
                    best.rho<- RES$rho[best.comb]
                }
            } else {
                break()
            }
            
        }
        
        RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
        print(paste(round((i/num.restarts)*100,3), "% finished"))
    }
    
    RES_TOT <- unique(RES_TOT[,3:5])
    
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
    } else {
        order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
    }
    rownames(order.by.best)<-NULL
    
    order.by.i.comb <- c()
    for(i in 1:length(selected.var)){
        f1 <- which(RES_TOT$n.var==i)
        f2 <- which.max(RES_TOT$rho[f1])
        order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
    }
    rownames(order.by.i.comb)<-NULL
    
    if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
    out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
    )
    out
    
}

bio.env <- function(fix.mat, var.mat,
fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
scale.fix=FALSE, scale.var=TRUE,
output.best=10,
var.max=ncol(var.mat)
){
    if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
    if(var.max > dim(var.mat)[2]){stop("var.max cannot be larger than the number of variables (columns) in var.mat")}
    
    require(vegan)
    
    combn.sum <- sum(factorial(ncol(var.mat))/(factorial(1:var.max)*factorial(ncol(var.mat)-1:var.max)))
    
    if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
    if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
    fix.dist <- vegdist(fix.mat, method=fix.dist.method)
    RES_TOT <- c()
    best.i.comb <- c()
    iter <- 0
    for(i in 1:var.max){
        var.comb <- combn(1:ncol(var.mat), i, simplify=FALSE)
        RES <- data.frame(var.incl=rep(NA, length(var.comb)), n.var=i, rho=0)
        for(f in 1:length(var.comb)){
            iter <- iter+1
            var.dist <- vegdist(as.matrix(var.mat[,var.comb[[f]]]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES$var.incl[f] <- paste(var.comb[[f]], collapse=",")
            RES$rho[f] <- temp$estimate
            if(iter %% 100 == 0){print(paste(round(iter/combn.sum*100, 3), "% finished"))}
        }
        
        order.rho <- order(RES$rho, decreasing=TRUE)
        best.i.comb <- c(best.i.comb, RES$var.incl[order.rho[1]])
        if(length(order.rho) > output.best){
            RES_TOT <- rbind(RES_TOT, RES[order.rho[1:output.best],])
        } else {
            RES_TOT <- rbind(RES_TOT, RES)
        }
    }
    rownames(RES_TOT)<-NULL
    
    if(dim(RES_TOT)[1] > output.best){
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)[1:output.best]
    } else {
        order.by.best <- order(RES_TOT$rho, decreasing=TRUE)
    }
    OBB <- RES_TOT[order.by.best,]
    rownames(OBB) <- NULL
    
    order.by.i.comb <- match(best.i.comb, RES_TOT$var.incl)
    OBC <- RES_TOT[order.by.i.comb,]
    rownames(OBC) <- NULL
    
    out <- list(
    order.by.best=OBB,
    order.by.i.comb=OBC,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(OBB$var.incl[1], ",")))], collapse=",") ,
    best.model.rho=OBB$rho[1]
    )
    out
}



heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
distfun = dist, hclustfun = hclust, dendrogram = c("both",
"row", "column", "none"), reorderfun = function(d, w) reorder(d,
w), symm = FALSE, scale = c("none", "row", "column"),
na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
col = "heat.colors", colsep, rowsep, sepcolor = "white",
sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
na.color = par("bg"), trace = c("column", "row", "both",
"none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks),
linecol = tracecol, margins = c(5,5,5,5), ColSideColors, RowSideColors, side.height.fraction=0.3,
cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,
NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
key = TRUE, keysize = 1.5, density.info = c("histogram",
"density", "none"), denscol = tracecol, symkey = any(x <
0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL,
key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL,
key.par = list(), main = NULL, xlab = NULL, ylab = NULL,
lmat = NULL, lhei = NULL, lwid = NULL, ColSideColorsSize = 1, RowSideColorsSize = 1, extrafun = NULL, ...)
{
    library(gtools)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
    "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
    "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
    else if (all(Colv == "Rowv"))
    Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 4)
    stop("`margins' must be a numeric vector of length 4")
    if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
            dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
            dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
        nr))
        stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        browser()
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd >
        nc))
        stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
    (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
    (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
    1) {
        if (missing(col) || is.function(col))
        breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
    col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
            stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
            stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[4]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            #            par(mar = c(0.5, 0, 0, margins[4]))
            par(mar = c(0.5, margins[2], 0, margins[4]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
    par(mar = margins)
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
    retval$rowDendrogram <- ddr
    if (exists("ddc"))
    retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
    }
    if (is.null(srtCol))
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
    offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
    padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol))
            adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
            tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
            strheight("M"), labels = labCol, adj = adjCol,
            cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        par(mar = c(margins[1L], 0, 0, margins[4L]))
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
        tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2,
            line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
            y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
            srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[4] - 1.25)
    if (!missing(add.expr))
    eval(substitute(add.expr))
    if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
    xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
    1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
    1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
    1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
    col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol,
                lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
    col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE,
        yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[4]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i",
        leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
            stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab))
        mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab))
        mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title))
        mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0)
        do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row")
            key.xlab <- "Row Z-Score"
            else if (scale == "column")
            key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) *
                0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Density"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
            key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title))
            title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
            key.ylab <- "Count"
            if (!is.na(key.ylab))
            mtext(side = 2, key.ylab, line = par("mgp")[1],
            padj = 0.5)
        }
        else if (is.null(key.title))
        title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun))
    extrafun()
    invisible(retval)
}


normalizeByRows <- function (df, rsum=1)
{
    while (any(abs((rowSums(df)-rsum))>1e-13)) {
        df <- rsum*(df / rowSums(df))
    }
    return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
    if (is.null(level)) {
        while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
            missing <- which(colSums(df)==0)
            df <- sweep(df, 2, colSums(df)/csum, "/")
            df[,missing] <- 0
        }
    } else {
        tmp <- df
        tmp$taxa <- rownames(tmp)
        tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
        names <- rownames(tmp)[order(tmp$splitter)]
        tmp <- ddply(tmp, .(splitter), function(x) {
            x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
            while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
                x <- sweep(x, 2, colSums(x)/csum, "/")
            }
            x
        })
        rownames(tmp) <- names
        df <- tmp[, setdiff(colnames(tmp), "splitter")]
    }
    return(df)
}





data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}




# Title: Geometric Mean of Pairwise Ratios (GMPR) for Microbiome Sequencing data normalization
# Version: 0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Date: 2017/02/07
# Description: The function calculates the normalizing factors for microbiome sequencing data or, more generally, zeroinflated sequencing data.
# The size factors can be used as offsets in count-based regression models or as divisors to produce normalized data


require(matrixStats)

GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'),
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n',
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}
Go_cleanMito <- function(psIN, project){
  
  
  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  
  #====== step 1 removing mitochondria reads
  seqtab.nochim <- as.matrix(otu_table(psIN))
  tax <- as.matrix(tax_table(psIN))

  # is.mito <- tax[,"Order"] %in% "Rickettsiales" 
  is.a <- tax[,"Order"] %in% "Rickettsiales" 
  seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
  tax.a <- tax[!is.a,]
  
  is.NA <- tax.a[,"Phylum"] %in% NA 
  seqtab.noNA <- seqtab.a[,!is.NA] ;dim(seqtab.noNA)
  tax.noNA <- tax.a[!is.NA,]
  # remove na column sum
  seqtab.noNA <-data.frame(t(seqtab.noNA))
  seqtab.noNA.sum <- seqtab.noNA[,colSums(seqtab.noNA) > 1];dim(seqtab.noNA.sum)
  seqtab.Nomito <- as.matrix(t(seqtab.noNA.sum))
  
  seqs <- getSequences(seqtab.Nomito)
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  
  write(fasta, file=sprintf("%s/%s.%s.No_mito.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  #====== step 2 merge phyloseq
  ps.Nomito <- phyloseq(otu_table(seqtab.Nomito, taxa_are_rows=FALSE), tax_table(tax.noNA));ps.Nomito
  sample_names(ps.Nomito)
  sample_names(ps.Nomito) <- gsub("X","",sample_names(ps.Nomito));sample_names(ps.Nomito)

  #sample_names(ps.Nomito) <- gsub("\\_.*","",sample_names(ps.Nomito));sample_names(ps.Nomito)
  
  
  
  
  #====== step 3 get the table
  otu <- as.data.frame(t(otu_table(ps.Nomito)));dim(otu)
  tax <- tax_table(ps.Nomito);dim(tax)
  
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s.%s.asvTable_No_mito.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  saveRDS(ps.Nomito, sprintf("%s/ps_No_mito.%s.%s.rds", rds, project,format(Sys.Date(), "%y%m%d")))
  
  
  print(psIN)
  print(ps.Nomito)  
}





Go_ancombc <- function(psIN,project, metaData, adjust,taxanames,filter,name){
  # outpur files
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  out_ancombs <- file.path(sprintf("%s_%s/table/ancombs",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_ancombs)) dir.create(out_ancombs)
  
  out_ancombs.Tab <- file.path(sprintf("%s_%s/table/ancombs/tab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_ancombs.Tab)) dir.create(out_ancombs.Tab)
  
  out_ancombs.ps <- file.path(sprintf("%s_%s/table/ancombs/ps",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_ancombs.ps)) dir.create(out_ancombs.ps)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  
  # taxa aggregation
  if(!is.null(taxanames)){
    psIN <- aggregate_taxa(psIN, taxanames)
  }else{
    psIN <- psIN
  }
  
  # map 정리
  mapping <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)
  
  
  
  
  
  for (mvar in rownames(subset(metadata.sel, Go_ancombc =="yes"))) {
    
    # NA remove
    mapping.sel[mapping.sel==""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    na.count <- length(mapping.sel.na)
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
    
    
    
    if (length(unique(mapping.sel[, mvar])) == 1) {
      next
    }
    # integer control
    if (class(mapping.sel.na.rem[,mvar]) == "character"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    
    # combination
    mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = orders)
    mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
    my_comparisons <- {}
    for(i in 1:ncol(cbn)){
      x <- cbn[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    # subset sample by combination
    for(i in 1:length(my_comparisons)){
      print(my_comparisons[i])
      combination <- unlist(my_comparisons[i]);combination
      basline <-combination[1]
      smvar <- combination[2]
      
      mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar)) 
      
      mapping.sel.cb[,mvar] <- factor(mapping.sel.cb[,mvar])
      psIN.cb <- psIN.na
      
      sample_data(psIN.cb) <- mapping.sel.cb;dim(mapping.sel.cb)
      
      psIN.cb <- Go_filter(psIN.cb, cutoff = filter) #0.00005
      
      unique(mapping.sel.cb[,mvar])
      
      summary(mapping.sel.cb[,mvar])
      
      
      
      if(!is.null(adjust)){
        out <- ancombc(phyloseq = psIN.cb, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                       formula = sprintf("%s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")), 
                       group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                       max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
      }else{
        #out <- ancombc(phyloseq = ps.taxa.rel, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
        #               formula = mvar, 
        #               group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
        #               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
        
        out <- ancombc(phyloseq = psIN.cb, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                       formula = mvar, 
                       group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                       max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
      }
      
      # data table with 
      res = out$res
      res.df <- as.data.frame(res)
      colnames(res.df) <- gsub(mvar,"", colnames(res.df))
      colnames(res.df)<-c("best","se","W","pval","qval","diff_abn")
      
      res.df$basline <- basline
      res.df$smvar <- smvar
      
      res.df <- cbind(as(res.df, "data.frame"), as(tax_table(psIN.cb)[rownames(res.df), ], "matrix"))
      
      res.or_p <- res.df[order(res.df$pval),]
      
      res.or_p.sel <- subset(res.or_p, res.or_p$diff_abn  == T);dim(res.or_p.sel)[1]
      taxa_sig <- rownames(res.or_p)[1:dim(res.or_p.sel)[1]]; summary(taxa_sig)
      
      
      # selected ps
      ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
      print(ps.taxa.sig)
      
      
      
      # save table and phyloseq object
      
      
      if(!is.null(taxanames)){
        if (!is.null(name)){
          write.csv(res.df, quote = FALSE,col.names = NA,#row.names = FALSE, 
                    file=sprintf("%s/ancombdTab.%s.%s.(%svs%s).T%s.%s.%s.%s.csv",out_ancombs.Tab,project,mvar,basline,smvar,dim(res.or_p.sel)[1],name,taxanames,format(Sys.Date(), "%y%m%d"), sep="/"))
          saveRDS(ps.taxa.sig, sprintf("%s/ps.ancom.sigTaxa.%s.%s.(%svs%s),T%s.%s.%s.%s.rds", out_ancombs.ps, project,mvar,basline,smvar,dim(res.or_p.sel)[1],name, taxanames,format(Sys.Date(), "%y%m%d"), sep="/"))
        }else{
          write.csv(res.df, quote = FALSE,col.names = NA,#row.names = FALSE, 
                    file=sprintf("%s/ancombdTab.%s.%s.(%svs%s).T%s.%s.%s.csv",out_ancombs.Tab,project,mvar,basline,smvar,dim(res.or_p.sel)[1],taxanames,format(Sys.Date(), "%y%m%d"), sep="/"))
          saveRDS(ps.taxa.sig, sprintf("%s/ps.ancom.sigTaxa.%s.%s.(%svs%s).T%s.%s.%s.rds", out_ancombs.ps, project,mvar,basline,smvar,dim(res.or_p.sel)[1], taxanames,format(Sys.Date(), "%y%m%d"), sep="/"))
        }
      }else{
        if (!is.null(name)){
          write.csv(res.df, quote = FALSE,col.names = NA,#row.names = FALSE, 
                    file=sprintf("%s/ancombdTab.%s.%s.(%svs%s).T%s.%s.%s.csv",out_ancombs.Tab,project,mvar,basline,smvar,dim(res.or_p.sel)[1],name,format(Sys.Date(), "%y%m%d"), sep="/"))
          saveRDS(ps.taxa.sig, sprintf("%s/ps.ancom.sigTaxa.%s.%s.(%svs%s),T%s.%s.%s.rds", out_ancombs.ps, project,mvar,basline,smvar,dim(res.or_p.sel)[1],name, format(Sys.Date(), "%y%m%d"), sep="/"))
        }else{
          write.csv(res.df, quote = FALSE,col.names = NA,#row.names = FALSE, 
                    file=sprintf("%s/ancombdTab.%s.%s.(%svs%s).T%s.%s.csv",out_ancombs.Tab,project,mvar,basline,smvar,dim(res.or_p.sel)[1],format(Sys.Date(), "%y%m%d"), sep="/"))
          saveRDS(ps.taxa.sig, sprintf("%s/ps.ancom.sigTaxa.%s.%s.(%svs%s).T%s.%s.rds", out_ancombs.ps, project,mvar,basline,smvar,dim(res.or_p.sel)[1], format(Sys.Date(), "%y%m%d"), sep="/"))
        }
      }
    }
  }
}


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

#' A Go_DA
#'



Go_DA <- function(psIN,  project, filter, taxanames=NULL, data_type = "other", 
                  cate.vars,  cate.conf=NULL, orders=NULL,
                  name=NULL, fdr=0.05){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/table/Differential_Abundance",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA)) dir.create(out_DA)
  
  out_DA.Tab <- file.path(sprintf("%s_%s/table/Differential_Abundance/tab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA.Tab)) dir.create(out_DA.Tab)
  
  out_DA.ps <- file.path(sprintf("%s_%s/table/Differential_Abundance/ps",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA.ps)) dir.create(out_DA.ps)
  
  
  # taxa aggregate
  if(!is.null(taxanames)){
    psIN <- aggregate_taxa(psIN, taxanames)
  }else{
    psIN <- psIN
  }
  mapping <- data.frame(sample_data(psIN))
  
  # start
  res <- {}
  for (mvar in cate.vars) {
    if (length(unique(mapping[, mvar])) == 1) {
      next
    }

    #na remove
    mapping.sel <- data.frame(sample_data(psIN))
    mapping.sel[mapping.sel==""] <- "NA"
    mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
    na.count <- length(mapping.sel.na)
    psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
    mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))

    if (length(unique(mapping.sel.na.rem[,mvar])) == 1 )
      next

   print(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    if (length(mapping.sel.na.rem[,mvar]) < 4){
      next
      print(sprintf("%s is removed because length(%s) less than 4", mvar, length(mapping.sel.na.rem[,mvar])))
    }


    
    
    # integer control
    if (class(mapping.sel.na.rem[,mvar]) == "character"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }
    if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
      mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
      sample_data(psIN.na) <- mapping.sel.na.rem
    }



    # for changing "-" to character 
    mapping.sel[,mvar] <- gsub("V-","Vn",mapping.sel[,mvar])

    # combination
    if(!is.null(orders)){
     mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = intersect(orders, mapping.sel[,mvar]))
    }else{
     mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    }
    
    # mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
    cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
    
    my_comparisons <- {}
    for(i in 1:ncol(cbn)){
      x <- cbn[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    # subset sample by combination
    for(i in 1:length(my_comparisons)){
    print(my_comparisons[i])
    combination <- unlist(my_comparisons[i]);combination
    basline <- combination[1]
    smvar <- combination[2]
    
    mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar)) # phyloseq subset은 작동을 안한다.
    
    psIN.cb <- psIN.na
    sample_data(psIN.cb) <- mapping.sel.cb
    
    psIN.cb <- Go_filter(psIN.cb, cutoff = filter); #0.00005

    #-- DESeq2 for phyloseq --#
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    if (length(cate.conf) >= 1) {
      form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")))
      print(form)
      dds = phyloseq_to_deseq2(psIN.cb, form)
    }    else {
      dds = phyloseq_to_deseq2(psIN.cb, as.formula(sprintf("~ %s", mvar)))
      print(sprintf("~ %s", mvar))
    }

    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
    dds = estimateDispersions(dds)
    vst = getVarianceStabilizedData(dds)
    dds = DESeq(dds, fitType="local")
    resultsNames(dds)
    
   
    
    #-- ANCOM-bc for phyloseq --#
    if(!is.null(cate.conf)){
      out <- ancombc(phyloseq = psIN.cb, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                     formula = sprintf("%s + %s", mvar, paste(setdiff(cate.conf, "SampleType"), collapse="+")), 
                     group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
    }else{
      #out <- ancombc(phyloseq = ps.taxa.rel, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
      #               formula = mvar, 
      #               group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
      #               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
      
      out <- ancombc(phyloseq = psIN.cb, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                     formula = mvar, 
                     group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
    }
    
    
    res.ancom = out$res
    
    res.ancom.df <- as.data.frame(res.ancom)
    colnames(res.ancom.df) <- gsub(mvar,"", colnames(res.ancom.df))

    
    if(!is.null(cate.conf)){
      names(res.ancom.df)[length(names(res.ancom.df))]<-"diff_abn" 
    }else{
      colnames(res.ancom.df)<-c("best","se","W","pval","qval","diff_abn")
    }
    

    # calculation
      print("pass2")
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))
      
      tmp$deseq2 <- ifelse(tmp$padj < fdr, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      # merge deseq2 + amcom 
      tmp$ancom <- factor(res.ancom.df$diff_abn[match(rownames(tmp), rownames(res.ancom.df))]);head(tmp$ancom)
      



      tmp$mvar <- mvar
      tmp$basline<-basline
      tmp$bas.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == basline))
      tmp$smvar <- smvar
      tmp$smvar.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == smvar))

      
      
      

      
      #-- give taxa name --#
      res <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
      print("pass3")
        taxaRanks <- c("Kingdon","Phylum","Class","Order","Family","Genus","Species")
        for(t in 2:length(taxaRanks)){
          
          if (!is.null(taxanames)) {
            if (taxanames == taxaRanks[t-1]){
              break
            }
          }

          res[,taxaRanks[t]] == "NA"
          res[,taxaRanks[t]]<- as.character(res[,taxaRanks[t]])
          res[,taxaRanks[t]][is.na(res[,taxaRanks[t]])] <- "__"
          
          for(i in 1:length(res[,taxaRanks[t]])){
            if (res[,taxaRanks[t]][i] == "s__" || res[,taxaRanks[t]][i] == "g__" || res[,taxaRanks[t]][i] == "f__" || res[,taxaRanks[t]][i] == "o__" || res[,taxaRanks[t]][i] == "c__"|| res[,taxaRanks[t]][i] == "p__"|| res[,taxaRanks[t]][i] == "__"){
              res[,taxaRanks[t]][i] <- ""
            }
          } 
        }
        
        
        print("pass4")
        res$TaxaName <- paste(res$Phylum,"",res$Class,"",res$Order,"",res$Family,"",res$Genus,"",res$Species)
        
        #res$ShortName <- paste(res$Phylum,res$Family," ",res$Genus," ",res$Species)
        
        
        res$Species[res$Species=="NA NA"] <- "  "
        
        if (!is.null(taxanames)) {
          if (data_type == "dada2" | data_type == "DADA2") {
            if(taxanames == "Species"){
              res$ShortName <- paste(res$Genus,"",res$Species)
            }else{
              res$ShortName <- paste(res[,taxanames],"",res$Species)
            }
            
          }
          else if (data_type == "Nephele" | data_type == "nephele") {
            res$ShortName <- paste(res[,taxanames],"",res$Species)
          }
          else if (data_type == "other" | data_type == "Other") {
            res$ShortName <- paste(res[,taxanames])
          }
          
        }else{
          if (data_type == "dada2" | data_type == "DADA2") {
            res$ShortName <- paste(res$Genus,"",res$Species)
          }
          else if (data_type == "Nephele" | data_type == "nephele") {
            res$ShortName <- paste(res$Genus,"",res$Species)
          }
          else if (data_type == "other" | data_type == "Other") {
            res$ShortName <- paste(res$Species)
          }
        }
        

        # use last taxa name
        for(taxa in c("Family", "Order", "Class","Phylum")){
          for(i in 1:length(res[,taxa])){
            if (res$ShortName[i] != "  "){
              next
            }      else if (res$ShortName[i] == "  " & res[,taxa][i] != ""){
              res$ShortName[i] <- paste(res[,taxa][i])
            }
          }
        }
      
      #--- give simple name to res---#
      #headers <- vector(dim(res)[2], mode="character")
      #for (i in 1:dim(res)[1]) {
      #  headers[i] <- paste("ASV", i, sep="_")
      #}
      headers <- rownames(res)
      
      
      res$taxa <- headers
      print("pass5")
      #-- create table --#
      res <- as.data.frame(res)
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$dir <- ifelse(res$padj < fdr, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
      
      
      # get ps objectonly significant taxa 
      res.sel <- subset(res, res$ancom  == T & !(res$deseq2  == "NS"));dim(res.sel)[1]
      taxa_sig <- rownames(res.sel)[1:dim(res.sel)[1]]; summary(taxa_sig)
      
      if(dim(res.sel)[1] == 0){
        ps.taxa.sig <- psIN.cb
      }else{
        ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
        print(ps.taxa.sig)
      }
      
      # "name definition
      if (class(name) == "function"){
        name <- NULL
      }

      # for changing "n" to "-" 
      res$basline <- gsub("Vn","V-",res$basline)
      res$smvar <- gsub("Vn","V-",res$smvar)

      res <- arrange(res, res$padj)

      write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).Sig%s.%s.%s%s%s%s.DA.csv",out_DA.Tab,
                                                               basline, 
                                                               smvar,
                                                               dim(res.sel)[1],
                                                               mvar,
                                                               ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                                               ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")), 
                                                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                               project, sep="/"))
      
      saveRDS(ps.taxa.sig,sprintf("%s/(%s.vs.%s).Sig%s.%s.%s%s%s%s.ancom.rds",out_DA.ps,
                                  basline, 
                                  smvar,
                                  dim(res.sel)[1],
                                  mvar, 
                                  ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                  ifelse(is.null(cate.conf), "", paste("with_confounder", ".", sep = "")), 
                                  ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                  project))
      
      
    }
  }
}

#' A Go_deseq2_volc
#'


Go_DA_plot <- function(project, file_path,files, type="taxonomy", plot = "volcano", fdr, fc, mycols=NULL, name, overlaps=10, font, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/pdf/DA_plot",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA)) dir.create(out_DA)
  
  # add input files
  path <- file_path
  filenames <- list.files(path, pattern=files);filenames

  print(path)
  print(filenames)
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }


  
  for (fn in 1:length(filenames)) {

    df <- read.csv(sprintf("%s/%s",path, filenames[fn]), row.names=NULL ,check.names=FALSE,quote = "")
    # remove NA
    df[df==""] <- "NA"
    df$deseq2 <- ifelse(df$padj < fdr & abs(df$log2FoldChange) > fc, ifelse(sign(df$log2FoldChange)==1, "up", "down"), "NS")
    df.na <- df[!is.na(df$dir), ]

    basline <- unique(df$basline)
    smvar <- unique(df$smvar)
    mvar <- unique(df$mvar)


    # colors and names
    df.na$deseq2<- gsub('down',basline, gsub('up',smvar, df.na$deseq2))
    
    df.na$deseq2 <- factor(df.na$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
    

   if(!is.null(mycols)){
    dircolors <- c(mycols[1], "grey",mycols[2]); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }else{
     dircolors <- c("#f8766d", "grey","#7cae00"); names(dircolors) <- c(as.character(basline), "NS", as.character(smvar))
    }
    
    legend.labs <- 
      c(paste(basline, " (n=", unique(df.na$bas.count),")",sep=""),
        paste("NS"),
        paste(smvar, " (n=", unique(df.na$smvar.count), ")",sep=""))

    

    df.na$ancom[is.na(df.na$ancom)] <- FALSE
    ancomshape <- c(18,5); names(ancomshape) <- c(TRUE, FALSE)# 16,1,1
    
    #------------------#
    #   plot style     #
    #------------------#
    if(plot == "volcano"){
      print("Generating Volcano plots.")
      p1 <- ggplot(data=df.na, aes(x=log2FoldChange, y=-log10(pvalue),colour=deseq2)) + 
        xlab("log2 fold change") + ylab("-log10 (p-value)")+ 
        geom_vline(xintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1) 
      
    } else if(plot == "maplot"){
      print("Generating M (log ratio) A (mean average)  plots.")
      p1 <-  ggplot(df.na, aes(x=log2(baseMean +1), y=log2FoldChange, colour=deseq2)) +
        xlab("Log2 mean of normalized counts") + ylab("Log2 fold change")+ 
        geom_hline(yintercept = c(-log2(fc), 0,log2(fc)),col = dircolors, linetype = "dotted", size = 1)
      
    } else if(plot == "forest"){
      print("Generating forest plots.")
      resSig <- as.data.frame(subset(df.na, padj < fdr)); resSig <- resSig[order(resSig$log2FoldChange),]
      resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
      if (dim(resSig)[1] == 0 | dim(resSig.top)[1] == 0 ){
        next
      }
      
      resSig$smvar <- factor(resSig$smvar)
      lims <- max(abs(resSig$log2FoldChange) + abs(resSig$lfcSE))*1.0
      resSig.top$deseq2<- gsub('down',basline, gsub('up',smvar, resSig.top$deseq2))
      resSig.top$deseq2 <- factor(resSig.top$deseq2, levels = c(as.character(basline), "NS", as.character(smvar)))
      resSig.top$ancom[is.na(resSig.top$ancom)] <- FALSE

      p1 <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=log2FoldChange, color=deseq2)) + 
        geom_hline(yintercept=0) + geom_point(aes(shape=ancom)) + coord_flip() + theme_classic() + 
        scale_color_manual(values=dircolors, labels=legend.labs) + scale_shape_manual(values = ancomshape) + #guides(shape = "none") +
        #theme_classic() +theme_bw() 
        geom_errorbar(aes(x=taxa, ymin=log2FoldChange-lfcSE, max=log2FoldChange+lfcSE), width=0.2)  + 
        ylim(c(-lims, lims))+ xlab("Taxa") + ylab("log2FoldChange")+
        theme(text = element_text(size=font), plot.title = element_text(size=font, hjust = 1),
              axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=font,face = "italic")) #hjust =1
      
      if(type == "taxonomy" | type == "taxanomy"){
        p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s__%s__%s",as.character(resSig$Phylum),as.character(resSig$Family), as.character(resSig$ShortName))) 
      } else if(type == "bacmet" ){
        p1 <- p1 + scale_x_discrete(breaks = as.character(resSig$taxa), labels = sprintf("%s",as.character(resSig$ShortName))) 
      }
    }
    




    
    if(plot == "volcano" |  plot == "maplot"){
      p1 <- p1 + theme_bw() + scale_color_manual(values=dircolors, labels=legend.labs) + theme(text = element_text(size=font+8),plot.title = element_text(size=font+8), legend.text=element_text(size=font+8),  legend.position="bottom",legend.justification = "left",legend.box = "vertical")  +
        geom_point(aes(shape=ancom), size=font-1.5) + scale_shape_manual(values = ancomshape)
      if(type == "taxonomy" | type == "taxanomy" |type == "bacmet" ){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(ShortName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(ShortName),'')), size=font, segment.fdr = 0.25, fontface="italic",max.overlaps = overlaps )
      }else if(type == "function"){
        p1 <- p1 +  geom_text_repel(aes(label=ifelse(KOName != "NA" & df.na$padj < fdr & abs(df.na$log2FoldChange) > fc, as.character(KOName),'')), size=font,max.overlaps = overlaps)
      }
    }

    
    if(!is.null(df.na$des)){
      des <- unique(df.na$des)
      p1 <- p1 + ggtitle(sprintf("%s-%s, (padj < %s,cutoff=%s) ", mvar, des, fdr, fc))
    }else{
      p1 <- p1 + ggtitle(sprintf("%s, (padj < %s,cutoff=%s) ", mvar,  fdr, fc))
    }



      pdf(sprintf("%s/%s%s.(%s.vs.%s).%s.%s(%s.%s).%s.pdf", out_DA, 
              ifelse(is.null(plot), "", paste(plot, ".", sep = "")), 
              mvar,
              basline, 
              smvar,
              project, 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr, 
              fc, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
    print(p1)
     dev.off()
  } 
  print("plot for volcano, maplot and forest")
}#' A Go_deseq2_heat
#' 
Go_DA_heat <- function(df, project, data_type, facet,groupby,font, 
                       addnumber=TRUE,
                       fdr,fc, orders, name, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  pdf(sprintf("%s/DA.heatmap.%s.%s%s(%s.%s).%s.pdf", out_path, 
              project, 
              ifelse(is.null(facet), "", paste(facet, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              fdr, 
              fc, 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  
  resSig <- as.data.frame(subset(df, padj < fdr)); resSig <- resSig[order(resSig$log2FoldChange),]
  
  
  if (length(subset(resSig, ancom == TRUE)) > 1){
    print(sprintf("Combination Deseq2(%s) and Ancom(%s)",length(resSig$deseq2),length(subset(resSig, ancom == TRUE))))
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
    resSig.top <- subset(resSig.top, ancom == TRUE)
    
  } else{
    print(sprintf("Use only Deseq2(%s)",length(resSig.top$deseq2)))
    resSig.top <- as.data.frame(subset(resSig, abs(resSig$log2FoldChange) > fc))
  }
  
  #print("c")
  #if (length(unique(resSig$smvar)) >=2 ){
  
  if (dim(resSig)[1] >= 1) {
    # re-order
    if (!is.null(orders)) {
      resSig.top[,groupby] <- factor(resSig.top[,groupby], levels = intersect(orders, resSig.top[,groupby]))
      resSig.top[,facet] <- factor(resSig.top[,facet], levels = intersect(orders, resSig.top[,facet]))
      print("Re-ordered")
    } else {
      resSig.top[,groupby] <- factor(resSig.top[,groupby])
      resSig.top[,facet] <- factor(resSig.top[,facet])
      print(2)
    }
    


    resSig.top$basline <- paste(resSig.top$basline," (n=",resSig.top$bas.count, ")",sep="")
    resSig.top$smvar <- paste(resSig.top$smvar," (n=", resSig.top$smvar.count, ")",sep="")
    
    
    # re-order using number
    new.orders <- c()
    for(i in orders){
      if(length(order <- grep(i, unique(resSig.top$smvar)))){
        order <- c(unique(resSig.top$smvar)[order])
      }
      new.orders <- c(new.orders, order)
    }
    
    tt <- try(resSig.top$smvar  <- factor(resSig.top[,facet], levels = intersect(new.orders, resSig.top$smvar)),T)
    
    if (class(tt) =="try-error"){
      resSig.top$smvar  <- factor(resSig.top[,groupby], levels = intersect(new.orders, resSig.top$smvar))
    } else{
      resSig.top$smvar  <- factor(resSig.top[,facet], levels = intersect(new.orders, resSig.top$smvar))
    }
    


    
    

    print(1)
    if (groupby == "smvar"){
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=smvar, color=smvar)) + theme_classic()+ coord_flip() #x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
 
    }  else {
      p <- ggplot(resSig.top, aes(x=reorder(taxa,log2FoldChange), y=mvar, color=mvar)) + theme_classic()+ coord_flip()#x=reorder(taxa,Estimate); 원래 x=factor(taxa). 값에 따라 정열 하기 위해x=reorder(taxa,Estimate)를 사용함
    }

    
    p = p + geom_tile(aes(fill = log2FoldChange), colour = "white") + 
      labs(y = "Comparison Group") +labs(x = NULL) +
      scale_fill_gradient2(low = "#1170aa", mid = "white", high = "#fc7d0b")+
      ggtitle(sprintf("%s baseline %s vs %s (padj < %s, cutoff=%s) ", unique(resSig$mvar), unique(resSig$basline), "All groups",  fdr,fc))  + 
      theme(plot.title = element_text(hjust = 0.5),legend.position= "right")+ #0.5
      theme(axis.text.x = element_text(angle=0, vjust=0.5, hjust=1, size=8),
             axis.text.y = element_text(angle=0, vjust=0.5, hjust=1, size=8,face = "italic"),
            plot.title=element_text(size=9,face="bold")) 
    
    
    print(2)
    if (data_type == "dada2" | data_type == "DADA2") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$Phylum, resSig$ShortName)))
    } else if (data_type == "Other" | data_type == "other") {
      p1 = p + scale_x_discrete(breaks = as.character(resSig$taxa), labels = as.character(paste(resSig$KOName)))
    }
    
    print(3)
    if (groupby == "smvar"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"smvar"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "smvar", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  smvar, scales="free_x", ncol = 10)
      }
    }else if (groupby == "des"){
      if (length(facet) == 1) {
        ncol <- length(unique(resSig.top[,facet]))*length(unique(resSig.top[,"des"]))
        p2 = p1 + facet_wrap(as.formula(sprintf("~ %s+%s", "des", facet)), scales="free_x", ncol = ncol)
      } else {
        p2 = p1 + facet_wrap(~  des, scales="free_x", ncol = 10)
      }
    }
    #print(4)
    #plotlist[[length(plotlist)+1]] <- p
    p3 = p2 + theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + theme(text = element_text(size=font), plot.title = element_text(hjust=1))
   # print(p3)
  }else{
    next
  }

  
  p4 <- ggplotGrob(p3)
  id <- which(p4$layout$name == "title")
  p4$layout[id, c("l","r")] <- c(1, ncol(p4))
  #grid.newpage()
  grid.draw(p4)
  dev.off()
}

#' A Go_groupBox
#'

Go_groupBox <- function(psIN, mainGroup, project, orders=NULL, top=NULL, name =NULL, rank, cutoff, mycols=NULL, ylim=NULL,flip,height, width){
  
  if(!is.null(dev.list())) dev.off()
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  set.seed(151) 
  
  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }

  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }

  
  if(!is.null(top)){
    Top = names(sort(taxa_sums(psIN), TRUE)[1:top])
    ps.top = prune_taxa(Top, psIN);ps.top
  }else{
    ps.top = psIN
  }
  
  
  ### decide for  log transformation
  
  if( max(data.frame(otu_table(ps.top))) < 1){
    ps.top.rel <- ps.top
  }else{
    #ps.top.rel <- transform_sample_counts(ps.top, function(x) x / log2(x)) # log(1+x) 를 하면 NaN가 많이 나온다. 
    ps.top.rel <- transform_sample_counts(ps.top, function(x) x / sum(x)) # percent
  }

  
  tab <- data.frame(otu_table(ps.top.rel));head(tab)
  
  nsamps_threshold <- 0.01 # fraction of relabund to call a sample positive
  filt_threshold <- 0.1 # fraction of samples that need to be positive to keep an OTU for association testing
  nperm <- 100000
  
  ### aggregation by rank
  otu.filt <- as.data.frame((otu_table(ps.top.rel))) # for dada2  t(otu_table(ps.relative)
  otu.filt$func <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(ps.top.rel), taxRanks =colnames(tax_table(psIN)), level= rank)
  agg <- aggregate(. ~ func, otu.filt, sum, na.action=na.pass);dim(agg)
  
  
  
  funcNames <- agg$func
  agg <- agg[,-1]
  rownames(agg) <- funcNames
  #ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
  # agg <- agg[intersect(ftk,ftk),]
  agg$Taxa <- rownames(agg)
  # agg_t <- t(agg)
  
  
  
  aggDf<-as.data.frame(agg, row.names = agg$Taxa)
  aggDf <- aggDf[,-1]
  agg_t <- t(aggDf)
  
  colnames(agg_t)
  rownames(agg_t)
  
  # Add grouping information
  map <- sample_data(ps.top.rel);dim(map)
  df <- data.frame(agg_t, Group = map[,mainGroup]) #, name = map$StudyID,  NoOfFMT= map$NoOfFMT );head(df)
  
  df[,mainGroup] <- as.character(df[,mainGroup]);df[,mainGroup]
  df[,mainGroup][df[,mainGroup]==""] <- "NA";df[,mainGroup]
  df.na <- subset(df, df[,mainGroup] != "NA");df.na[,mainGroup]  # subset 를 사용한 NA 삭제
  df.na[,mainGroup] <- as.factor(df.na[,mainGroup]);df.na[,mainGroup]  
  
  

  group_1 <- as.factor(df.na[,mainGroup]); group_1
  
  # N <- 20
  # taxaname<-colnames(df)[1:N];taxaname
  
  
  #df.na1 <- subset(df.na, select = c(length(colnames(agg_t))) )
  #df.na <- df.na[1:N];df.na
  
  
  kruskal.wallis.table <- data.frame()
  for (i in 1:dim(df.na)[2]) {
    ks.test <- kruskal.test(df.na[,i], g=group_1)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                  data.frame(id=names(df.na)[i],
                                             p.value=ks.test$p.value
                                  ))
    # Report number of values tested
    cat(paste("Kruskal-Wallis test for ",names(df.na)[i]," ", i, "/", 
              dim(df.na)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
  }
  
  
  
  kw <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing = FALSE), ] # decreasing, increasing
  
  kw.sig <- kw[which(kw$p.value < cutoff),];dim(kw.sig)[1]
  
  cat(paste(dim(kw.sig)[1]," was p < ", cutoff,".","\n", sep=""))
  
  kw.mat <- as.matrix(kw.sig);dim(kw.mat)
  
  funcNames.sig <- kw.mat[,1];length(funcNames.sig)
  
  df.sel <- df.na[funcNames.sig]
  df.sel <- data.frame(df.sel, Group = map[,mainGroup]) 
  

  
  
  df.sel.melt <- melt(df.sel, id.vars = mainGroup, measure.vars = funcNames.sig)
  df.sel.melt$value <- as.numeric(df.sel.melt$value)
  df.sel.melt.clean <- subset(df.sel.melt, variable != "Group" &  value > 0)
  
  
  
  
  if (!is.null(orders)) {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup], levels = rev(orders))
  } else {
    df.sel.melt.clean[,mainGroup] <- factor(df.sel.melt.clean[,mainGroup])
  }
  
  df.sel.melt.clean$variable <- as.character(df.sel.melt.clean$variable)
  
  df.sel.melt.clean <- df.sel.melt.clean[order(df.sel.melt.clean$variable ,  decreasing = F), ]
  
  print(unique(df.sel.melt.clean$variable))
  p <- ggplot(df.sel.melt.clean, aes_string(x="variable", y="value", fill=mainGroup)) +  geom_boxplot(outlier.shape = NA,lwd=0.3) + 
    theme_bw() + theme(strip.background = element_blank()) + 
    labs(y="Relative abundance", x= NULL) + ggtitle(sprintf("kruskal wallis p < %s",cutoff))
  
  # + stat_compare_means(aes_string(group = mainGroup),label = "p.format") + 
    
  #+ scale_x_discrete(limits = rev)
  
  
  if(!is.null(mycols)){
    p <- p + scale_fill_manual(values = mycols)
  }else{
    p <- p
  }

  if(flip == T){
    p <- p+ coord_flip()
  }else{
    p <- p + theme(text=element_text(size=9), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  }
  

  
  # tt <- try(mycols, T)
  # if(class(tt) == "try-error"){
  #  p <- p
  # }else{
  #   p <- p + scale_fill_manual(values = mycols)
  # }

  if(!is.null(ylim)){
    p = p + ylim(ylim[1] , ylim[2])
  }else{
    p=p
  }
  
  #=== image size ===#
  #height <- 0.4*length(unique(df.sel.melt.clean[,mainGroup])) + 0.4*dim(kw.sig)[1];height
  #width <- log((max(nchar(funcNames.sig)))*max(nchar(as.character(unique(df.sel.melt.clean[,mainGroup])))));width
  print(p)
  
  pdf(sprintf("%s/groupBox.%s.%s.%s%s%s.pdf", out_path, 
              project, 
              mainGroup,
              ifelse(is.null(rank), "", paste(rank, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  print(p)
  
  dev.off()
}









Go_myCols <- function(custumCols=NULL, RColorBrewer=NULL, piratepal=NULL) {
  # reset colors
  # rm(list = ls()[grep("mycols", ls())])
  
  if(!is.null(dev.list())) dev.off()
# https://bookdown.org/ndphillips/YaRrr/more-colors.html

library("yarrr")
  # for custom
  if(is.null(custumCols) & is.null(RColorBrewer) & is.null(piratepal)){
    
    cat("#=== Please select your colors. If not, basic color would be used. ===#","\n","\n", sep=" ")
    
    cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC") # Tableau10
    cols2 <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
    
    cat("Custum Colors1: cols1","\n", cols1, "\n","\n", sep=" ")
    cat("Custum Colors2: cols2","\n", cols2, "\n","\n", sep=" ")

    cat("RColorBrewer:","\n", "Set3     Set2    Set1   Pastel2", "\n", "Pastel1  Paired  Dark2  Accent", "\n","\n", sep=" ")
    cat("piratepal:","\n", "basel    pony    Xmen    decision    southpark", "\n", "google   eternal   evildead    usualsuspects   ohbrother", 
    "\n", "appletv    brave    bugs    cars    nemo", "\n", "rat    up    espresso    ipod   info    info2","\n",sep=" ")

    display.brewer.all(type = "qual")
    readline(prompt="Press [enter] to see other colors.")
    piratepal("all")
  }
  
  if(!is.null(custumCols)){
    if(custumCols == "cols1"){
      cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC") # Tableau10
      mycols <- cols1
      
      barplot(rep(1,length(cols1)), col=cols1, main= custumCols,yaxt="n")
      return(mycols)
    }  else if(custumCols == "cols2"){
      cols2 <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
      mycols <- cols2
      barplot(rep(1,length(cols2)), col=cols2, main= custumCols,yaxt="n")
      return(mycols)
    }
  }  # for RColorBrewer
  else if(!is.null(RColorBrewer)){
    if(RColorBrewer == "Set1"){
      mycols <- brewer.pal(9, RColorBrewer)
    }else if(RColorBrewer == "Set2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Set3"){
      mycols <- brewer.pal(12, RColorBrewer)
    }else if(RColorBrewer == "Pastel2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Pastel1"){
      mycols <- brewer.pal(9, RColorBrewer)
    }else if(RColorBrewer == "Paired"){
      mycols <- brewer.pal(12, RColorBrewer)
    }else if(RColorBrewer == "Dark2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Accent"){
      mycols <- brewer.pal(8, RColorBrewer)
    }
  barplot(rep(1,length(mycols)), col=mycols, main= RColorBrewer,yaxt="n")
  cat(sprintf("mycols was set as [%s].\n.\n",RColorBrewer))
  return(mycols)
  }   # for piratepal
  else if(!is.null(piratepal)){ 
    if(piratepal == "basel"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "pony"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "Xmen"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "decision"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "southpark"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "google"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "eternal"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "evildead"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "usualsuspects"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "ohbrother"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "appletv"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "brave"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "bugs"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "cars"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "nemo"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "rat"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "up"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "espresso"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "ipod"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "info"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "info2"){
      mycols <- piratepal(palette =piratepal)
    }

 
  
  cols <- as.data.frame(mycols)
  mycols <- cols$mycols
  barplot(rep(1,length(mycols)), col=mycols, main= piratepal,yaxt="n")
  cat(sprintf("mycols was set as [%s].\n.\n",piratepal))
  return(mycols)
  }
}
#' Go_function2ps
#'
#'
#' @param function2ps
#' @keywords function2ps
#' @export
#' @examples
#' function2ps

Go_function2ps <- function(tabPath,
                           project=NULL, func.type,
                           name=NULL){
  # Read tab
  func.tab <- read_tsv(tabPath,col_types = cols())
  NumOFsample <- dim(func.tab)[2]
  
  # split a data frame
  otu <- as.matrix(func.tab[,3:NumOFsample])
  tax <- as.matrix(func.tab[,1:2])
  
  rownames(otu) <- tax[,1]
  rownames(tax) <- tax[,1]
  
  # define kegg or pathway
  if (any(grepl("K0", rownames(tax)))){
    colnames(tax) <- c("KO","KO.des")
    func <- "KEGG"
  }else if (any(grepl("PWY", rownames(tax)))){
    colnames(tax) <- c("pathway","path.des")
    func <- "pathway"
  }
  
  #merge phyloseq
  ps <- phyloseq(otu_table(otu, taxa_are_rows=T), tax_table(tax));ps
  
  print(ps)
  
  
  tt<- try(class(func.type),T)
  
  if(class(tt) == "try-error"){
    print("Please define the data type. PICRUSt or Hummann")
    break
  }else if(any(grepl(func.type, c("picrust","Picrust","Picrust2","PICRUSt","PICRUSTt2")))){
    func.type <- "PICRUSTt2"
  }else if(any(grepl(func.type, c("Human","Humann","humann","human","Human2","Humann2","humann2","human2")))){
    func.type <- "Humann2"
  }
  
  
  # saving file
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  saveRDS(ps, sprintf("%s/ps.%s.%s.%s.%s%s.rds",rds,
                      func.type,
                      func,
                      project,
                      ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                      format(Sys.Date(), "%y%m%d")))
}
