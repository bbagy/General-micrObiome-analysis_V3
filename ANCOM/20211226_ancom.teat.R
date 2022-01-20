##############################################
#----           ANCOM tes             ------#
##############################################
# 20211226

# https://github.com/FrederickHuangLin/ANCOM

library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(readr)
library(tidyverse)
source("~/Dropbox/04_scripts/R_source/microbiome2/ANCOM/scripts/ancom_v2.1.R")


setwd("~/Dropbox/04_scripts/R_source/microbiome2/ANCOM/")


otu_data = read_tsv("data/moving-pics-table.tsv", skip = 1)
otu_id = otu_data$`feature-id`
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data = read_tsv("data/moving-pics-sample-metadata.tsv")[-1, ]
meta_data = meta_data %>% rename(Sample.ID = SampleID)



# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "SampleID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Subject"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res$out, "outputs/res_moving_pics.csv")

# Step 3: Volcano Plot

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  



##############################################
#----           ANCOMBC test           ------#
##############################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
#https://www.yanh.org/2021/01/01/microbiome-r/#ancom-bc


#====== code start =====#

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




Go_ancombc(psIN=ps1.sel, project, metaData="3_map/metadata.CCM_phase1_ASVs.Set4.211219.csv", 
           adjust=NULL,taxanames=NULL,filter=0.0005,name=NULL)






#==========heatmap========#

Go_pheatmap <- function(psIN,project, title, group1=NULL, group2=NULL,Ntax, name=NULL,
                        show_rownames = T,show_colnames = F,
                        cluster_rows = T,cluster_cols = T, 
                        height, width){
  # BiocManager::install("ComplexHeatmap")
  # install.packages("Cairo")
  if(!is.null(dev.list())) dev.off()
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  Tableau10 = c("#B07AA1","#FF9DA7","#76B7B2", "#59A14F","#EDC948", "#9C755F", "#BABOAC","#1170aa", "#fc7d0b","#E15759") 
  
  psIN <- ps.sig
  psIN.prune = prune_samples(sample_sums(psIN) > 1, psIN);psIN.prune
  
  ps.rel <- transform_sample_counts(psIN.prune, function(x) x/sum(x)*100);ps.rel
  
  ps.rel.sel <- prune_taxa(names(sort(taxa_sums(ps.rel),TRUE)[1:Ntax]), ps.rel);ps.rel.sel
  
  
  
  matrix <- data.frame(t(otu_table(ps.rel.sel)))
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  
  
  
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
  phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Dark2")
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
    pdf(sprintf("%s/pheatmap.%s.(%s).%s.%s.pdf", out_path, project,title, name,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }   else {
    pdf(sprintf("%s/pheatmap.%s.(%s).%s.pdf", out_path, project,title, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  print(p)
  
  dev.off()
}



ps.sig <- readRDS("CCM_phase1_ASVs_211228/table/ancombs/ps/ps.ancom.sigTaxa.CCM_phase1_ASVs.HPV_CIN_3.(HPVLp_CIN1vsHPVHp_CIN2+).T21.211228.rds");ps.sig


Go_pheatmap(psIN=ps.sig,project, title="ancom", group1="hivstatus", group2="HPV_CIN_3", Ntax=21, name=NULL,
                        show_rownames = T,show_colnames = F,
                        cluster_rows = T,cluster_cols = T, 
                        height=6, width=8)
dev.off()
