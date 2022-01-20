

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
