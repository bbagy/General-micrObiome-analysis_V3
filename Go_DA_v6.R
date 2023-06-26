

# ancombd fixed version 20230510
Go_DA <- function(psIN,  project, filter, taxanames=NULL, data_type = "other", 
                  cate.vars,  cate.conf=NULL, cont.conf=NULL, orders=NULL,
                  name=NULL, fdr=0.05){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_DA <- file.path(sprintf("%s_%s/table/DA",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA)) dir.create(out_DA)
  
  out_DA.Tab <- file.path(sprintf("%s_%s/table/DA/tab",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_DA.Tab)) dir.create(out_DA.Tab)
  
  out_DA.ps <- file.path(sprintf("%s_%s/table/DA/ps",project, format(Sys.Date(), "%y%m%d"))) 
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

    ### categorical and continuous confounder control
    if (length(cate.conf) >= 1) {
      for(cate in cate.conf){
        mapping.sel.cb[,cate] <- as.factor(mapping.sel.cb[,cate])
        sample_data(psIN.cb) <- mapping.sel.cb
      }
    }

    if (length(cont.conf) >= 1) {
      for(cont in cont.conf){
        mapping.sel.cb[,cont] <- as.numeric(mapping.sel.cb[,cont])
        sample_data(psIN.cb) <- mapping.sel.cb
      }
    }
    
    if (!is.null(cate.conf) | !is.null(cont.conf)) {
      confounder <- c(cate.conf,cont.conf)
      
      form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(confounder, "SampleType"), collapse="+")))
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
    if (!is.null(cate.conf) | !is.null(cont.conf)) {
      confounder <- c(cate.conf,cont.conf)
      
      formula_str <- sprintf("%s + %s", mvar, paste(setdiff(confounder, "SampleType"), collapse=" + "))
      out <- ancombc2(
        data = psIN.cb,
        p_adj_method = "holm",
        lib_cut = 1000,
        fix_formula = formula_str,
        group = mvar,
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        global = TRUE,
        em_control = list(tol = 1e-5, max_iter = 100)
      )
    }else{
      #out <- ancombc(phyloseq = ps.taxa.rel, p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
      #               formula = mvar, 
      #               group = mvar, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
      #               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
      
      out <- ancombc2(
        data = psIN.cb,
        p_adj_method = "holm",
        lib_cut = 1000,
        fix_formula = mvar,
        group = mvar,
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        global = TRUE,
        em_control = list(tol = 1e-5, max_iter = 100)
      )
    }
    
    
    res.ancom = out$res
    
    
    df_mvar = res.ancom %>%
      dplyr::select(taxon, contains(mvar))
    
    
    #rownames(df_mvar) <- df_mvar$taxon 
    #names(df_mvar)[length(names(df_mvar))]<-"diff_abn" 
    names(df_mvar)<- c("taxa", "lfc_ancombc", "se_ancombc", "W_ancombc", "pvalue_ancombc", "qvalue_ancombc",  "diff_abn" )
    
    # calculation
      print("pass2")
      tmp <- as.data.frame(results(dds, contrast = c(mvar, smvar, basline)))
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))
      
      
      

      # merge deseq2 + amcom 
      #tmp$ancom <- factor(df_mvar$diff_abn[match(rownames(tmp), df_mvar$taxon)]);head(tmp$ancom)
      

      # Merge the data frames
      tmp <- merge(tmp, df_mvar, by = "taxa", all.x = TRUE)
      
      
      
      
      tmp$mvar <- mvar
      tmp$basline<-basline
      tmp$bas.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == basline))
      tmp$smvar <- smvar
      tmp$smvar.count <-  sum(with(mapping.sel.cb, mapping.sel.cb[,mvar] == smvar))

      
      #-- give taxa name --#
      res <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[tmp$taxa, ], "matrix"))
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
      res$deseq2 <- ifelse(res$padj < fdr, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
      res$ancom <- ifelse(res$diff_abn == T, ifelse(sign(res$lfc_ancombc)==1, "up", "down"), "NS")
      
      # get ps objectonly significant taxa 
      #res.sel <- subset(res, res$ancom  == T & !(res$deseq2  == "NS"));dim(res.sel)[1]
      res.sel <- subset(res, res$deseq2  %in% c("up","down"));dim(res.sel)[1]
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

