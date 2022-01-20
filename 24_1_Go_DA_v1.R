#' A Go_DA
#'



Go_DA <- function(psIN, metaData, project, order,type, filter, taxanames, data_type, adjust, des, name, alpha=0.05){

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
  
  
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  

  
  # map 정리
  mapping <- data.frame(sample_data(psIN))
  sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
  metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
  mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)

  dim(mapping.sel)
  
  
  if(!is.null(taxanames)){
    psIN <- aggregate_taxa(psIN, taxanames)
  }else{
    psIN <- psIN
  }

  
  
  # 최근 버전 for unstrafied (20210112 확인)
   if(type == "function"){
    # remove colume sum 0 and psIN 재구성(20201027)
    a <- data.frame(otu_table(psIN))*10000
    a.ceiling <- ceiling(a[-c(99),])
    b <- a.ceiling[, -which(numcolwise(sum)(a.ceiling) < 1)]
    if (length(b) == 0){
      OTU.sta <- otu_table(a, taxa_are_rows = TRUE);head(OTU.sta)
      colnames(OTU.sta) <- gsub("X", "", colnames(OTU.sta))
      otu_table(psIN) <-  OTU.sta
    }else if(length(b) > 1){
      OTU.sta <- otu_table(b, taxa_are_rows = TRUE);head(OTU.sta)
      colnames(OTU.sta) <- gsub("X", "", colnames(OTU.sta))
      otu_table(psIN) <-  OTU.sta
    }
  }else if(type == "taxanomy"){
    psIN <- psIN
  }else if(type == "bacmet"){
    psIN <- psIN
  }

  
  
  # start
  res <- {}
  for (mvar in rownames(subset(metadata.sel, Go_deseq2 =="yes"))) {
    if (length(unique(mapping.sel[, mvar])) == 1) {
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

    if (length(des) == 1) {
      print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                    des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))

    } else{
      print(sprintf("##-- %s (total without NA: %s/%s) --##",
                    mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
    }

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
    
    mapping.sel.cb <- subset(mapping.sel, mapping.sel[[mvar]] %in% c(basline, smvar)) # phyloseq subset은 작동을 안한다.
    psIN.cb <- psIN.na
    sample_data(psIN.cb) <- mapping.sel.cb
    
    psIN.cb <- Go_filter(psIN.cb, cutoff = filter); #0.00005

    #-- DESeq2 for phyloseq --#
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    if (length(adjust) >= 1) {
      form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
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
    
    
    res.ancom = out$res
    
    res.ancom.df <- as.data.frame(res.ancom)
    colnames(res.ancom.df) <- gsub(mvar,"", colnames(res.ancom.df))
    colnames(res.ancom.df)<-c("best","se","W","pval","qval","diff_abn")
    
    
    

    # calculation
      print("pass2")
      tmp <- results(dds, contrast = c(mvar, smvar, basline))
      tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
        tmp <- unlist(strsplit(x, ";"))
        tmp[length(tmp)]
      }))
      
      tmp$deseq2 <- ifelse(tmp$padj < alpha, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
      # merge deseq2 + amcom 
      tmp$ancom <- factor(res.ancom.df$diff_abn[match(rownames(tmp), rownames(res.ancom.df))]);head(tmp$ancom)
      
      tmp$mvar <- mvar
      tmp$basline<-basline
      tmp$smvar <- smvar
      if (length(des) == 1) {
        tmp$des <- des
      }
      
      
      

      
      #-- give taxa name --#
      res <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
      print("pass3")
      


      if(type == "taxonomy"){
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
        
        
        res$Species[res$Species=="NA NA"] <- ""
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
      } else if(type == "function"){
        for(taxa in c("KO", "KO.des","Path","Path.des")){
          res[,taxa] == "NA"
          res[,taxa]<- as.character(res[,taxa])
          res[,taxa][is.na(res[,taxa])] <- "__"
          for(i in 1:length(res[,taxa])){
            if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
              res[,taxa][i] <- ""
            }
          }
        }
        print("pass4")
        res$KOName <- paste(res$Path,"",res$KO)
        res$ShortName <- paste(res$Path.des,"",res$KO.des)
        
        
        unique(res$ShortName)
        
        # use last taxa name
        for(taxa in c("KO", "KO.des","Path","Path.des")){
          for(i in 1:length(res[,taxa])){
            if (res$ShortName[i] != "  "){
              next
            }      else if (res$ShortName[i] == "  " & res[,taxa][i] != ""){
              res$ShortName[i] <- paste(res[,taxa][i])
            }
          }
        }
      }else if(type == "bacmet"){
        for(taxa in c("Gene",	"Organism",	"Compound",	"NCBI_annotation")){
          res[,taxa] == "NA"
          res[,taxa]<- as.character(res[,taxa])
          res[,taxa][is.na(res[,taxa])] <- "__"
          for(i in 1:length(res[,taxa])){
            if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
              res[,taxa][i] <- ""
            }
          }
        }
        print("pass4")
        res$TaxaName <- paste(res$Compound,"",res$Gene,"",res$Organism)
        res$ShortName <- paste(res$Compound,"",res$Gene,"",res$Organism)
      }
      
      #--- give simple name to res---#
      headers <- vector(dim(res)[2], mode="character")
      for (i in 1:dim(res)[1]) {
        headers[i] <- paste("ASV", i, sep="_")
      }
      
      
      res$taxa <- headers
      print("pass5")
      #-- create table --#
      res <- as.data.frame(res)
      res$padj <- p.adjust(res$pvalue, method="fdr")
      res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
      
      
      # get ps objectonly significant taxa 
      res.sel <- subset(res, res$ancom  == T & !(res$deseq2  == "NS"));dim(res.sel)[1]
      taxa_sig <- rownames(res.sel)[1:dim(res.sel)[1]]; summary(taxa_sig)
      
      if(dim(res.sel)[1] == 0){
        ps.taxa.sig <- psIN.cb
      }else{
        ps.taxa.sig <- prune_taxa(taxa_sig, psIN.cb)
        print(ps.taxa.sig)
      }
      
      
      write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).Sig%s.%s.%s%s%s%s.DA.csv",out_DA.Tab,
                                                               basline, 
                                                               smvar,
                                                               dim(res.sel)[1],
                                                               mvar,
                                                               ifelse(is.null(des), "", paste(des, ".", sep = "")), 
                                                               ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                                               ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                                               project, sep="/"))
      
      saveRDS(ps.taxa.sig,sprintf("%s/(%s.vs.%s).Sig%s.%s.%s%s%s%s.ancom.rds",out_DA.ps,
                                  basline, 
                                  smvar,
                                  dim(res.sel)[1],
                                  mvar, 
                                  ifelse(is.null(des), "", paste(des, ".", sep = "")), 
                                  ifelse(is.null(taxanames), "", paste(taxanames, ".", sep = "")), 
                                  ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                                  project))
      
      
    }
  }
}

