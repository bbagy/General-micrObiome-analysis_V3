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
      agg <- normalizeByCols(agg)
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
            }
          }else{
            reg <- "GLM"
            form <- as.formula(sprintf("value ~ %s  %s", mvar, 
                                       ifelse(is.null(cate.conf), "", paste("+",setdiff(cate.conf, "SampleType"), collapse=""))));form
            print(form)
            mod <- glm(form, data=df,  family = binomial(link='logit'))
          }
            


          #exp(coef(mod))

          ## lmer에서 control=은 "number of levels of each grouping~~" 오류가 있을때만 사용한다.
          ##
          # mod2 <- lmer(form, data=df)
          
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
      # res.sel <- as.data.frame(subset(res, pvalue < pval))
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

