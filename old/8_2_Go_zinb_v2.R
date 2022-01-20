#' A Go_zinb
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
filt_threshold= 0.1 

Go_zinb <- function(psIN, metaData, StudyID, project, ranks, alpha, nsamps_threshold, taxRanks, data_type,name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/zinb",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  ranks <- taxRanks
  taxaname <- ranks
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

  # compute GMPR normalization
  sizeFactors <- GMPR(as.matrix(as.data.frame(otu_table(psIN))), intersect.no=9)
  # fix column types
  mapping <- data.frame(sample_data(psIN))
  for (mvar in  rownames(subset(metadata, Go_zinb=="yes" | Go_zinbConfounder=="yes"))) {
    if (metadata[mvar, "type"] == "factor") {
      mapping[,mvar] <- factor(mapping[,mvar])
      if (!(is.na(metadata[mvar, "baseline"])) && metadata[mvar, "baseline"] != "") {
        mapping[,mvar] <- relevel(mapping[,mvar], metadata[mvar, "baseline"])
      }
    } else if (metadata[mvar, "type"] == "numeric") {
      mapping[,mvar] <- as.numeric(as.character(mapping.sel[,mvar]))
    } else if (metadata[mvar, "type"] == "date") {
      mapping[,mvar] <- as.Date(sprintf("%06d", mapping.sel[,mvar]), format="%m%d%y")
      mapping[,mvar] <- factor(as.character(mapping[,mvar]), levels=as.character(unique(sort(mapping.sel[,mvar]))))
    }
  }
  

  for(i in 1:length(taxaname)){
    # dada2 or nephele
    if (data_type == "dada2" | data_type == "DADA2") {
      otu.filt <- as.data.frame(t(otu_table(psIN)))
    }
    else if (data_type == "Nephele" | data_type == "nephele") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }
    else if (data_type == "other" | data_type == "Other") {
      otu.filt <- as.data.frame(otu_table(psIN))
    }

    # continue
    otu.filt[,taxaname[i]] <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=taxRanks,level=taxaname[i])

    if (dim(otu.filt)[2] == 2){
      next
    }
    agg <- aggregate(as.formula(sprintf(". ~ %s" , taxaname[i])), otu.filt, sum, na.action=na.pass)
    genera <- agg[,taxaname[i]]
    agg <- agg[,-1]
    agg <- normalizeByCols(agg)
    rownames(agg) <- genera
    dim(agg)
    ftk <- names(which(unlist(apply(agg, 1, function(x) length(which(x>=nsamps_threshold)))) > ceiling(filt_threshold*ncol(agg))))
    agg <- agg[intersect(ftk,ftk),]
    # control data set after filter
    if (dim(agg)[1] == 0)
      next

    agg[,taxaname[i]] <- rownames(agg)

    # metatable에서 useForNB는 여러개의 yes가 가능 하지만, useAsConfounder 는 그렇지 않다.
    ## baseline등을 관리 하려면 다음이 필요하다.
    mapping <- data.frame(sample_data(psIN))
    sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
    metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
    mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)

    dim(mapping.sel)

    #-------------------------------------------------------------------#
    #--------------    ZINB+GMPR, emmeans for estimates    -------------#
    #-------------------------------------------------------------------#
    res <- {}; res.interaction <- {}
    for (f in agg[,taxaname[i]]) {
      # clean bacteria name
      if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
        next
      }
      
      df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)

      nsamps_detected <- length(which(df$value>=nsamps_threshold))
      
      for (con in rownames(subset(metadata, Go_zinb=="yes"))) {
        df[, con] <- mapping.sel[df$SampleID, con]
      }
      for (var in rownames(subset(metadata, Go_zinbConfounder=="yes"))) {
        df[, var] <- mapping.sel[df$SampleID, var]
      }
      
      df$size_factor <- log(sizeFactors[df$SampleID])
      for (mvar in rownames(subset(metadata, Go_zinb=="yes"))) {
        tt <- try(m <- zeroinfl(as.formula(sprintf("value ~ %s * %s + offset(size_factor) | 1",mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*"))), data = df, dist = "negbin", EM = F, maxit=100), silent=T) # using EM=TRUE causes certain models to hang...
        
        if (class(tt) == "zeroinfl") {
          coef <- summary(m)$coefficients$count # rows are [(Intercept), comparisons, Log(theta)], columns are [Estimate, SE, Z, pval]
          coef <- coef[grep(".", rownames(coef)),,drop=F] #:
          colnames(coef) <- c("Estimate", "SE", "t", "pval")
          res.interaction <- rbind(res.interaction, cbind(f, nsamps_detected, mvar, "ZINB", rownames(coef), coef))
          
          if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
            form <- as.formula(sprintf("pairwise ~ %s | %s", mvar, mvar))
            emm <- emmeans(m, form, adjust="none")
            coef <- as.data.frame(emm$contrasts)
            coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
            res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
          }else{
            form <- as.formula(sprintf("pairwise ~ %s | %s", mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")))
            emm <- emmeans(m, form, adjust="none")
            coef <- as.data.frame(emm$contrasts)
            coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
            #coef$confounder <- paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")
            res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
          }
          
        } else if (class(tt) == "try-error") {
          tt <- try(m <- glm.nb(as.formula(sprintf("value ~ %s * %s + offset(size_factor)",mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*"))), data = df), silent=T)
          
          if (class(tt)[1] == "negbin") {
            coef <- summary(m)$coefficients # rows are [(Intercept), comparisons], columns are [Estimate, SE, Z, pval]
            coef <- coef[grep(".", rownames(coef)),,drop=F]#:
            colnames(coef) <- c("Estimate", "SE", "t", "pval")
            res.interaction <- rbind(res.interaction, cbind(f, nsamps_detected, mvar, "NB", rownames(coef), coef))
            
            print(1)
            
            if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
              form <- as.formula(sprintf("pairwise ~ %s | %s", mvar, mvar))
              emm <- emmeans(m, form, adjust="none")
              coef <- as.data.frame(emm$contrasts)
              coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
              res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
            }else{
              form <- as.formula(sprintf("pairwise ~ %s | %s", mvar,paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")))
              emm <- emmeans(m, form, adjust="none")
              coef <- as.data.frame(emm$contrasts)
              coef$contrast <- sapply(as.character(coef$contrast), function(x) paste(rev(unlist(strsplit(x, " - "))), collapse=" - ")); coef$estimate <- -1*coef$estimate; coef$z.ratio <- -1*coef$z.ratio # reverse order of contrast
              #coef$confounder <- paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")
              res <- rbind(res, cbind(f, nsamps_detected, mvar, "NB", coef))
              }
          }
        }
      }
    }
    print(1)
    
    #form <- as.formula(sprintf("pairwise ~ %s", mvar))
    #emm <- emmeans(m, form, adjust="none")
    
    #-- tidy up results and create table --#
    if (dim(subset(metadata, Go_zinbConfounder=="yes"))[1] == 0){
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "Estimate", "SE", "df", "t", "pval")
    } else if (mvar == paste(setdiff(rownames(subset(metadata, Go_zinbConfounder=="yes")), "SampleType"), collapse="*")){
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "Estimate", "SE", "df", "t", "pval")
    }    else{
      colnames(res) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "contrast", "condition", "Estimate", "SE", "df", "t", "pval")
    }


    res$padj <- p.adjust(res$pval, method="fdr")
    res <- res[order(res$pval, decreasing=F),]
    print(2)
    
    res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
    res.interaction <- as.data.frame(res.interaction)
    colnames(res.interaction) <- c("Genus", "nsamps_detected", "metadata_variable", "model", "coefficient", "Estimate", "SE", "t", "pval"); rownames(res.interaction) <- {}
    res.interaction$Estimate <- as.numeric(as.character(res.interaction$Estimate))
    res.interaction$SE <- as.numeric(as.character(res.interaction$SE))
    res.interaction$pval <- as.numeric(as.character(res.interaction$pval))
    res.interaction$padj <- p.adjust(res.interaction$pval, method="fdr")
    res.interaction$sigstr <- ifelse(res.interaction$padj < alpha,"*", "")
  }
  
  
  
  
  if (length(name) == 1) {
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.%s.zinb.res.%s.csv",out_table, name, project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
    write.csv(res.interaction, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.%s.zinb.res.interaction.%s.csv",out_table, name, project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
  }else{
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.zinb.res.%s.csv",out_table,  project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
    write.csv(res.interaction, quote = FALSE,col.names = NA,file=sprintf("%s/%s.%s.zinb.res.interaction.%s.csv", out_table,  project,taxaname[i], format(Sys.Date(), "%y%m%d"), sep="/"))
  }
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$res <- res
    results$interaction <-res.interaction
    return(results) 
  }
  cat("\n")
  print("$res and $interaction are returned.")
  functionReturningTwoValues()
}

