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


Go_lmem <- function(psIN, metaData, StudyID, project, nsamps_threshold, filt_threshold, taxRanks, data_type, des, name){

  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/lmem",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)

  ranks <- taxRanks
  taxaname <- ranks
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

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
    print(1)

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


    print(2)
    #--------------    lmer    -------------#
    res <- {}
    for (f in agg[,taxaname[i]]) {
      # clean bacteria name
      if (f == "s__" || f == "g__" || f == "f__" || f == "o__" || f == "c__"|| f == "p__"){
        next
      }

      df <- melt(agg[f,]); colnames(df) <- c("Genus", "SampleID", "value"); df$SampleID <- as.character(df$SampleID)

      df$StudyID <- mapping.sel[df$SampleID, "StudyID"]
      
      for (cvar in rownames(subset(metadata, Go_lmemConfounder =="yes"))) {
        df[, cvar] <- mapping.sel[df$SampleID, cvar]
      }
      for (mvar in rownames(subset(metadata, Go_lmem =="yes"))) {

        # na remove
        mapping <- data.frame(sample_data(psIN))
        mapping[mapping==""] <- "NA"
        mapping.na <- mapping[!is.na(mapping[,mvar]), ]
        na.count <- length(mapping.na)
        if (length(unique(mapping.na[,mvar])) == 1)
          next

        
        #------------ fix column types------------#
        if (metadata[mvar, "type"] == "factor") {
          mapping.na[,mvar] <- factor(mapping.na[,mvar])
          if (length(unique(mapping.na[,mvar])) ==1 ){
            next
          }
          if (metadata[mvar, "baseline"] != "") {
            mapping.na[,mvar] <- relevel(mapping.na[,mvar], metadata[mvar, "baseline"])
          }
        } else if (metadata[mvar, "type"] == "numeric") {
          mapping.na[,mvar] <- factor(mapping.na[,mvar])
        }
          
        print(3)
        
        
        # na count
        if (length(des) == 1) {
          print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                        des,mvar, dim(mapping.na)[1], dim(mapping)[1]))
        } else{
          print(sprintf("##-- %s (total without NA: %s/%s) --##",
                        mvar, dim(mapping.na)[1], dim(mapping)[1]))
        }
        print(4)

        df[,mvar] <- mapping.na[df$SampleID, mvar]
        #form <- as.formula(sprintf("value ~ %s + %s + (1 | StudyID)", mvar, paste(rownames(subset(metadata, Go_lmemConfounder=="yes")), collapse="-")))
        form <- as.formula(sprintf("value ~ %s +  (1 | StudyID)", mvar))

        mod <- lmer(form, data=df,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
        ## lmer에서 control=은 "number of levels of each grouping~~" 오류가 있을때만 사용한다.
        ##
        # mod2 <- lmer(form, data=df)

        coef <- summary(mod)$coefficients
        coef <- coef[grep(mvar, rownames(coef)),,drop=F]

        res <- rbind(res, cbind(f, mvar, rownames(coef), coef))
        dim(res)
      }
    }


    #-- create table --#
    res <- as.data.frame(res)
    colnames(res) <- c("taxa", "metadata", "coefficient", "Estimate", "SE", "df", "t", "pvalue")
    res$pvalue <- as.numeric(as.character(res$pvalue))
    res$Estimate <- as.numeric(as.character(res$Estimate))
    res$SE <- as.numeric(as.character(res$SE))
    res$padj <- p.adjust(res$pvalue, method="fdr")
    res <- res[order(res$pvalue),]
    if (length(des) == 1) {
      res$des <- des
    }

    if (length(des) == 1) {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE, col.names = NA, file=sprintf("%s_%s/table/lmem/%s.%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], name, des, project, format(Sys.Date(), "%y%m%d"), sep="/"))
      }
      else {
        write.csv(res, quote = FALSE,col.names = NA, sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i],des, project, format(Sys.Date(), "%y%m%d"), sep="/"))
      }
    }
    else {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], name, project,format(Sys.Date(), "%y%m%d"), sep="/"))
      }
      else{
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s_%s/table/lmem/%s.%s.lmem.%s.csv",project, format(Sys.Date(), "%y%m%d"),  taxaname[i], project,format(Sys.Date(), "%y%m%d"), sep="/"))
      }
    }
  }
}
