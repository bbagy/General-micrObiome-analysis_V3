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
