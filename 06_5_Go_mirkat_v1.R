
Go_mirkat<- function(psIN, metaData, project, orders,name=NULL){
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
  
  
  # install package 
  packages <- c("data.table","CompQuadForm","devtools","ecodist","GUniFrac","GLMMMiRKAT","lme4","MASS","Matrix","MiRKAT","permute") 
  for (package in packages){
    if(!package %in% installed.packages()){
      install.packages(package)
    }else{library(package, character.only = TRUE)}
  }
  
  
  # install package from install_github
  github <- c("GLMMMiRKAT")
  if(!github %in% installed.packages()){
    install_github("hk1785/GLMM-MiRKAT", force=T)
  }else{library(package, character.only = TRUE)}
  
  #install.packages("data.table", version = "1.13.0")
  
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/mirkat",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  # check by metadata
  mapping <- data.frame(sample_data(psIN))
  for (mvar in  rownames(subset(metadata, Go_mirkat=="yes" | Confounder=="yes"))) {
    if (metadata[mvar, "type"] == "factor") {
      mapping[,mvar] <- factor(mapping[,mvar])
    } else if (metadata[mvar, "type"] == "numeric") {
      mapping[,mvar] <- as.numeric(as.character(mapping[[mvar]]))
    } else if (metadata[mvar, "type"] == "date") {
      mapping[,mvar] <- as.Date(sprintf("%06d", mapping[,mvar]), format="%m%d%y")
      mapping[,mvar] <- factor(as.character(mapping[,mvar]), levels=as.character(unique(sort(mapping[,mvar]))))
    }
  }
  
  sample_data(psIN) <- mapping
  
  
  res <- {}
  for (mvar in rownames(subset(metadata, Go_mirkat =="yes"))) {
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
    
    print(count(group.cbn))
    
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
      
      for (covar in rownames(subset(metadata, Confounder =="yes"))) {
        if (metadata[covar, "type"] == "factor") {
          if (mvar == covar){
            next
          }else{
            df.covar[,covar] <-as.numeric(mapping.cbn[,covar]  == as.character(unique(mapping.cbn[,covar] )[1]))
          }
          
        } else if (metadata[covar, "type"] == "numeric") {
          df.covar[,covar] <- mapping.cbn[,covar]
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
      
      
      if (length(rownames(subset(metadata, Confounder =="yes"))) >=1){
        for (covar in rownames(subset(metadata, Confounder =="yes"))) {
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
        
      }else if(length(rownames(subset(metadata, Confounder =="yes"))) ==0){
        
        permutation <- MiRKAT(y = df[,mvar], Ks = Ks, X = NULL, out_type = "D", method = "davies", 
                              omnibus = "permutation", returnKRV = FALSE, returnR2 = FALSE)
      }
      
      # create table
      per.df <- data.frame(unlist(permutation));per.df
      
      colnames(per.df) <- paste(basline,"vs",smvar);per.df
      per.df.t <- as.data.frame(t(per.df));per.df.t
      
      per.df.t$mvar <- mvar
      class(per.df.t)
      
      
      if (length(rownames(subset(metadata, Confounder =="yes"))) >=1){
        covars <- rownames(subset(metadata, Confounder =="yes"))[mvar != rownames(subset(metadata, Confounder =="yes"))]
        per.df.t$Confounder <- paste(setdiff(rownames(subset(metadata, Confounder=="yes")), "SampleType"), collapse="+")
        per.df.t$covar <- paste(setdiff(colnames(df.covar), "SampleType"), collapse="+")
        
      }else if(length(rownames(subset(metadata, Confounder =="yes"))) ==0){
        per.df.t <- per.df.t
      }
      
      
      res <- rbind(res, per.df.t);res
      print(res)
    }
  }
  
  if (length(name) == 1) {
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/mirkat.%s.%s.%s.csv",out_table, project,name, format(Sys.Date(), "%y%m%d"),sep="/"))
  }else{
    write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/mirkat.%s.%s.csv",out_table, project,format(Sys.Date(), "%y%m%d"),sep="/"))
  }
}
