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
}