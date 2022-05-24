#' Go_function2ps
#'
#'
#' @param function2ps
#' @keywords function2ps
#' @export
#' @examples
#' function2ps

Go_function2ps <- function(tabPath,
                           project=NULL, func.type,
                           name=NULL){
  # Read tab
  func.tab <- read_tsv(tabPath,col_types = cols())
  NumOFsample <- dim(func.tab)[2]
  
  # split a data frame
  otu <- as.matrix(func.tab[,3:NumOFsample])
  tax <- as.matrix(func.tab[,1:2])
  
  rownames(otu) <- tax[,1]
  rownames(tax) <- tax[,1]
  
  # define kegg or pathway
  if (any(grepl("K0", rownames(tax)))){
    colnames(tax) <- c("KO","KO.des")
    func <- "KEGG"
  }else if (any(grepl("PWY", rownames(tax)))){
    colnames(tax) <- c("pathway","path.des")
    func <- "pathway"
  }
  
  #merge phyloseq
  ps <- phyloseq(otu_table(otu, taxa_are_rows=T), tax_table(tax));ps
  
  print(ps)
  
  
  tt<- try(class(func.type),T)
  
  if(class(tt) == "try-error"){
    print("Please define the data type. PICRUSt or Hummann")
    break
  }else if(any(grepl(func.type, c("picrust","Picrust","Picrust2","PICRUSt","PICRUSTt2")))){
    func.type <- "PICRUSTt2"
  }else if(any(grepl(func.type, c("Human","Humann","humann","human","Human2","Humann2","humann2","human2")))){
    func.type <- "Humann2"
  }
  
  
  # saving file
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)
  saveRDS(ps, sprintf("%s/ps.%s.%s.%s.%s%s.rds",rds,
                      func.type,
                      func,
                      project,
                      ifelse(is.null(name), "", paste(name, ".", sep = "")), 
                      format(Sys.Date(), "%y%m%d")))
}
