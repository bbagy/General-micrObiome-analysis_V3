

Go_deseq2fishtaco <- function(psIN, project, file_path){
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)
  out_fishtaco <-file.path(sprintf("%s_%s/table/deseq2/funTabforFishtaco",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_fishtaco)) dir.create(out_fishtaco)
  
  
  
  
  path <- file_path
  print(path)
  filenames <- list.files(path, pattern="Forfishtaco.csv")
  print(filenames)
  
  for(f in filenames){
    top.deseq2 <- read.csv(sprintf("%s/%s",path,f),row.names=1,check.names=FALSE)
    print(head(top.deseq2))
    ranks <-c("KO", "KO.des","Genus","Species","Path","Path.des")
    # otu table
    otu.filt <- as.data.frame(otu_table(psIN)) 
    
    otu.filt$ko <- getTaxonomy(otus=rownames(otu.filt), tax_tab=tax_table(psIN), taxRanks=ranks,  level="KO")
    agg <- aggregate(. ~ ko, otu.filt, sum)
    
    
    agg.clean <- agg[agg$ko %in% top.deseq2$x, ]# agg[!agg$ko %in% top.deseq2$x, ]  top.deseq2$x만 제거
    
    write.csv(agg.clean, quote = FALSE,col.names = NA,file=sprintf("%s/%s",out_fishtaco,f,sep="/"))
  }
}







