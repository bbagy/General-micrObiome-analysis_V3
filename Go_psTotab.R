Go_psTotab <- function(psIN, project){
  
  
  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  #====== step 1 read ps object
  seqtab.nochim <- as.matrix(otu_table(psIN))
  tax <- as.matrix(tax_table(psIN))
  
  
  #====== step 2 extract fna
  seqs <- getSequences(t(seqtab.nochim))
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  

  
  if (nchar(headers[1]) < 100){
    seqs <- getSequences(seqtab.nochim)
    headers <- paste(">", seqs, sep="")
    fasta <- c(rbind(headers, seqs))
    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  }else{
    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  }
  
  

  
  #====== step 3 get the table
  otu <- as.data.frame(t(otu_table(psIN)));dim(otu)
  tax <- tax_table(psIN);dim(tax)
  

  tt <- try( otuTable <- cbind(otu,tax),T)
  if(class(tt) == "try-error"){
    otu <- as.data.frame(otu_table(psIN));dim(otu)
    tax <- tax_table(psIN);dim(tax)
    otuTable <- cbind(otu,tax)
    
    seqs <- getSequences(t(seqtab.nochim))
    headers <- paste(">", seqs, sep="")
    fasta <- c(rbind(headers, seqs))
    
    write(fasta, file=sprintf("%s/%s.%s.psTotab.seqs.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
    
  }else{
    otuTable <- cbind(otu,tax)
  }
  
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s.%s.psTotab.ASVs.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
}



