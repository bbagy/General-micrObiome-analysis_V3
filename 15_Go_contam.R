#' A Go_contam
#'
#'
#' @param Go_contam
#' @keywords Go_contam
#' @export
#' @examples
#' Go_contam


Go_contam <- function(psIN, project, 
                      col.name = "StudyID", 
                      neg.con = "neg",
                      taxa=NULL) {


  out <- file.path("1_out") 
  if(!file_test("-d", out)) dir.create(out)
  
  rds <- file.path("2_rds") 
  if(!file_test("-d", rds)) dir.create(rds)

  map <- data.frame(sample_data(psIN)) 
  
  #====== step 1 define contaminated ps and ucontaminated ps
  contams <- subset_samples(psIN, map[,col.name] %in% neg.con);contams
  samples <- subset_samples(psIN, !(map[,col.name] %in% neg.con));samples

  #====== step 2 substrate ucontaminated ps -  contaminated ps
  contams.n1 <- prune_taxa(taxa_sums(contams)>=1, contams) 
  samples.n1 <- prune_taxa(taxa_sums(samples)>=1, samples) 
  allTaxa <- names(sort(taxa_sums(psIN),TRUE))
  negtaxa <- names(sort(taxa_sums(contams.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  ps.decontam <- prune_taxa(taxa.noneg,samples.n1)

  if(!is.null(taxa)){
    seqtab.nochim <- as.matrix(otu_table(ps.decontam))
    tax <- as.matrix(tax_table(ps.decontam))

    is.a <- tax[,"Order"] %in% taxa 
    seqtab.a <- seqtab.nochim[,!is.a];dim(seqtab.a)
    tax.a <- tax[!is.a,]

    ps.decontam <- phyloseq(otu_table(seqtab.a, taxa_are_rows=FALSE), tax_table(tax.a));ps.decontam
    sample_names(ps.decontam)
    sample_names(ps.decontam) <- gsub("X","",sample_names(ps.decontam));sample_names(ps.decontam)
  }else{
    ps.decontam <- ps.decontam
  }


  #====== step 3 get decontam.fna
  seqs <- getSequences(ps.decontam)
  headers <- paste(">", seqs, sep="")
  fasta <- c(rbind(headers, seqs))
  
  write(fasta, file=sprintf("%s/%s.%s.decontam.fna",out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
  #====== step 4 get the table
  otu <- as.data.frame(t(otu_table(ps.decontam)));dim(otu)
  tax <- tax_table(ps.decontam);dim(tax)
  
  otuTable <- cbind(otu,tax)
  
  write.csv(otuTable, quote = FALSE,col.names = NA,#row.names = FALSE, 
            file=sprintf("%s/%s.%s.asvTable_decontam.csv",out,project,format(Sys.Date(), "%y%m%d"), sep="/"))
  
  saveRDS(ps.decontam, sprintf("%s/ps_decontam.%s.%s.rds", rds, project,format(Sys.Date(), "%y%m%d")))
  
  
  print(psIN)
  print(ps.decontam)  
  
  return(ps.decontam)
}