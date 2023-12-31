

#===== Go_krakenLog

Go_krakenLog <- function(project,
                         hm.log,
                         st.log,
                         mpa){
  
  #===== Human read log
  hm.log.tab <- read.table(hm.log, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="");head(hm.log.tab)
  colnames(hm.log.tab) <- c("Total_reads","Human_reads", "Non_human_reads","Human_reads_%", "Non_human_reads_%")
  
  #===== Standeard read log
  st.log.tab <- read.table(st.log, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="");head(st.log.tab)
  colnames(st.log.tab) <- c("Hm_unclassified", "ST_classified", "ST_unclassified", "ST_classified_%", "ST_unclassified_%")
  
  for(col.name in c("Hm_unclassified")){
    st.log.tab[,col.name] <- NULL
  }
  
  log.tab <- merge(hm.log.tab, st.log.tab, by = "row.names")
  
  # Extract suffixes and get the unique names for log.tab
  suffixes <- gsub(".*_", "", log.tab$Row.names)
  most_common_suffix <- names(sort(table(suffixes), decreasing = TRUE))[1]
  cleaned_names <- gsub(paste0("_", most_common_suffix), "", log.tab$Row.names)
  cleaned_names <- make.unique(cleaned_names)
  rownames(log.tab) <- cleaned_names
  rownames(log.tab) <- gsub("-",  ".", rownames(log.tab))
  
  #===== read kraken mpa
  mpatable <- read.table(mpa, header=T, as.is=T, sep="\t", comment.char="", row.names=1, quote="");head(mpatable)
  
  mpatable.sel <- subset(mpatable, rownames(mpatable) %in% c("d__Archaea", "d__Bacteria", "d__Eukaryota|k__Fungi", "d__Viruses"))
  mpatable.sel.t <- t(mpatable.sel)
  
  # Extract suffixes and get the unique names for log.tab
  suffixes <- gsub(".*_", "", rownames(mpatable.sel.t))
  most_common_suffix <- names(sort(table(suffixes), decreasing = TRUE))[1]
  cleaned_names <- gsub(paste0("_", most_common_suffix), "", rownames(mpatable.sel.t))
  cleaned_names <- make.unique(cleaned_names)
  rownames(mpatable.sel.t) <- cleaned_names
  
  
  
  #===== final log.tab
  final.log.tab <- merge(log.tab,mpatable.sel.t, by = "row.names")
  
  # Set row names
  rownames(final.log.tab) <- final.log.tab$Row.names
  
  # Remove the "Row.names" column
  for(i in c("Row.names","Row.names")){
    final.log.tab[,i] <- NULL
  }
  
  out <- file.path(sprintf("%s", "1_out")) 
  if(!file_test("-d", out)) dir.create(out)
  
  write.csv(final.log.tab, quote = F,col.names = NA, row.names = T,
            file=sprintf("%s/%s.final.log.tab.%s.csv", out, project, format(Sys.Date(), "%y%m%d"),sep="/"))
  
}
  
  


