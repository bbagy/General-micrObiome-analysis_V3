Go_emptyMap <- function(psIN, project){
  
  # out dir
  map <- file.path("3_map") 
  if(!file_test("-d", map)) dir.create(map)
  
  
  # empty map table
  SampleID <- sample_names(psIN)
  StudyID <- sample_names(psIN)
  emptyMap <- data.frame(SampleID, StudyID)
  
  cat(sprintf("empty map is saved in %s.\n",map))
  cat("                                                       \n")
  write.csv(emptyMap, quote = FALSE, col.names = NA, row.names = F,
            file=sprintf("%s/emptyMap.%s.%s.csv",map, project,format(Sys.Date(), "%y%m%d"),sep="/"))
  
  
  # empty metadata table
  column.names <- c("StudyID", "Variation1", "Variation2","etc")
  col.count <- length(column.names)
  
  # 	"Go_overview","Go_ancombc","Go_deseq2","Go_box","Go_bdiv",	"Go_barchart","Go_linear","Go_clme","Go_perm",
  analysis <- c("type",	"baseline",	"Go_reg", "Go_mirkat", "Go_lmem","Confounder")

  row.count <- length(analysis)
  
  emptyMetadata <- data.frame(matrix(ncol = col.count, nrow = row.count))
  colnames(emptyMetadata) <- column.names
  rownames(emptyMetadata) <- analysis


  for(an in analysis){
    if (an == "type"){
      emptyMetadata[c(an), ] <- c("", "factor", "numeric", "factor")
    }else if(an == "baseline"){
      emptyMetadata[c(an), ] <- c("", "control", "before", "male")
    }else{
      emptyMetadata[c(an), ] <- c("no", "no", "yes", "yes")
    }
  }
  
  #cat(sprintf("empty metadata is saved in %s.\n",map))
  #cat("                                                       \n")
  #write.csv(emptyMetadata, quote = FALSE, col.names = NA,  row.names = T,
  #          file=sprintf("%s/emptyControlpanel.%s.%s.csv",map, project,format(Sys.Date(), "%y%m%d"),sep="/"))
} 








