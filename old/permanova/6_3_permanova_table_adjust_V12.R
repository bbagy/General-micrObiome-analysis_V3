#' A Go_perm
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Beta diversity Adonis test (PERMANOVA)
#' @export
#' @examples
#' Mar 07 2020
#' adjsted 기능을 추가 하였다.
#' 분석 할때 마다 수치가 조금 변하는 것을 수정 하였다.set.seed(1)
#' dm를 따로 분리 하여 시간을 단축 하였고, dm를 다른 방법으로 분석 할수 있게 되었다.
#' Go_perm()


Go_perm <- function(psIN, metaData, project, distance,distance_metrics, adjust,des, name){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s/table",out)) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_perm <- file.path(sprintf("%s/perm",out_path)) 
  if(!file_test("-d", out_perm)) dir.create(out_perm)
  
  out_distance <- file.path(sprintf("%s/distance",out_path)) 
  if(!file_test("-d", out_distance)) dir.create(out_distance)
  
  
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  # Run
  if (length(des) == 1) {
    # Uni
    print(sprintf("#--- Running Paired-PERMANOVA (%s) ---#", des))
  }
  else {
    print("#--- Running Paired-PERMANOVA  ---#")
  }
  
  set.seed(1)
  mapping.sel <-data.frame(sample_data(psIN))
  res.pair <-{}
  for (mvar in rownames(subset(metadata, Go_perm =="yes"))) {
    if (length(unique(as.character(mapping.sel[,mvar]))) == 1) {
      next
    }
    for (distance_metric in distance_metrics) {
      mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
      if (length(unique(mapping.sel.na[,mvar])) == 1)
        next
      
      psIN.sel <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
      
      
      distance <- Go_dist(psIN = psIN.sel, project = project, distance_metrics = distance_metric)
      # save distance
      distance.tab <- as.matrix(distance(psIN, method=distance_metric, type="samples"))
  
      if (length(des) == 1) {
        if (length(name) == 1) {
          print(1)
          write.csv(distance.tab, quote = FALSE,col.names = NA,file=sprintf("%s/distanceTab.%s.%s.%s.%s.%s.csv",out_distance, project,mvar, des, name, distance_metric,format(Sys.Date(), "%y%m%d"),sep="/"))
        }
        else {
          write.csv(distance.tab, quote = FALSE,col.names = NA,file=sprintf("%s/distanceTab.%s.%s.%s.%s.csv",out_distance, project, mvar,des, distance_metric,format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      }
      else{
        if (length(name) == 1) {
          print(2)
          write.csv(distance.tab, quote = FALSE,col.names = NA,file=sprintf("%s/distanceTab.%s.%s.%s.%s.csv",out_distance, project,mvar,name, distance_metric,format(Sys.Date(), "%y%m%d"),sep="/"))
        }
        else {
          write.csv(distance.tab, quote = FALSE,col.names = NA,file=sprintf("%s/distanceTab.%s.%s.%s.csv",out_distance, project,mvar,distance_metric, format(Sys.Date(), "%y%m%d"),sep="/"))
        }
      }
      
      
      
      
      
      
      
      
      
      

      if (length(adjust) >= 1) {
        if(mvar == adjust )
          next
        
        # integer control (add 20201029)
        if (class(mapping.sel.na[,mvar]) == "character"){
          mapping.sel.na[,mvar] <- factor(mapping.sel.na[,mvar])
          sample_data(psIN.sel) <- mapping.sel.na
        }
        if (class(mapping.sel.na[,mvar]) == "integer" | class(mapping.sel.na[,mvar]) == "numeric"){
          mapping.sel.na[,mvar] <- factor(mapping.sel.na[,mvar])
          sample_data(psIN.sel) <- mapping.sel.na
        }
        
        
      pair.ado <- pairwise.adonis.ad(x=as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar], adjust=adjust, map=mapping.sel.na, mvar=mvar)
      
      
      } else{
        pair.ado <- pairwise.adonis(as.dist(distance[[distance_metric]]), factors = mapping.sel.na[,mvar],map=mapping.sel.na, mvar=mvar)
      }
      
      tmp <- as.data.frame(pair.ado)
      tmp$distance_metric <- distance_metric
      tmp$mvar <- mvar
      tmp$adjusted <- paste(setdiff(adjust, "SampleType"), collapse="+")
      res.pair <- rbind(res.pair, tmp)
    }
  }
  
  
  if (length(adjust) >= 1) {
    if (length(des) == 1) {
      if (length(name) == 1) {
        print(1)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.%s.csv",out_perm, project, des, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.csv",out_perm, project, des, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
    else{
      if (length(name) == 1) {
        print(2)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.%s.csv",out_perm, project,name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.adjusted.%s.%s.csv",out_perm, project, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  } else{
    if (length(des) == 1) {
      if (length(name) == 1) {
        print(3)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.%s.csv",out_perm, project, des, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.csv",out_perm, project, des, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
    
    else{
      if (length(name) == 1) {
        print(4)
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.%s.csv",out_perm, project, name, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
      else {
        write.csv(res.pair, quote = FALSE,col.names = NA,file=sprintf("%s/pair_permanova.%s.%s.csv",out_perm, project, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  }

  return(res.pair)
}




