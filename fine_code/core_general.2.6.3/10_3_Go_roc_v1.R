

Go_roc <- function(model, project,map, categorical, name, height, width){
    
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  if (length(name) == 1) {
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.%s.pdf", out_path, name,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }else{
    pdf(sprintf("%s/10_3.ROC.%s.%s.%s.pdf", out_path,categorical, project,format(Sys.Date(), "%y%m%d")),height = height, width=width)
  }
  
  rf_model <- readRDS(model)
  
  
  pred <- predict(rf_model, type="prob")
  
  # map 정리 
  sel <- intersect(rownames(map), rownames(pred))
  map <- map[sel,, drop=F]
  #dim(pred);dim(map); class(pred);class(map);
  
  
  #dim(map.sel.sel)
  par(mfrow = c(1,length(colnames(pred))))
  for (i in 1:length(colnames(pred))){
    print(i)
    pred2 <- prediction(pred[,i], as.matrix(map[,categorical]))# or oedered(map[,categorical])
    perf <- performance(pred2, "tpr", "fpr")
    perf.auc <- performance(pred2, "auc")
    pred_df <- data.frame(SampleID=rownames(pred), predicted=colnames(pred)[apply(pred, 1, function(x) which.max(x))], true = map[,categorical], stringsAsFactors=F); 
    pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
    confusion_matrix <- table(pred_df[, c("true", "predicted")])
    accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
    vec.pred <- as.numeric(pred_df$predicted)-2; vec.true <- as.numeric(pred_df$true)-2 # it was "-1"
    mccvalue <- mcc(vec.pred, vec.true)
    if (i == 1){
      color = "red"
    }else if(i == 2){
      color = "blue"
    }else if(i == 3){
      color = "green"
    }
    
    if (length(name) == 1) {
      plot(perf, main=sprintf("RF_ROC %s (%s)", colnames(pred)[i], name), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }else{
      plot(perf, main=sprintf("RF_ROC %s", colnames(pred)[i]), col = color) + text(x=0.7, y=0.1, label=sprintf("mean AUC=%.4g\n accuracy=%.2f%%\n MCC=%.4g", unlist(perf.auc@y.values), accuracy, mccvalue))
      abline(0,1, col="grey")
    }
  }
  dev.off()
}




