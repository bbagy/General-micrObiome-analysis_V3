
Go_importance_plot <- function(psIN, model, project,title,aggregate, MDA, name, bySample, height, width){
    
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  RF_model <- readRDS(model)
  tmp <- data.frame(RF_model$importance)
  
  if (aggregate == "YES" | aggregate == "Yes"| aggregate == "yes") {
    tmp$features <- rownames(tmp)
    tmp.sorted <- arrange(tmp, desc(MeanDecreaseAccuracy))
    RF.sorted.sel <- as.data.frame(subset(tmp.sorted, MeanDecreaseAccuracy > MDA))
    RF.sorted.sel$ShortName <- RF.sorted.sel$features
    
  } else if (aggregate == "NO" | aggregate == "No" | aggregate == "no") {
    RF <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF[,taxa] == "NA"
      RF[,taxa]<- as.character(RF[,taxa])
      RF[,taxa][is.na(RF[,taxa])] <- "__"
      for(i in 1:length(RF[,taxa])){
        if (RF[,taxa][i] == "s__" || RF[,taxa][i] == "g__" || RF[,taxa][i] == "f__" || RF[,taxa][i] == "o__" || RF[,taxa][i] == "c__"|| RF[,taxa][i] == "__"){
          RF[,taxa][i] <- ""
        }
      }
    }
    
    RF$ShortName <- paste(RF$Genus,"",RF$Species)
    for(taxa in c("Family", "Order", "Class","Phylum")){
      for(i in 1:length(RF[,taxa])){
        if (RF$ShortName[i] != "   "){
          next
        }      else if (RF$ShortName[i] == "   " & RF[,taxa][i] != ""){
          RF$ShortName[i] <- paste(RF[,taxa][i])
        }
      }
    }
    RF$features <- rownames(RF)
    RF.sorted <- arrange(RF, desc(MeanDecreaseAccuracy)  )
    RF.sorted.sel <- as.data.frame(subset(RF.sorted, MeanDecreaseAccuracy > MDA))
    
    for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
      RF.sorted.sel[,taxa] = NULL

    }
  }
  
  RF.sorted.sel.melt <- melt(RF.sorted.sel, id.vars=c("features","MeanDecreaseAccuracy","MeanDecreaseGini","ShortName"))
  

  if (bySample == "NO" |bySample == "No" | bySample == "no"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  } else if (bySample == "Yes" | bySample == "Yes"){
    p <- ggplot(data=RF.sorted.sel.melt, aes(x=reorder(features,MeanDecreaseAccuracy), y=value, fill=variable))+ geom_bar(stat="identity",position=position_dodge()) + geom_hline(yintercept=0)+ theme_classic()+ coord_flip() + scale_fill_brewer(palette="Set1")+scale_x_discrete(breaks = as.character(RF.sorted.sel.melt$features), labels = sprintf("%s",as.character(RF.sorted.sel.melt$ShortName))) + xlab("Taxa (Important features)") +ylab("Importance of the features") + theme(plot.title = element_text(hjust = 0.5))
  }


  
  p<- p+ ggtitle(sprintf("RF %s (MDA>%s) ", title, MDA)) 
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_1.RF.%s.%s.(%s).%s.pdf",out_path, project, name, MDA,format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_1.RF.%s.(%s).%s.pdf", out_path, project, MDA, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}





