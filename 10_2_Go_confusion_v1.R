

Go_confusion <- function(model, project,legend, name, height, width){
    
  if(!is.null(dev.list())) dev.off()
   
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  rf_model <- readRDS(model)
  
  # run
  conf_mat <- rf_model$confusion[, -ncol(rf_model$confusion)]
  melt_conf_mat <- reshape2::melt(conf_mat, na.rm = TRUE)
  table <- data.frame(melt_conf_mat);table
  
  plotTable <- table %>%
    mutate(goodbad = ifelse(table$Var1 == table$Var2, "good", "bad")) %>%
    group_by(Var2) %>%
    mutate(prop = value/sum(value))
  
  
  # plot
  p <- ggplot(data = plotTable, mapping = aes(x = Var2, y = Var1, fill = goodbad, alpha = prop)) +
    geom_tile() +
    geom_text(aes(label = value), vjust = .5, fontface  = "bold", alpha = 1) +
    scale_fill_manual(values = c(good = "blue", bad = "red")) +
    theme_bw() + xlab("True class") + ylab("Predicted class") + theme(legend.position=legend, plot.title = element_text(hjust = 0.5))
    #xlim(rev(levels(table$Var1))) 

  p <- p+ ggtitle(sprintf("confusionMatrix(%s) ",name)) 
  
  # out file
  if (length(name) == 1) {
    pdf(sprintf("%s/10_2.confusion.%s.%s.%s.pdf",out_path, project, name, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  } else {
    pdf(sprintf("%s/10_2.confusion.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  print(p)
  dev.off()
}




