Go_dist_plot <- function(psIN, project, distance_metrics, distance, group,orders, name,height,width,plot = TRUE) {
    
  if(!is.null(dev.list())) dev.off()
   
  colorset = "Paired" # Dark1 Set1
  

  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  if (plot == TRUE) {
    if (length(name) == 1) {
      pdf(sprintf("%s/5_distplot.%s.%s.%s.pdf", out_path, project, name,format(Sys.Date(), "%y%m%d")),height = height, width=width)
    } else{
      pdf(sprintf("%s/5_distplot.%s.%s.pdf", out_path, project, format(Sys.Date(), "%y%m%d")),height = height, width=width)
    }
  
  # calc distances
  for (dist in distance_metrics){
    wu <- as.dist(distance[[dist]])
    wu.m <- melt(as.matrix(wu))
    
    # remove self-comparisons
    wu.m <- wu.m %>%
      filter(as.character(Var1) != as.character(Var2)) %>%
      mutate_if(is.factor, as.character)
    
    
    wu.m  = wu.m  %>%
      rowwise() %>%      # for each row
      mutate(Samples = paste(sort(c(Var1, Var2)), collapse = "-")) %>%  # sort the teams alphabetically and then combine them separating with -
      ungroup()
    
    wu.m.sel  = distinct(wu.m, Samples, .keep_all=T)
    
    
    # get sample data (S4 error OK and expected)
    mapping <- data.frame(sample_data(psIN))
    mapping$ID <- as.character(rownames(mapping))
    mapping[,group] <- as.character(mapping[,group])

    
    sd <- mapping %>%
      select("ID", group) %>%
      mutate_if(is.factor,as.character)
    
    # combined distances with sample data
    # sample1
    colnames(sd) <- c("Var1", "Type1")
    wu.m.sel$Var1 <- factor(wu.m.sel$Var1)
    sd$Var1 <- factor(sd$Var1)
    wu.sd <- left_join(wu.m.sel, sd, by = "Var1")
    # sample2
    wu.sd$Var2 <- as.factor(wu.sd$Var2)
    colnames(sd) <- c("Var2", "Type2")
    sd$Var2 <- factor(sd$Var2)
    wu.sd <- left_join(wu.sd, sd, by = "Var2")
    
    wu.sd$Type3 <- ifelse(wu.sd$Type1 == wu.sd$Type2, wu.sd$Type1,"across_group")
    wu.sd$Type3 <- factor(wu.sd$Type3)
    # make a combination for stat
    wu.sd.sel <- subset(wu.sd, Type3 != "across_group")
    wu.sd.sel$Type3 <- factor(wu.sd.sel$Type3)

    # cbn <- combn(x = levels(wu.sd$Type3), m = 2)
    baseline <- "across_group"
    cbn <-{}
    for (x in levels(wu.sd.sel$Type3)){
      cbn <- cbind(cbn, c(baseline, x))
    }
    
    cbn.sel<-cbn[, !duplicated(t(cbn))]
    
    my_comparisons <- {}
    for(i in 1:ncol(cbn.sel)){
      x <- cbn.sel[,i]
      my_comparisons[[i]] <- x
    };my_comparisons
    
    
    # plot
    if (length(orders) >= 1) {
      wu.sd$Type2 <- factor(wu.sd$Type2, levels = orders)
    }       else {
      wu.sd$Type2 <- factor(wu.sd$Type2)
    }
    
    p <- ggplot(wu.sd, aes(x = Type3, y = value, colour=Type3)) + theme_classic()+
      geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16, alpha = 0.3,position=position_jitter(0.2)) +
      scale_color_brewer(palette=colorset) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      #facet_wrap(~ Type2, scales = "free_x") +
      ylab(dist) + xlab(NULL) +  theme(legend.position="none")+
      stat_compare_means(method= "wilcox.test", label = "p.format", comparisons = my_comparisons, size = 2.5)
    
    if (length(name) == 1) {
      p<- p+ ggtitle(sprintf("%s_%s",dist,name))
    } else{
      p <- p+ ggtitle(sprintf("%s",dist))
    }
    print(p)
  }
    dev.off()
    
  } else {
    return(wu.sd)
  }
}
