#' A Go_clme
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords Taxa barplots
#' @export
#' @examples
#' Go_clme()


Go_clme <- function(psIN, metaData, project, paired, node, decreasing, height,timepoint,ID, orders,xangle, name, width, plotCols, plotRows){
    
    
    if(!is.null(dev.list())) dev.off()
    
  
  alpha_metrics = c("Chao1","Shannon")
  
  colorset = "Dark2" # Dark1 Set1 Paired
  
  # Descriptions 분석 하고자 하는 variation에 subgroup
  # paired 환자나 같은 사람 ID
  # node 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # decreasing 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  
  
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))
  
  if (length(name) == 1) {
    pdf(sprintf("%s_%s/pdf/3_clme.%s.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,name,node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  else {
    pdf(sprintf("%s_%s/pdf/3_clme.%s.(%s.%s).%s.pdf",project, format(Sys.Date(), "%y%m%d"),project,node,decreasing, format(Sys.Date(), "%y%m%d")), height = height, width = width)
  }
  
  
  
  # adiv
  adiv <- estimate_richness(psIN, measures=alpha_metrics)
  mapping <-data.frame(sample_data(psIN))
  rownames(adiv) <- gsub("^X", "", rownames(adiv))
  adiv$SampleID <- rownames(adiv)
  rownames(adiv) <- rownames(mapping)
  adiv <- merge(adiv, mapping, by="row.names"); rownames(adiv) <- adiv$SampleID

  if (length(orders) >= 1) {
    adiv[,timepoint] <- factor(adiv[,timepoint], levels = orders)
  }

  
  
  # clme
  cons <- list(order = "umbrella" , node=node, decreasing = decreasing) 
  # 전반적인 패턴을 보고 가장 높은 time point 에 node를 설정  
  # 패턴에 증가 하는지 감소 하는지 판단 하고 decreazing = true and false 를 판단, mean and median 으로 판단
  
  print(cons)
  
  plotlist <- list()
  for (mvar in rownames(subset(metadata, Go_clme == "yes"))) {
    print(mvar)
    
    if (length(unique(adiv[,mvar])) < 2){
      next
    }
    
    
    
    # Na 제거

    adiv[,mvar] <- data.frame(adiv[,mvar]);adiv[,mvar]
    adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
    adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]
    
    
    # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
    adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
    adiv <- adiv.na
    
    for (des in unique(adiv[,mvar])){
      
      if(dim(subset(adiv, adiv[,mvar] == des))[1] < 3){
        next
      }
      
      if(timepoint == mvar){
        next
      }
      
      print(des)
      for (am in alpha_metrics){
        form <-as.formula(sprintf("%s ~ %s + (1|%s)" , am, timepoint, paired))
        
        clme.mod <- clme(form, data = adiv[adiv[,mvar] == des,], constraints = cons, seed = 2, nsim = 1000)
        clme.sum <- summary(clme.mod, seed=2)
        clme.globalp <- function(model) { label <- substitute(
          italic(p) == globalp,
          list(globalp <- model$p.value) )
        as.character(as.expression(format(globalp, nsmall=3))) 
        }
        
        clme.globalp <- paste("CLME P=",clme.globalp(clme.sum))
        
        # plot
        p <- ggplot(adiv[adiv[,mvar]==des,], mapping = aes_string(x=timepoint, y=am, color=timepoint, group=paired)) + geom_line(color="grey") + geom_point(size = 1.25) + xlab(timepoint) + ylab(sprintf("%s Index\n", am)) + ggtitle(sprintf("%s-%s \n (%s) ", mvar, des, clme.globalp))  + scale_color_brewer(palette=colorset)+theme_bw() +theme(title=element_text(size=8), axis.text.x=element_text(angle=xangle,hjust=1,vjust=0.5)) + theme(legend.position= "NONE" )
        
        if (length(ID) == 1) {
          p= p + geom_text_repel(aes_string(label = ID), size = 2)
        }
        plotlist[[length(plotlist)+1]] <- p
      }
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}


