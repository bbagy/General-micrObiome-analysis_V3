Go_regression <- function(df, metaData, project, orders, outcomes, pvalue=0.05, des,name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/regression",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  
  #meta data
  metadataInput <- read.csv(sprintf("%s",metaData),header=T,as.is=T,row.names=1,check.names=F)
  metadata <- as.data.frame(t(metadataInput))

  
  # data control
  
  # fix column types
  adiv <- data.frame(df)

  
  for (mvar in  rownames(subset(metadata, Go_reg=="yes" | Go_regConfounder=="yes"))) {
    if (metadata[mvar, "type"] == "factor") {
      adiv[,mvar] <- factor(adiv[,mvar])
      if (!(is.na(metadata[mvar, "baseline"])) && metadata[mvar, "baseline"] != "") {
        adiv[,mvar] <- relevel(adiv[,mvar], metadata[mvar, "baseline"])
      }
    } else if (metadata[mvar, "type"] == "numeric") {
      adiv[,mvar] <- as.numeric(as.character(adiv[[mvar]]))
    } else if (metadata[mvar, "type"] == "date") {
      adiv[,mvar] <- as.Date(sprintf("%06d", adiv[,mvar]), format="%m%d%y")
      adiv[,mvar] <- factor(as.character(adiv[,mvar]), levels=as.character(unique(sort(adiv[,mvar]))))
    }
  }


  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes){
    
    if (metadata[outcome, "type"] == "factor") {
      adiv[,outcome] <- factor(adiv[,outcome])
      if (!(is.na(metadata[outcome, "baseline"])) && metadata[outcome, "baseline"] != "") {
        adiv[,outcome] <- relevel(adiv[,outcome], metadata[outcome, "baseline"])
      }
    } else if (metadata[outcome, "type"] == "numeric") {
      adiv[,outcome] <- as.numeric(adiv[,outcome])
    } else if (metadata[outcome, "type"] == "date") {
      adiv[,outcome] <- as.Date(sprintf("%06d", adiv[,outcome]), format="%m%d%y")
      adiv[,outcome] <- factor(as.character(adiv[,outcome]), levels=as.character(unique(sort(adiv[,outcome]))))
    }
    
    res <- {}
    for (mvar in rownames(subset(metadata, Go_reg =="yes"))) {
      if (outcome == mvar | outcome == "Chao1" & mvar == "Shannon" | outcome == "Shannon" & mvar == "Chao1") {
        next
      }
      

      # NA 제거
      adiv[,mvar] <- as.character(adiv[[mvar]]);adiv[,mvar]
      adiv[,mvar][adiv[,mvar]==""] <- "NA";adiv[,mvar]
      #adiv[,mvar]<- as.factor(adiv[,mvar]);adiv[,mvar]
      # adiv.na <- adiv[!(is.na(adiv[,mvar])), ];adiv.na[,mvar] 틀린건 없는 거 같은데 지워지지 않는다. 
      adiv.na <- subset(adiv, adiv[,mvar] != "NA");adiv.na[,mvar]  # subset 를 사용한 NA 삭제
      
      print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                    mvar, dim(adiv.na)[1], dim(adiv)[1]))
      
      if (length(unique(adiv.na[,mvar])) ==1) {
        next
      }
      
      # column 정리

      if (metadata[mvar, "type"] == "factor") {
        adiv[,mvar] <- factor(adiv[,mvar])
        if (metadata[mvar, "baseline"] != "") {
          adiv[,mvar] <- relevel(adiv[,mvar], metadata[mvar, "baseline"])
        } 
      } else if (metadata[mvar, "type"] == "numeric") {
        adiv[,mvar] <- as.numeric(as.character(adiv[,mvar]))
      } else if (metadata[mvar, "type"] == "date") {
        adiv[,mvar] <- as.Date(sprintf("%06d", adiv[,mvar]), format="%m%d%y")
        adiv[,mvar] <- factor(as.character(adiv[,mvar]), levels=as.character(unique(sort(adiv[,mvar]))))
      }
      
      if (length(rownames(subset(metadata, Go_regConfounder =="yes"))) >= 1){
        regConfounder <- rownames(subset(metadata, Go_regConfounder =="yes"))[mvar != rownames(subset(metadata, Go_regConfounder =="yes"))]
        form <- as.formula(sprintf("%s ~ %s + %s", outcome, mvar, paste(setdiff(regConfounder, "SampleType"), collapse="+")))
        print(form)
        print(1)
      } else{
        form <- as.formula(sprintf("%s ~ %s", outcome, mvar))
        print(form)
        print(3)
      }
      
      
      if (class(adiv[,outcome]) == "numeric"){
        mod <- lm(form, adiv)  # lm or glm or lmer
        m <- "lm"
      } else if (class(adiv[,outcome]) == "factor"){
        mod <- glm(form, adiv[adiv[,outcome] %in% levels(adiv[,outcome]),],  family = binomial(link='logit'))
        m <- "glm"
      }
      
      summary(mod)
      coef <- as.data.frame(summary(mod)$coefficients)
      coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
      colnames(coef) <- c("Estimate", "SE", "t", "pval")

      if (dim(coef)[1] == 0){
        next
      }
      
      coef$outcome <- outcome
      coef$mvar <- mvar
      coef$model <- m
      
      if (length(rownames(subset(metadata, Go_regConfounder =="yes"))) >= 1){
        coef$multi <- paste(setdiff(rownames(subset(metadata, Go_regConfounder=="yes")), "SampleType"), collapse="+")
        type <-"multi"
      }else{
        type <-"uni"
      }
      
      
      
      res <- rbind(res, coef)
      
    }
    
    if (length(des) == 1) {
      res$des <- des
    }
    res$padj <- p.adjust(res$pval, method="fdr")
    #res <- res[order(res$time_point),]
    res$comp <- factor(rownames(res), levels=rownames(res))
    res$dir <- ifelse(res$pval < pvalue, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
    
    if (length(des) == 1) {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.%s.csv",out_table, project,outcome, des,name, type, format(Sys.Date(), "%y%m%d"),sep="/"))
      }else{
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.csv",out_table, project,outcome, des,type, format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }else{  
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.%s.csv",out_table, project,outcome, name, type, format(Sys.Date(), "%y%m%d"),sep="/"))
      } else{ 
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/5_6_regression_%s.%s.%s.csv",out_table, project,outcome,type,format(Sys.Date(), "%y%m%d"),sep="/"))
      }
    }
  }
}






