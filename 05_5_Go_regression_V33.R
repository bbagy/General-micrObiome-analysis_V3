
Go_regression <- function(data, project,  
                          outcomes, 
                          uni.vars=NULL, 
                          mul.vars=NULL,
                          interaction=NULL, 
                          orders, pvalue=0.05, name=NULL){
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)
  out_table <- file.path(sprintf("%s_%s/table/regression",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_table)) dir.create(out_table)
  
  
  # data control
  
  # fix column types
  if(!is.null(uni.vars) | !is.null(mul.vars)){
    
    varis <- unique(c(uni.vars,mul.vars))
    
    for (mvar in  varis) {
      if (class(data$mvar) == "character") {
        data[,mvar] <- factor(data[,mvar])
        # setting baseline
        data[,mvar] <- factor(data[,mvar], levels = intersect(orders, data[,mvar]))
        # NA 제거
        data[,mvar] <- as.character(data[[mvar]]);data[,mvar]
        data[,mvar][data[,mvar]==""] <- "NA";data[,mvar]
        #data.na <- subset(data, data[,mvar] != "NA");data.na[,mvar]  # subset 를 사용한 NA 삭제
        
      } else if (class(data$mvar) == "numeric" | min(mvar) < 0) {
        # NA 제거
        data[,mvar] <- as.character(data[[mvar]]);data[,mvar]
        data[,mvar][data[,mvar]==""] <- "NA";data[,mvar]
        #data.na <- subset(data, data[,mvar] != "NA");data.na[,mvar]  # subset 를 사용한 NA 삭제
        data[,mvar] <- as.numeric(as.character(data[[mvar]]))
      } 
    }
  }



  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes){
    if (length(unique(data[,outcome])) > 3){
      
      
      
    }else{
      if (class(outcome) == "character") {
        data[,outcome] <- factor(data[,outcome])
        data[,outcome] <- factor(data[,outcome], levels = intersect(orders, data[,outcome]))
        
        # NA 제거
        data[,outcome] <- as.character(data[[outcome]]);data[,outcome]
        data[,outcome][data[,outcome]==""] <- "NA";data[,outcome]
        #data.na <- subset(data, data[,outcome] != "NA");data.na[,outcome]  # subset 를 사용한 NA 삭제
        # set the baseline for outcome
        
        if(length(unique(data[,outcome])) == 2){
          
          data[,outcome] <- factor(data[,outcome])
          out <- levels(data[,outcome])[1]
          
          
          data[,outcome] <- factor(ifelse(data[,outcome]== levels(data[,outcome])[1],0,1), levels=c(0,1), 
                                   labels = levels(data[,outcome]))
          
          
          print(levels(data[,outcome]))
        }
      } else if (class(outcome)  == "numeric") {
        
        # NA 제거
        data[,outcome] <- as.character(data[[outcome]]);data[,outcome]
        data[,outcome][data[,outcome]==""] <- "NA";data[,outcome]
        # data <- subset(data, data[,mvar] != "NA");data[,outcome]  # subset 를 사용한 NA 삭제
        data[,outcome] <- as.numeric(as.character(data[[outcome]]))
        
        data[,outcome] <- as.numeric(data[,outcome])
      } 
      
      res <- {}
      
      for (mvar in uni.vars) {
        
        if (outcome == mvar) {
          next
        }
        
        # print(sprintf("##-- %s (total without NA: %s/%s) --##", 
        #              mvar, dim(data.na)[1], dim(data)[1]))
        
        if (length(unique(data[,mvar])) ==1) {
          next
        }
        
        # get formula
        if(!is.null(uni.vars) & is.null(mul.vars) & is.null(interaction)){
          form <- as.formula(sprintf("%s ~ %s", outcome, mvar))
          print("Univariate anaysis")
          type <- "uni"
        } else if (!is.null(mul.vars)){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(mul.vars, interaction)
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+")))
            print("Multivariate anaysis")
            type <- "multi-interantion"
            
          }else{
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars, "SampleType"), collapse="+")))
            
            print("Multivariate anaysis with interaction effect")
            type <- "multi"
          }
        } 
        
        print(form)
        
        
        # model
        if (class(data[,outcome]) == "numeric"){
          m <- "Regression (glm-poisson)"
          mod <- glm(form, data=data,  family = poisson(link='log'))
        } else if (length(unique(data[,outcome])) == 2){
          m <- "Logistic regression (glm-binomial)"
          mod <- glm(form, data=data,  family=binomial("logit"))
        } 
        print(m)
        
        
        # out for the model
        coef <- as.data.frame(summary(mod)$coefficients)
        coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
        colnames(coef) <- c("Estimate", "SE", "t", "pval")
        
        if (dim(coef)[1] == 0){
          next
        }
        
        
        
        # out for the confidence interval 
        conf <- data.frame(confint(mod))
        conf <- conf[setdiff(rownames(conf), "(Intercept)"),,drop=F]
        conf.na <- na.omit(conf) 
        
        coef$`2.5 %` <- conf.na$`2.5 %`
        coef$`97.5 %` <- conf.na$`97.5 %`
        
        coef$outcome <- outcome
        coef$mvar <- mvar
        coef$model <- m
        coef$deviance <- pchisq(q=mod$null.deviance-mod$deviance,df=mod$df.null-mod$df.residual, lower.tail = FALSE)
        
        # get formula
        if (!is.null(mul.vars)){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(mul.vars, interaction)
            mul.inter.form <- sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+"))
            coef$formula <- mul.inter.form
          }else{
            mul.form <- sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars, "SampleType"), collapse="+"))
            coef$formula <- mul.form
          }
        } else{
          coef$mvar <- mvar
        }
        
        
        res <- rbind(res, coef)
        
        
        # stop looing for multivariate analysis
        if(!is.null(mul.vars) | !is.null(interaction)){
          break
        }
      }
      
      
      res$padj <- p.adjust(res$pval, method="fdr")
      #res <- res[order(res$time_point),]
      res$comp <- factor(rownames(res), levels=rownames(res))
      res$dir <- ifelse(res$pval < pvalue, ifelse(sign(res$Estimate)==1, "up", "down"), "NS")
      
      print(res)
      
      write.csv(res, quote = FALSE, col.names = NA,file=sprintf("%s/regression_%s.%s.%s.%s%s.csv",out_table,
                                                                project,
                                                                outcome,
                                                                type,
                                                                ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                                                                format(Sys.Date(), "%y%m%d"), sep="/"))
      # return model
      if(!is.null(mul.vars) | !is.null(interaction)){
        saveRDS(mod,sprintf("%s/regression_%s.%s.%s.%s%s.rds",out_table,
                            project,
                            outcome,
                            type,
                            ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                            format(Sys.Date(), "%y%m%d"), sep="/")) 
        
      }
    }
  }
}






