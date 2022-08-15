
Go_regression <- function(data, project,  
                          outcomes, 
                          cate.vars=NULL, 
                          con.vars=NULL,
                          mul.vars=F,
                          interaction=NULL, 
                          randomEff=NULL,
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
  if(!is.null(cate.vars)){
    for (cate in  cate.vars) {
      data[,cate] <- factor(data[,cate], levels = intersect(orders, data[,cate]))
    }
  } 
  if(!is.null(con.vars)){
    for (con in  con.vars) {
      # NA 제거
      data[,con] <- as.character(data[[con]]);data[,con]
      data[,con][data[,con]==""] <- "NA";data[,con]
      data[,con] <- as.numeric(data[[con]])
    }
  }


  


  if(!is.null(cate.vars) & !is.null(con.vars)){
    varis <- unique(c(cate.vars, con.vars))
  }

  
  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes){
    if (class(data[,outcome]) == "character") {
      if (class(data[,outcome]) == "character" | length(unique(data[,outcome])) > 3){
        print("Not able to analysis using")
        break
      }
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
    } else if (class(data[,outcome])  == "numeric") {
      
      # NA 제거
      data[,outcome] <- as.character(data[[outcome]]);data[,outcome]
      data[,outcome][data[,outcome]==""] <- "NA";data[,outcome]
      # data <- subset(data, data[,mvar] != "NA");data[,outcome]  # subset 를 사용한 NA 삭제
      data[,outcome] <- as.numeric(as.character(data[[outcome]]))
      
      data[,outcome] <- as.numeric(data[,outcome])
    } 
    
    res <- {}
    
    for (mvar in varis) {
      if (outcome == mvar) {
        next
      }
      
      # print(sprintf("##-- %s (total without NA: %s/%s) --##", 
      #              mvar, dim(data.na)[1], dim(data)[1]))
      
      if (length(unique(data[,mvar])) ==1) {
        next
      }
      
      # get formula
      if(!is.null(randomEff)){
        if(mul.vars == FALSE){
          form <- as.formula(sprintf("%s ~ (1 | %s) + %s", outcome, randomEff,mvar))
          print("Univariate anaysis")
          type <- "LMEM_uni"
        } else if (mul.vars == TRUE){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(varis, interaction)
            form <- as.formula(sprintf("%s ~ (1 | %s) + %s", outcome, randomEff,paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+")))
            print("Multivariate anaysis with interaction effect")
            type <- "LMEM_multi-interantion"
          }else{
            form <- as.formula(sprintf("%s ~ (1 | %s) + %s", outcome, randomEff,paste(setdiff(varis, "SampleType"), collapse="+"))) # , "SampleType"
            print("Multivariate anaysis")
            type <- "LMEM_multi"
          }
        } 
      }else{
        if(mul.vars == FALSE){
          form <- as.formula(sprintf("%s ~ %s", outcome, mvar))
          print("Univariate anaysis")
          type <- "uni"
        } else if (mul.vars == TRUE){
          if (!is.null(interaction)){
            mul.vars.interaction <- c(varis, interaction)
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+")))
            print("Multivariate anaysis with interaction effect")
            type <- "multi-interantion"
          }else{
            form <- as.formula(sprintf("%s ~ %s", outcome, paste(setdiff(varis, "SampleType"), collapse="+"))) # , "SampleType"
            print("Multivariate anaysis")
            type <- "multi"
          }
        } 
      }
      

      
      print(form)
      
      #=============#
      #    model    #
      #=============#
      
      if(!is.null(randomEff)){
        m <- "Regression (LMEM)"
        mod <- lmer(form, data=data, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      }else{
        if (class(data[,outcome]) == "numeric"){
          m <- "Regression (glm-poisson)"
          mod <- glm(form, data=data,  family = poisson(link='log'))
        } else if (length(unique(data[,outcome])) == 2){
          m <- "Logistic regression (glm-binomial)"
          mod <- glm(form, data=data,  family=binomial("logit"))
        } 
      }

      print(m)
      
      
      
      # out for the model
      coef <- as.data.frame(summary(mod)$coefficients)
      coef <- coef[setdiff(rownames(coef), "(Intercept)"),,drop=F]
      colnames(coef) <- c("Estimate", "SE", "t", "pval")
      
      if(!is.null(randomEff)){
        colnames(coef) <- c("Estimate", "SE", "df","t", "pval")
      }else{
        colnames(coef) <- c("Estimate", "SE", "t", "pval")
      }
      
      
      
      if (dim(coef)[1] == 0){
        next
      }
      
      
      
      # out for the confidence interval 
      conf <- data.frame(confint(mod))
      conf <- conf[setdiff(rownames(conf), "(Intercept)"),,drop=F]
      conf.na <- na.omit(conf) 
      colnames(conf.na) <- c("2.5 %", "97.5 %")
      
      if(!is.null(randomEff)){
        coef <- coef
      }else{
        coef$`2.5 %` <- conf.na$`2.5 %`
        coef$`97.5 %` <- conf.na$`97.5 %`
      }

      
      coef$outcome <- outcome
      coef$mvar <- mvar
      coef$model <- m
      
      if(!is.null(randomEff)){
        coef <- coef
      }else{
        coef$deviance <- pchisq(q=mod$null.deviance-mod$deviance,df=mod$df.null-mod$df.residual, lower.tail = FALSE)
      }

      
      # get formula

    
      if (mul.vars == TRUE){
        if (!is.null(interaction)){
          mul.vars.interaction <- c(varis, interaction)
          mul.inter.form <- sprintf("%s ~ %s %s", outcome, ifelse(is.null(randomEff), "", sprintf("(1 | %s) +", randomEff)), 
                                    paste(setdiff(mul.vars.interaction, "SampleType"), collapse="+"))
          coef$formula <- mul.inter.form
        }else{
          mul.form <- sprintf("%s ~ %s %s", outcome, ifelse(is.null(randomEff), "", sprintf("(1 | %s) +", randomEff)),
                              paste(setdiff(varis, "SampleType"), collapse="+"))
          coef$formula <- mul.form
        }
      } else{
        coef$mvar <- mvar
      }
      
      
      
      
      
      
      ifelse(!is.null(randomEff), "", sprintf("(1 | %s)", randomEff))
      
      res <- rbind(res, coef)
      
      
      # stop looing for multivariate analysis
      if(mul.vars == TRUE | !is.null(interaction)){
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
    if(mul.vars == TRUE | !is.null(interaction)){
      saveRDS(mod,sprintf("%s/regression_%s.%s.%s.%s%s.rds",out_table,
                          project,
                          outcome,
                          type,
                          ifelse(is.null(name), "", paste(name, ".", sep = "")),  
                          format(Sys.Date(), "%y%m%d"), sep="/")) 
      
    }
  }
}

