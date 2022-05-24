
Go_regression <- function(df, project,  
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
      if (class(df$mvar) == "character") {
        df[,mvar] <- factor(df[,mvar])
        # setting baseline
        df[,mvar] <- factor(df[,mvar], levels = intersect(orders, df[,mvar]))
        # NA 제거
        df[,mvar] <- as.character(df[[mvar]]);df[,mvar]
        df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
        #df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
        
      } else if (class(df$mvar) == "numeric" | min(mvar) < 0) {
        # NA 제거
        df[,mvar] <- as.character(df[[mvar]]);df[,mvar]
        df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
        #df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
        df[,mvar] <- as.numeric(as.character(df[[mvar]]))
      } 
    }
  }



  #----------------------------------------------------#
  #--------------    regression model     -------------#
  #----------------------------------------------------#
  set.seed(1)
  for (outcome in outcomes){
    if (class(outcome) == "character") {
      df[,outcome] <- factor(df[,outcome])
      df[,outcome] <- factor(df[,outcome], levels = intersect(orders, df[,outcome]))
      
      # NA 제거
      df[,outcome] <- as.character(df[[outcome]]);df[,outcome]
      df[,outcome][df[,outcome]==""] <- "NA";df[,outcome]
      #df.na <- subset(df, df[,outcome] != "NA");df.na[,outcome]  # subset 를 사용한 NA 삭제
      # set the baseline for outcome

      if(length(unique(df[,outcome])) == 2){
        
        df[,outcome] <- factor(df[,outcome])
        out <- levels(df[,outcome])[1]
        

        df[,outcome] <- factor(ifelse(df[,outcome]== levels(df[,outcome])[1],0,1), levels=c(0,1), 
                               labels = levels(df[,outcome]))
        
        
        print(levels(df[,outcome]))
      }
    } else if (class(outcome)  == "numeric") {
      
      # NA 제거
      df[,outcome] <- as.character(df[[outcome]]);df[,outcome]
      df[,outcome][df[,outcome]==""] <- "NA";df[,outcome]
      # df <- subset(df, df[,mvar] != "NA");df[,outcome]  # subset 를 사용한 NA 삭제
      df[,outcome] <- as.numeric(as.character(df[[outcome]]))
      
      df[,outcome] <- as.numeric(df[,outcome])
    } 
    
    res <- {}
    
    for (mvar in uni.vars) {
      if (outcome == mvar) {
        next
      }
      
     # print(sprintf("##-- %s (total without NA: %s/%s) --##", 
      #              mvar, dim(df.na)[1], dim(df)[1]))
      
      if (length(unique(df[,mvar])) ==1) {
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
      if (class(df[,outcome]) == "numeric"){
        m <- "Regression (glm-poisson)"
        mod <- glm(form, data=df,  family = poisson(link='log'))
      } else if (length(unique(df[,outcome])) == 2){
        m <- "Logistic regression (glm-binomial)"
        mod <- glm(form, data=df,  family=binomial("logit"))
      } else if (length(unique(df[,outcome])) > 3){
        m <- "Multinomial logistic regression"
        mod <- multinom(form, data=df)
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
      return(mod)
    }
  }
}






