#' A Go_linear
#'
#' This function help to generation of linear regression plots
#' @param test test
#' @keywords basic statistic and simple plots
#' @export
#' @examples

Go_linear <- function(df, cont.vars, project, outcomes, method="lm",
                      mycols =NULL, maingroup=NULL, orders=NULL, name=NULL, 
                      height, width, plotCols, plotRows){
    
  if(!is.null(dev.list())) dev.off()
    
print("Method option: glm, lm, loess, gam")
  # out dir
  out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out)) dir.create(out)
  out_path <- file.path(sprintf("%s_%s/pdf",project, format(Sys.Date(), "%y%m%d"))) 
  if(!file_test("-d", out_path)) dir.create(out_path)

  # out file
  # "name" definition
  if (class(name) == "function"){
    name <- NULL
  }
  
  tt <- try(mycols,T)
  if(class(tt) == "try-error"){
    print("mycols is not defined.")
    mycols <- NULL
  }
  tt <- try(orders,T)
  if(class(tt) == "try-error"){
    print("orders is not defined.")
    orders <- NULL
  }
  pdf(sprintf("%s/linear.%s.%s.%s%s%s.pdf", out_path, 
              project, 
              method,
              ifelse(is.null(maingroup), "", paste(maingroup, ".", sep = "")), 
              ifelse(is.null(name), "", paste(name, ".", sep = "")), 
              format(Sys.Date(), "%y%m%d")), height = height, width = width)
  
  

  my.formula <- y ~ x

  
  # plot
  plotlist <- list()
  for (mvar in cont.vars) {
    if (!is.null(maingroup)){
      if (mvar == maingroup){
        next
      }
    }
    # Na 제거
    df[,mvar] <- as.numeric(as.character(df[[mvar]]))
    #df[,mvar] <- as.numeric(as.character(df$percenttwl))
    df[,mvar][df[,mvar]==""] <- "NA";df[,mvar]
    df[,mvar]<- as.numeric(df[,mvar]);df[,mvar]
    
    df.na <- subset(df, df[,mvar] != "NA");df.na[,mvar]  # subset 를 사용한 NA 삭제
    
    print(sprintf("##-- %s (total without NA: %s/%s) --##", 
                  mvar, dim(df.na)[1], dim(df)[1]))
    if (length(unique(df.na[,mvar])) ==1) {
      next
    }
    summary.df.na <- summary(df.na[,mvar])
    
    # na제거 in the maingroup
    if (!is.null(maingroup)) {
      df.na[,maingroup] <- as.character(df.na[,maingroup]);df.na[,maingroup]
      df.na[,maingroup][df.na[,maingroup]==""] <- "NA";df.na[,maingroup]
      df.na[,maingroup]<- as.factor(df.na[,maingroup]);df.na[,maingroup]
      df.na.na <- subset(df.na, df.na[,maingroup] != "NA");df.na.na[,maingroup]
      
      if (!is.null(orders)) {
        df.na.na[,maingroup] <- factor(df.na.na[,maingroup], levels = orders)
      }else{
        df.na.na[,maingroup] <- factor(df.na.na[,maingroup])
      }

    }

    
    
    
    for(i in 1:length(outcomes)){
      
      if (outcomes[i] == mvar | outcomes[i] == "Chao1" & mvar == "Shannon" | outcomes[i] == "Shannon" & mvar == "Chao1") {
        print(sprintf("Stop function bacause out was %s and mvar was %s", outcomes[i], mvar))
        next
      }
      
      df.na[,outcomes[i]] <- as.numeric(df.na[,outcomes[i]])
      
      print(outcomes[i])
      
      if (!is.null(maingroup)) {
        df.na.na[,outcomes[i]] <- as.numeric(df.na.na[,outcomes[i]])
        p<- ggplot(df.na.na, aes_string(x=mvar, y=outcomes[i], group= maingroup, color=maingroup, linetype = maingroup))
          
      }else {
        p <- ggplot(df.na, aes_string(x=mvar, y=outcomes[i]))+theme_classic() + geom_point(size = 0.5)
        
      }
      
     p <- p + theme_classic() + geom_point(size = 0.5)
       # scale_colour_brewer(palette = colorset) + 

       if(method == "glm"){
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5, method.args = list(family = "poisson"))   #method.args = list(family = "poisson")
       }else if(method == "lm"){
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5) + 
           #stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE, size = 3) +
           stat_fit_glance(method.args = list(formula = my.formula), method = method, 
                           #geom = 'text', 공식이 한쪽으로 정리가 되지 않고, 라인에 수치가 붙는다.
                           aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g', 
                                               stat(r.squared), stat(p.value))), 
                           parse = TRUE, size = 3)
       }else{
         p <- p + geom_smooth(method = method, formula = my.formula, linetype="solid", fill="lightgrey", se=T, size=0.5)
       }
     
     p <- p + ggtitle(sprintf("%s with %s", mvar, outcomes[i])) + theme(title=element_text(size=10)) + labs(x = NULL)+
       theme(title=element_text(size=10),
             axis.text.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = .5, face = "plain"),
             axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))
     
      
      if(!is.null(mycols)){
        p <- p + scale_color_manual(values = mycols)
      }else{
        p <- p
      }



      plotlist[[length(plotlist)+1]] <- p
    }
  }
  multiplot(plotlist=plotlist, cols=plotCols, rows=plotRows)
  dev.off()
}
