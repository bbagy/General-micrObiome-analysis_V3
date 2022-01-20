pairwise.adonis2 <- function(x,factors, map,mvar, sim.function = 'vegdist', sim.method = 'bray', adjust, p.adjust.m ='bonferroni'){
  set.seed(1)
  
  co <- combn(unique(as.character(map[,mvar])),2)
  for(elem in 1:ncol(co)){
    x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                    factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    
    # run
    map.pair <- subset(map, map[,mvar] %in% c(co[1,elem],co[2,elem]))
    
    if (count(map.pair[,mvar])[1,2] <=2 | count(map.pair[,mvar])[2,2] <=2){
      next
    }
    
    if (!is.null(adjust)) {
      form <- as.formula(sprintf("x1 ~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
      print(form)
    }else{
      form <- as.formula(sprintf("x1 ~ %s", mvar))
      print(form)
    }
    
    ad <- adonis2(form, data = map.pair, permutations=999, by="margin")# "terms"  "margin" NULL
    
    pairs <- c()
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c()
    Df <- c(Df,ad[1,1])
    SumsOfSqs <- c()
    SumsOfSqs <- c(SumsOfSqs, ad[1,2])
    R2 <- c()
    R2 <- c(R2,ad[1,3])

    F.Model <- c()
    F.Model <- c(F.Model,ad[1,4]);
    p.value <- c()
    p.value <- c(p.value,ad[1,5])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  sig[p.adjusted <= 0.0001] <-'****'
  
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,R2,F.Model,p.value,p.adjusted,sig)
  
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}

