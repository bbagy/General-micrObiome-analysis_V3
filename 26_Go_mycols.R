
Go_myCols <- function(custumCols=NULL, RColorBrewer=NULL, piratepal=NULL) {
  # reset colors
  # rm(list = ls()[grep("mycols", ls())])
  
  if(!is.null(dev.list())) dev.off()
# https://bookdown.org/ndphillips/YaRrr/more-colors.html

library("yarrr")
  # for custom
  if(is.null(custumCols) & is.null(RColorBrewer) & is.null(piratepal)){
    
    cat("#=== Please select your colors. If not, basic color would be used. ===#","\n","\n", sep=" ")
    
    cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC") # Tableau10
    cols2 <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
    
    cat("Custum Colors1: cols1","\n", cols1, "\n","\n", sep=" ")
    cat("Custum Colors2: cols2","\n", cols2, "\n","\n", sep=" ")

    cat("RColorBrewer:","\n", "Set3     Set2    Set1   Pastel2", "\n", "Pastel1  Paired  Dark2  Accent", "\n","\n", sep=" ")
    cat("piratepal:","\n", "basel    pony    Xmen    decision    southpark", "\n", "google   eternal   evildead    usualsuspects   ohbrother", 
    "\n", "appletv    brave    bugs    cars    nemo", "\n", "rat    up    espresso    ipod   info    info2","\n",sep=" ")

    display.brewer.all(type = "qual")
    readline(prompt="Press [enter] to see other colors.")
    piratepal("all")
  }
  
  if(!is.null(custumCols)){
    if(custumCols == "cols1"){
      cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC") # Tableau10
      mycols <- cols1
      
      barplot(rep(1,length(cols1)), col=cols1, main= custumCols,yaxt="n")
      return(mycols)
    }  else if(custumCols == "cols2"){
      cols2 <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
      mycols <- cols2
      barplot(rep(1,length(cols2)), col=cols2, main= custumCols,yaxt="n")
      return(mycols)
    }
  }  # for RColorBrewer
  else if(!is.null(RColorBrewer)){
    if(RColorBrewer == "Set1"){
      mycols <- brewer.pal(9, RColorBrewer)
    }else if(RColorBrewer == "Set2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Set3"){
      mycols <- brewer.pal(12, RColorBrewer)
    }else if(RColorBrewer == "Pastel2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Pastel1"){
      mycols <- brewer.pal(9, RColorBrewer)
    }else if(RColorBrewer == "Paired"){
      mycols <- brewer.pal(12, RColorBrewer)
    }else if(RColorBrewer == "Dark2"){
      mycols <- brewer.pal(8, RColorBrewer)
    }else if(RColorBrewer == "Accent"){
      mycols <- brewer.pal(8, RColorBrewer)
    }
  barplot(rep(1,length(mycols)), col=mycols, main= RColorBrewer,yaxt="n")
  cat(sprintf("mycols was set as [%s].\n.\n",RColorBrewer))
  return(mycols)
  }   # for piratepal
  else if(!is.null(piratepal)){ 
    if(piratepal == "basel"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "pony"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "Xmen"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "decision"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "southpark"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "google"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal == "eternal"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "evildead"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "usualsuspects"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "ohbrother"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "appletv"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "brave"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "bugs"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "cars"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "nemo"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "rat"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "up"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "espresso"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "ipod"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "info"){
      mycols <- piratepal(palette =piratepal)
    }else if(piratepal== "info2"){
      mycols <- piratepal(palette =piratepal)
    }

 
  
  cols <- as.data.frame(mycols)
  mycols <- cols$mycols
  barplot(rep(1,length(mycols)), col=mycols, main= piratepal,yaxt="n")
  cat(sprintf("mycols was set as [%s].\n.\n",piratepal))
  return(mycols)
  }
}
