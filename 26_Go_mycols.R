
Go_myCols <- function(custumCols, presetCols) {
  # reset colors
  # rm(list = ls()[grep("mycols", ls())])
  
  if(!is.null(dev.list())) dev.off()

  # for custom
  if(is.null(custumCols) & is.null(presetCols)){
    
    cat("#=== Please select your colors. If not, basic color would be used. ===#","\n","\n", sep=" ")
    
    cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC") # Tableau10
    cols2 <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
    
    cat("Custum Colors1: cols1","\n", cols1, "\n","\n", sep=" ")
    cat("Custum Colors1: cols2","\n", cols2, "\n","\n", sep=" ")

    cat("Preset Colors:","\n", "Set3     Set2    Set1   Pastel2", "\n", "Pastel1  Paired  Dark2  Accent", "\n", sep=" ")
    
    display.brewer.all(type = "qual")
    
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
  } 






    
  # for preset colors
  if(is.null(custumCols) & !is.null(presetCols)){
    if(presetCols == "Set1"){
      mycols <- brewer.pal(9, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Set2"){
      mycols <- brewer.pal(8, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Set3"){
      mycols <- brewer.pal(12, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Pastel2"){
      mycols <- brewer.pal(8, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Pastel1"){
      mycols <- brewer.pal(9, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Paired"){
      mycols <- brewer.pal(12, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Dark2"){
      mycols <- brewer.pal(8, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }else if(presetCols == "Accent"){
      mycols <- brewer.pal(8, presetCols)
      barplot(rep(1,length(mycols)), col=mycols, main= presetCols,yaxt="n")
      cat(sprintf("mycols was set as [%s].\n.\n",presetCols))
      return(mycols)
    }
  } else{
    display.brewer.all(type = "qual")
  }
}
