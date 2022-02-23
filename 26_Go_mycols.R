
Go_mycols <- function(custumCols, presetCols) {
  # reset colors
  # rm(list = ls()[grep("mycols", ls())])
  
  # for custom
  if(is.null(custumCols) & is.null(presetCols)){
    
    cat("#=== Please select your colors. If not, basic color would be used. ===#","\n","\n", sep=" ")
    
    cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BABOAC") # Tableau10
    cat("Custum Colors1: cols1","\n", cols1, "\n","\n", sep=" ")
    cat("Preset Colors:","\n", "Set3     Set2    Set1   Pastel2", "\n", "Pastel1  Paired  Dark2  Accent", "\n", sep=" ")
    
    display.brewer.all(type = "qual")
    
  }
  
  if(!is.null(custumCols)){
    if(custumCols == "cols1"){
      cols1 <- c("#1170aa", "#fc7d0b",  "#76B7B2", "#B07AA1", "#E15759","#59A14F","#EDC948", "#FF9DA7", "#9C755F", "#BABOAC") # Tableau10
      mycols <- cols1
      return(mycols)
    }
  }

    
  # for preset colors
  if(is.null(custumCols) & !is.null(presetCols)){
    if(presetCols == "Set1"){
      mycols <- brewer.pal(9, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Set2"){
      mycols <- brewer.pal(8, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Set3"){
      mycols <- brewer.pal(12, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Pastel2"){
      mycols <- brewer.pal(8, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Pastel1"){
      mycols <- brewer.pal(9, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Paired"){
      mycols <- brewer.pal(12, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Dark2"){
      mycols <- brewer.pal(8, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }else if(presetCols == "Accent"){
      mycols <- brewer.pal(8, presetCols)
      cat(sprintf("mycols was set as [%s].\n",presetCols))
      return(mycols)
    }
  } else{
    display.brewer.all(type = "qual")
  }
}
