# R source file

#-- to make color table---#
# 190529
# fishtaco에서 색깔을 가져 왔다.
# Arsenic 분석을 예시로 하였고, 작동 된다.
# 200524 function으로 제작 

#-------- give color to phylum --------#

Go_color <- function(cdf, taxaName){
  cPalette <-cdf
  num_of_final_phyla = length(unique(cdf$PhylumCol)); num_of_final_phyla
  cPalette$h = 0
  cPalette$s = 0
  cPalette$v = 0
  num_Actinobacteria = length(grep("Actinobacteria",cPalette$PhylumCol));num_Actinobacteria
  num_Bacteroidetes = length(grep("Bacteroidetes", cPalette$PhylumCol));num_Bacteroidetes
  num_Firmicutes = length(grep("Firmicutes", cPalette$PhylumCol));num_Firmicutes
  num_Proteobacteria = length(grep("Proteobacteria", cPalette$PhylumCol));num_Proteobacteria
  num_Fusobacteria = length(grep("Fusobacteria", cPalette$PhylumCol));num_Fusobacteria
  num_Verrucomicrobia = length(grep("Verrucomicrobia", cPalette$PhylumCol));num_Verrucomicrobia
  num_TM7 = length(grep("TM7", cPalette$PhylumCol));num_TM7
  
  # Synergistetes
  
  number_of_other_phyla = num_of_final_phyla - ((num_Actinobacteria > 0) + (num_Bacteroidetes > 0) + (num_Firmicutes > 0) + (num_Proteobacteria > 0) + (num_Fusobacteria > 0)+ (num_Verrucomicrobia > 0)+ (num_TM7 > 0))
  
  #print(number_of_other_phyla)
  # Actinobacteria = green_pallete
  cPalette[grep("Actinobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0.4, s=seq(0.3,1,length.out=num_Actinobacteria), v=0.9)
  
  # Fusobacteria = orange_pallete
  cPalette[grep("Fusobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0.2, s=seq(0.3,1,length.out=num_Fusobacteria), v=0.9)
  
  # Bacteroidetes = purple_pallete
  cPalette[grep("Bacteroidetes", cPalette$PhylumCol), -1] = expand.grid(h=0.8, s=seq(0.3,1,length.out=num_Bacteroidetes), v=0.9)
  
  # Firmicutes = blue_pallete
  cPalette[grep("Firmicutes", cPalette$PhylumCol), -1] = expand.grid(h=0.6, s=seq(0.3,1,length.out=num_Firmicutes), v=0.9)
  
  # Proteobacteria = red_pallete
  cPalette[grep("Proteobacteria", cPalette$PhylumCol), -1] = expand.grid(h=0, s=seq(0.3,1,length.out=num_Proteobacteria), v=0.9)
  
  # Verrucomicrobia = brown_pallete
  cPalette[grep("Verrucomicrobia", cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_Verrucomicrobia), v=1)
  
  # TM7 = yellow_pallete
  cPalette[grep("TM7", cPalette$PhylumCol), -1] = expand.grid(h=0.1, s=seq(0.3,1,length.out=num_TM7), v=1)
  
  #print(cPalette)
  #print(number_of_other_phyla)
  
  
  print(cPalette)
  
  
  
  # Add other and species name
  cPalette$PhylumCol <- taxaName
  other<-data.frame("[1_#Other]",0,0,0.75)
  names(other)<-c("PhylumCol", "h","s","v")
  color_table <- rbind(other, cPalette)
  
  
  ## hsv to color code ##
  ## taxa위치 변경을 용이 하게 하려면 hsv to color code를 해야 한다. 
  taxa_vs_color = cbind(color_table, apply(color_table[,-1], 1, function(x) hsv(x[1],x[2],x[3])))[,c(1,5)];taxa_vs_color
  colnames(taxa_vs_color) <- c("Taxa", "Color");taxa_vs_color
  class(taxa_vs_color)
  coloring <- as.character(taxa_vs_color$Color) 
  names(coloring) <- taxa_vs_color$Taxa;coloring
  
  
  functionReturningTwoValues <- function() { 
    results <- list()
    results$color_table <- color_table
    results$coloring <-coloring
    return(results) 
  }
  cat("\n")
  print("$color_table and $coloring are returned.")
  
  functionReturningTwoValues()
}
