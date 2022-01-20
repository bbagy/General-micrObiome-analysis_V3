##############################################
#---- install package  by bioconductor ------#
##############################################
# 210322
# install package and reads library is combined

# version 1
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioconductors <- c("dada2","DESeq2", "dplyr","ggpubr","ggfortify", "ggpmisc",
                   "illuminaio","msa","phyloseq","rstatix","useful","DECIPHER")

for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor)
  }else{library(bioconductor, character.only = TRUE)}
}

##############################################
#----          install package         ------#
##############################################
packages <- c("ape", "car","cluster","CLME","cowplot","crayon", "caret","colorspace","e1071",
           "digest","data.table", "devtools","doParallel","ellipse", "emmeans",
           "gplots","ggplot2","grid","gridExtra","gplots","ggrepel",
           "Hmisc","huge","irlba","igraph","irr","lme4","lmerTest",
           "Matrix","magrittr","MASS","missForest","nlme","phangorn","plot3D",
           "pheatmap","pkgconfig","plyr","parallel","pscl","plotly","rfUtilities",
           "rlang","randomForest","readxl","RColorBrewer","ROCR","reshape","reshape2",
           "stringi","S4Vectors","ShortRead","tidyverse","vegan","VGAM") #"venneuler",
# version 1
#for (pack in packs){install.packages(sprintf("%s",pack))}
# version 2 (better version)
for (package in packages){
  if(!package %in% installed.packages()){
    install.packages(package)
  }else{library(package, character.only = TRUE)}
}

#for (package in packages){library(sprintf("%s",package), character.only = TRUE)}

