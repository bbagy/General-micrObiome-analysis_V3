##############################################
#----          install package         ------#
##############################################

# for ignore yes or no
options(needs.promptUser = FALSE) 
packages <- c("ape", "Boruta","car","cluster","CLME","compositions","cowplot","crayon", "caret","colorspace",
           "digest","data.table", "devtools","doParallel","ellipse", "emmeans","e1071",
           "gplots","ggplot2","grid","gridExtra","gplots","ggrepel","doRNG","ggalluvial","ggforce",
           "Hmisc","huge","irlba","igraph","irr","lme4","lmerTest","nnet","MLmetrics",
           "Matrix","magrittr","MASS","missForest","nlme","phangorn",#"plot3D",
           "pheatmap","pkgconfig","plyr","parallel","pscl","plotly","rfUtilities",
           "rlang","randomForest","readxl","RColorBrewer","ROCR","reshape","reshape2","yarrr",
           "stringi","S4Vectors","tidyverse","vegan","VGAM","ShortRead","picante") #"venneuler",
# version 1
#for (pack in packs){install.packages(sprintf("%s",pack))}
# version 2 (better version)
for (package in packages){
  if(!package %in% installed.packages()){
    install.packages(package)
  }else{library(package, character.only = TRUE)}
}

#for (package in packages){library(sprintf("%s",package), character.only = TRUE)}



##############################################
#---- install package  by bioconductor ------#
##############################################
# 210322
# install package and reads library is combined

# version 1
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")



bioconductors <- c("phyloseq","microbiome","Rhtslib","dada2", "dplyr","ggpubr","ggfortify","genefilter", "ggpmisc", "DESeq2", "ANCOMBC", # "msa",
                   "illuminaio","rstatix","useful","DECIPHER","ComplexHeatmap")

# 


for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor, force = TRUE)
  }else{library(bioconductor, character.only = TRUE)}
}


##############################################
#----            Introduction          ------#
##############################################

cat(blue("#--------------------------------------------------------------# \n"))
cat(blue("#------       General analysis Of microbiome (Go)        ------# \n"))
cat(blue("#------    Quick statistics and visualization tools      ------# \n"))
cat(blue("#--------------------------------------------------------------# \n"))
cat(red("                                      Version: Go_tools.4.0.6 \n"))
cat("                                              Write by Heekuk \n")
cat(yellow("All the required packages were installed.\n"))
cat(yellow("All the required packages were loaded.\n"))
cat(blue("#--------------------------------------------------------------# \n"))

