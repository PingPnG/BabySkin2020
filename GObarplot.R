#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/80-bar-plots-and-modern-alternatives/
#https://medium.com/swlh/beautiful-charts-with-r-and-ggpubr-c94122d6b7c6
#author Ping Hu
rm(list=ls())
library(ggplot2)
library(ggpubr)
#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="AdultvsBaby.Up.top.txt"


A<-read.csv(filename, sep="\t", header=TRUE)
dim(A)
A$NegativeLog10Qvalue=-log10(A$qvalue)
jpeg(paste0(filename,".bar.jpg"), width = 1000, height = 1000)
ggbarplot(A, x = "Description", y = "NegativeLog10Qvalue",
          fill = "GOCategory",               # change fill color by cyl
          color = "GOCategory",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,
          group="GOCategory",
          #sorting="asending",
          rotate=TRUE,# Sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
)
dev.off()
