#https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/ComplexHeatmap/inst/doc/s5.legend.html
#http://bioconductor.riken.jp/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html#toc_22
#https://mcr.aacrjournals.org/content/molcanres/suppl/2018/08/09/1541-7786.MCR-18-0114.DC1/196554_3_supp_4905330_pc49zq.pdf
#http://biogps.org/#goto=welcome
#http://biogps.org/downloads/
#https://excellittleknownsecrets.blogspot.com/2017/07/creating-simple-joyplot-in-excel.html
#https://cran.r-project.org/web/packages/ggridges/vignettes/gallery.html
#https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
#https://stackoverflow.com/questions/45384281/ggjoy-facet-with-ggtree
#https://wencke.github.io/ ######GOMAP Fancy
rm(list=ls())

library(devtools)
install_github("jokergoo/ComplexHeatmap")

library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################

#setwd("/Users/ping/Desktop/BabySkin/gene_fig")
filename="lipid_data"


A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);

B=log10(A[1:d[1], 5:d[2]])
rownames(B)=paste(A[,2],A[,1])
#B$syn=A[,2]
Cname=colnames(A)[5:d[2]]
#RM<-rowMeans(B)
#x=data.frame(cbind(rowMeans(B),A[,2],A[,1]))
#colnames(x)=c("RM","Syn","Probe")
#library(dplyr)
#x1 <- x %>% group_by(Syn) %>% filter(RM == max(RM)) 
#B$syn=A[,2]
#C<-B[rownames(B) %in% x1$Probe]

Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)
ID=rep(NA,Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    ID[mm]=splitname[[mm]][4]
}
#https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#png(filename=paste0(filename,".png"), height=2800, width=3600, res=300)
#install.packages("pheatmap")
library(pheatmap)
#pheatmap(B)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(B, 1, cal_z_score))
#pheatmap(data_subset_norm)

my_hclust_gene <- hclust(dist(data_subset_norm), method = "complete")
#install.packages("dendextend")
library(dendextend)

#as.dendrogram(my_hclust_gene) %>%
#  plot(horiz = TRUE)
#dev.off()
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
#set.seed(1984)
#my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
#my_gene_col$random <- my_random

#head(my_gene_col)
my_sample_col =data.frame(cbind(Group1,Group2, Group3))
# change order
#my_sample_col$sample <- factor(my_sample_col$Group1, levels = c("Baby", "Adult"))
#change color
#my_colour = list(
#  sample = c(normal = "#5977ff", tumour = "#f74747"),
#  random = c(random1 = "#82ed82", random2 = "#9e82ed"),
#  cluster = c(cluster1 = "#e89829", cluster2 = "#cc4ee0")
#)

row.names(my_sample_col) <- colnames(B)
png(filename=paste0(filename,".png"), height=4800, width=2800, res=120)

p1<-pheatmap(data_subset_norm,
         #annotation_colors = my_colour,
         annotation_row = my_gene_col,
         annotation_col = my_sample_col,
         cutree_rows = 2,
         cutree_cols = 2)
print(p1)
dev.off()

####retrieve result from the heatmap
# use silent = TRUE to suppress the plot
#my_heatmap <- pheatmap(data_subset_norm, silent = TRUE)
# results are stored as a list
#class(my_heatmap)
#[1] "list"
#names(my_heatmap)
#[1] "tree_row" "tree_col" "kmeans"   "gtable"
#my_heatmap$tree_row %>%
#  as.dendrogram() %>%
#  plot(horiz = TRUE)

#save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
#  png(filename, width = width, height = height, res = res)
#  grid::grid.newpage()
#  grid::grid.draw(x$gtable)
#  dev.off()
#}

#save_pheatmap_png(my_heatmap, "my_heatmap.png")