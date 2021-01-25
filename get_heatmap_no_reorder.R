#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#filename <- args[1]
#outname <- args[2]
#rm(args)
########################################
#draw boxplot, first line format: name.position.color...
#probe\tgenesymbol\tgenename\tdata
#https://www.biostars.org/p/18211/
#https://www.youtube.com/watch?v=561MxCXX-DE
##########################################
rm(list = ls())
filename="PT2_PT5_FC1.2_data2.txt"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
outname="PT2_PT5_FC1.2"
####cleanup column information
Cname=colnames(A)[4:17]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
#strain=rep("probe", Clen)
#boxOrder=rep(0, Clen)
Ttype=rep("NA", Clen)
ID=rep("NA", Clen)
tissue=rep("NA",Clen)
for(mm in  1:Clen ){
      # strain[mm]=splitname[[mm]][1]
      # boxOrder[mm]=as.numeric(splitname[[mm]][2])
       Ttype[mm]=splitname[[mm]][3]
       ID[mm]=splitname[[mm]][7]
       tissue[mm]=splitname[[mm]][5]
}
#C<-unique(cbind(strain, boxOrder, boxColor))
#Ulen=length(C[,1])
########################



library("RColorBrewer")
library("amap")
library("pvclust")
library("gplots")

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
test_data <- as.matrix(A[,4:17])
rownames(test_data)=paste0(A[,2],".", A[,1])
#class(test_data)
#sampleType=data.frame(Ttype)
#my_data=t(test_data)
#fontsize <-0.6

scal_rows <-function (x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
my_data =scal_rows(test_data)

Heatmap(my_data, 
        cluster_columns = FALSE,
        clustering_distance_rows = "maximum",
        column_names_gp =  gpar(fontsize=6),
        clustering_method_rows = "complete" , #"ward.D"
        row_names_gp = gpar(fontsize=0.6)
        #show_row_names = TRUE,
        #row_dend_side = "left",
        #row_names_side = "left"
)
library(dendextend)
dend=hclust(dist(my_data, method="maximum"), method="ward.D")
Cdend=hclust(dist(t(my_data), method="maximum"), method="ward.D")

Heatmap(my_data, 
        cluster_columns = FALSE,
        #clustering_distance_rows = "maximum",
        column_names_gp =  gpar(fontsize=6),
        #clustering_method_rows = "complete" , #"ward.D"
        row_names_side = "left",
        row_names_gp = gpar(fontsize=0.6),
        cluster_rows = color_branches(dend, k=3),
        #cluster_columns = color_branches(Cdend, k=2)
)
TTYPE=data.frame(Ttype)
Tcolor=c(rep("red",7), rep("blue", 7))
names(Tcolor)<-Ttype
png(filename=paste0(outname,".heat.png"), height=4800, width=1600, res=300)


ha=HeatmapAnnotation(Ttype, col = list(Ttype = Tcolor))
Heatmap(my_data, 
        cluster_columns = FALSE,
        clustering_distance_rows = "maximum",
        column_names_gp =  gpar(fontsize=8),
        clustering_method_rows = "complete" , #"ward.D"
        row_names_side = "left",
        row_names_gp = gpar(fontsize=2),
        #cluster_rows = color_branches(dend, k=3),
        #cluster_columns = color_branches(Cdend, k=2)
        km=5, 
        bottom_annotation = HeatmapAnnotation(TTYPE, col=Tcolor, show_legend = FALSE)
)
dev.off()

##Gene correlation graph
#c<-cor(t(test_data))
#h<-Heatmap(c, show_column_names = FALSE, show_row_names = FALSE)
#print(h)

C<-cor(my_data)
h<-Heatmap(C, show_column_names = TRUE, show_row_names = FALSE)
png(filename=paste0(outname,".sampleCor.png"), height=800, width=800, res=120)
print(h)
dev.off()

#h<-hclust(d=dist(C))
#hc <- as.factor(cutree(h,k=4))
#hh=ComplexHeatmap::Heatmap(C[h$order, h$order],cluster_rows=FALSE, cluster_columns = FALSE)
#print(hh)
###another way for transformation
#Heatmap(log10(1+test_data))
####different color## can choose colorblind safe
#library(RColorBrewer)
#Heatmap(C, col=brewer.pal(7,"BuGn"))##can use Greys for black and white
#display.brewer.all(colorblindFriendly = TRUE)

#colnames(test_data)=Ttype
Tcolor=c(rep("red",7), rep("blue", 7))
row_name=paste0(A[,2],".",A[,1])	
png(filename=paste0(outname,".heat.png"), height=1600, width=1800, res=120)
heatmap.2(test_data, Colv='FALSE',Rowv='FALSE',
labCol=Ttype,
labRow=row_name,scale='row',cexRow=0.5,
col=redgreen(75), dendrogram='none',main=paste0(outname, ".heatmap"),
trace='none', density='none', ColSideColors=Tcolor);
dev.off()
#png(filename=paste0(outname,".notnormalized.png"), height=1600, width=1800, res=120)
#heatmap.2(test_data, trace='none', density='none')
#dev.off()