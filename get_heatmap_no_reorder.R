args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
outname <- args[2]
rm(args)
########################################
#draw boxplot, first line format: name.position.color...
#probe\tgenesymbol\tgenename\tdata
#
##########################################

#filename="GSS1429.collagen.data"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)

####cleanup column information
Cname=colnames(A)[4:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
strain=rep("probe", Clen)
boxOrder=rep(0, Clen)
boxColor=rep("yellow", Clen)
for(mm in  1:Clen ){
       strain[mm]=splitname[[mm]][1]
       boxOrder[mm]=as.numeric(splitname[[mm]][2])
       boxColor[mm]=splitname[[mm]][3]
       
}
C<-unique(cbind(strain, boxOrder, boxColor))
Ulen=length(C[,1])
########################

test_data <- data.matrix(A[1:d[1],4:d[2]]);

library("RColorBrewer")
library("amap")
library("pvclust")
library("gplots")


colnames(test_data)=strain

row_name=A[,2]	
png(filename=paste0(outname,".heat.png"), height=2800, width=3600, res=300)
heatmap.2(as.matrix(test_data), Colv='FALSE',Rowv='FALSE',
labCol=strain,labRow=row_name,scale='row',cexRow=0.5,
col=redgreen(75), dendrogram='none',main=paste0(outname, ".heatmap"),
trace='none', density='none', margins=c(20,12),ColSideColors=boxColor);

dev.off()

