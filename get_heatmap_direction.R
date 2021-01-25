args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]

rm(args)
########################################
#draw boxplot, first line format: name.position.color...
#probe\tgenesymbol\tgenename\tdata
#
##########################################

#filename="GSS1429.collagen.data"
#filename="AMP.GSS1412.xls"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)

####cleanup column information
Cname=colnames(A)[6:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")

Treat=rep(0, Clen)
for(mm in  1:Clen ){
      
       Treat[mm]=splitname[[mm]][3]
       
}


test_data <- data.matrix(A[1:d[1],6:d[2]]);
d2=dim(test_data)
#change_data <-matrix(0,nrow=d2[1], ncol=d2[2])
#change_data[test_data<=0.1 &test_data>0]=1
#change_data[test_data>=-0.1 &test_data<0]=-1


library("RColorBrewer")
library("amap")
library("pvclust")
library("gplots")


colnames(test_data)=Treat
rownames(test_data)=A[,2]
rowname=paste(A[,2],as.integer(A[,4]), sep=".")	
png(filename=paste0(filename,".png"), height=2800, width=3600, res=300)
heatmap.2(as.matrix(test_data), Colv='FALSE',Rowv='FALSE',
labCol=Treat,labRow=rowname,scale='none',
#cellnote=round(test_data,2),
col=greenred(16), dendrogram='none',main=filename,
trace='none', density='none', margins=c(20,12));
dev.off()

