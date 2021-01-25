rm(list=ls())
#options(echo=TRUE) # if you want see commands in output file
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename ="L6RAWhite_lipid_data.txt"
#outname ="test"
#filename=args[1]
#outname=args[2]
#rm(args)

A<-read.table(filename, sep="\t", header=TRUE)

d <- dim(A);

test_data <- data.matrix(A[,2:(d[2]-2)])
rownames(test_data) <-A[,d[2]]
#####only select average more than 0.005's taxon for correlation analysis####
#Select= (rowMeans(test_data) >=0.005)
#new_data=test_data[Select,]
xx<-t(test_data)


#colnames(xx) =A[Select,1]
#rownames(xx)=colnames(A)[2:d[2]]

Cname=rownames(xx)
Clen=length(Cname) ##there are 3 annotation columns
Site=rep("NA", Clen)
AgeGrp=rep("NA", Clen)
Age=rep("NA", Clen)
SID=rep("NA", Clen)
Order=rep("NA", Clen)
SiteAgeGrp=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
       Site[mm]=splitname[[mm]][1]
       AgeGrp[mm]=splitname[[mm]][3]
       Age[mm]=splitname[[mm]][4]
       SID[mm]=splitname[[mm]][5]
       Order[mm]=splitname[[mm]][6]
       SiteAgeGrp[mm]=paste0(Site[mm], AgeGrp[mm])
}


library("Hmisc")
#install.packages("corrplot")
library(corrplot)
mydata.cor = cor(as.matrix(xx[1:5,1:5]), method = c("spearman"))
#rCOR=rcorr(as.matrix(xx))
jpg(filename="lipid-all-cor", width=3800, height=3800,res=120)
p1<-corrplot(mydata.cor)
#print(p1)
dev.off()
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ff<-flattenCorrMatrix(rCOR$r, rCOR$P)
write.table(ff,file=paste0(outname,'.lipid.corFlat.xls'),sep="\t")
for(yyy in unique(Site)){
   data_i=xx[Site == yyy,]
   rCOR=rcorr(as.matrix(data_i))
   ff<-flattenCorrMatrix(rCOR$r, rCOR$P)
   write.table(ff,file=paste0(outname,yyy,'.lipid.corFlat.xls'),sep="\t")
}