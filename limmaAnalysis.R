#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
#GSS<-args[2]
#############################################################################
#author: Ping Hu
#Date:9-16-2014
#Function: GSS1423, GSS1424, GSS1429 
#calculate the stats agains all S_aureus(pathogen)
#signal data column name first part is the strain or chemical name
#############################################################################
filename="GSS20.Mom_v_Daughter.signal"
library("Biobase")
library(genefilter)
library(limma)

options(stringsAsFactors=F)
#A=read.csv(filename,header=TRUE,as.is=TRUE)
A<-read.table(filename, sep="\t", header=TRUE)
probe=A[,1]
####This is my control columns
MomCol=grep("Mother", as.character(colnames(A)))
DCol=grep("Daughter", as.character(colnames(A)))

Mmean=rowMeans(A[,c(MomCol)])
Dmean=rowMeans(A[,c(DCol)])
###calculate fc##########
fc=Mmean/Dmean
B<-log2(A[,c(MomCol, DCol)])
####regular two sample ttest#####
B_type<-c(rep("Mom",length(MomCol)),rep("Daughter",length(DCol)))
Pair=c(c(1:9), c(1:9))

tt=rowttests(data.matrix(B), factor(B_type))
######limma moderated ttest(they call anova######


design  = cbind (Mom = as.numeric(B_type == "Mom"), child=as.numeric(B_type == "Daughter") )


contrast=makeContrasts(Mom - child, levels = design)

fit=lmFit(B, design)
fit2= contrasts.fit(fit, contrast)

fit2 = eBayes(fit2)                                       
x=cbind(probe, tt$p.value, fit2$p.value, fc)

#makeContrasts(Tumor - Normal, levels = design)
#fit <- contrasts.fit(fit, contrast)


colnames(x)[2]="pval.two_sample_ttest.MvD"
colnames(x)[3]="pval.limma_moderated_ttest.MvD"
colnames(x)[4]="fold.MvD"
write.csv(x, file="GSS20.MvD.limma", row.names=FALSE)

library(statmod)
corfit = duplicateCorrelation(B, design, block=Pair)
fit = lmFit(B, design=design, block=Pair, correlation=corfit$consensus.correlation)
contrast=makeContrasts(Mom - child, levels = design)
fit2 = contrasts.fit(fit, contrast)
fit2 = eBayes(fit2)
x=cbind(probe, tt$p.value, fit2$p.value, fc)
colnames(x)[2]="pval.two_sample_ttest.MvD"
colnames(x)[3]="pval.limma_paired_ttest.MvD"
colnames(x)[4]="fold.MvD"
write.csv(x, file="GSS20.MvD.pairedlimma", row.names=FALSE)
