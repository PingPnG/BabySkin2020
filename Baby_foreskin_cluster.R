rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
##########################
#Ping Hu
#1-25-2021
###################################
filename="BabyOnly_No_outlier.txt"


A<-read.table(filename, sep="\t", header=TRUE)
d<-dim(A)
B=log10(A[1:d[1], 4:d[2]])
rownames(B)= paste0(A[,1], ".",A[,2])
Cname=colnames(A)[4:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)

GSSID=rep(NA, Clen)
for(mm in  1:Clen ){
    Group1[mm]=splitname[[mm]][1]
    Group2[mm]=splitname[[mm]][2]
    Group3[mm]=splitname[[mm]][3]
    GSSID[mm]=splitname[[mm]][4]
}
Group1[Group1=='BabyAdult']='BabyOutlier'
Group1[Group1=='NotUsed']='BabyForeskin'
Group1[Group1=='X20s']='Adult20s'
Group1[Group1=='X60s']='Adult60s'
Group1[Group1=='baby']='Baby'
C=t(B)
df_pca <- prcomp(C)

df_out <- as.data.frame(df_pca$x)
df_out$group2=Group2
df_out$group1=Group1
df_out$group3=Group3
df_out$name=Cname

rm(A)

library(ggplot2)
library(grid)
library(gridExtra)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


####Generate a dot plot##########
library(RColorBrewer)

####with label
df_pca <- prcomp(C)
df_out <- as.data.frame(df_pca$x)
df_out$group2=Group2
df_out$group1=Group1
df_out$group3=Group3
df_out$name=Cname
#install.packages("ggpmisc")
library("ggpmisc")
library('ggpubr')
#p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+ geom_point(size=6)+scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1"))+stat_chull(aes(color=Group1, fill=Group1), alpha=0.1, geom="polygon") 
    #+ stat_ellipse()
#p1

p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group1, label=name ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(Group1)), "Set1"))+theme_bw()
png(paste0(filename,".PCA_babyForskin.png"), width=800, height=800, res=120)
p1
dev.off()

df_pca2 <- prcomp(C[Group1 !="BabyOutlier",])
df_out2 <-as.data.frame(df_pca2$x)
group=Group1[Group1 !="BabyOutlier"]
name2=Cname[Group1 !="BabyOutlier"]
p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_NoOutlier.png"), width=800, height=800, res=120)
p1
dev.off()

p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_bw()
png(paste0(filename,".PCA_NoOutlier_BW.png"), width=800, height=800, res=120)
p1
dev.off()

df_pca2 <- prcomp(C[Group1 !="BabyOutlier" & Group1 !="BabyForeskin",])
df_out2 <-as.data.frame(df_pca2$x)
group=Group1[Group1 !="BabyOutlier" & Group1 != "BabyForeskin"]
name2=Cname[Group1 !="BabyOutlier"& Group1 != 'BabyForeskin']
p1<-ggplot(df_out2,aes(x=PC1,y=PC2,color=group, label=name2 ))+ geom_point(size=5)+scale_color_manual(values=brewer.pal(nlevels(as.factor(group)), "Set1"))+stat_chull(aes(color=group, fill=group), alpha=0.1, geom="polygon") +theme_pubclean()
png(paste0(filename,".PCA_NoForeskin.png"), width=800, height=800, res=120)
p1
dev.off()


##########################All the following not used######################
