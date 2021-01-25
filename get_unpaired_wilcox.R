args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################
#filename="/mnt/G6_2D/project/BabySkin/org/GSS2656_data"
filename="GSS2656_data"
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=log10(A[1:d[1], 2:d[2]])
rownames(B)==A[,1]
Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[_]")
Group=rep(NA, Clen)
GT=rep(NA,Clen)
BT=rep(NA, Clen)
Site=rep(NA,Clen)
Gender=rep(NA,Clen)
ID=rep(NA,Clen)
GSSID=rep(NA, Clen)
for(mm in  1:Clen ){
    Group[mm]=splitname[[mm]][1]
    GT[mm]=splitname[[mm]][2]
    BT[mm]=splitname[[mm]][3]
    Site[mm]=splitname[[mm]][4]
    Gender[mm]=splitname[[mm]][5]
    ID[mm]=splitname[[mm]][6]
    GSSID[mm]=splitname[[mm]][7]

}
Group1=Group
Group1[Group=="FT2"|Group=="FT5"]="FT"
Group1[Group=="LPT2"|Group=="LPT5"]="LPT"
Group1[Group=="VPT2"|Group=="VPT5"]="VPT"
Group0=Group1
Group0[Group1=="FT"|Group1=="LPT"|Group1=="VPT"]="Baby"

df_pca <- prcomp(t(B))
#plot(df_pca$x[,1], df_pca$x[,2])
df_out <- as.data.frame(df_pca$x)
df_out$group=Group
df_out$group1=Group1
df_out$group0=Group0
df_out$name=paste0(BT,Site,GT,Gender,ID)
library(ggplot2)
library(grid)
library(gridExtra)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))


####Generate a dot plot##########
library(RColorBrewer)
p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group, label=name ))+geom_point(size=3)+ scale_color_manual(values=brewer.pal(nlevels(as.factor(Group)), "Dark2"))

####with label
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Group, label=name ))+geom_point() + scale_color_manual(values=brewer.pal(nlevels(as.factor(Group)), "Dark2")) + geom_text(size=8)+theme  #+stat_ellipse()

png("PCA_name.png", width=2800, height=2800, res=120)
print(p)
dev.off()

#p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=group))+geom_point()
#p
#p1<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=group1))+geom_point()
#p1
#p0<-ggplot(df_out,aes(x=PC1,y=PC2,color=group,shape=group0))+geom_point()
#p0

#x=rowMean(B)

#boxplot(B)
#install.packages('limma', repos='http://cran.us.r-project.org')

#install.packages("BiocManager", repos='http://cran.rstudio.org')

#install.packages("remote", repos='http://cran.rstudio.org')
#BiocManager::install("limma",site_repository='http://cran.rstudio.org' )

#BiocManager::install("Glimma",site_repository='http://cran.rstudio.org' )
library("cluster")
library("factoextra")
library("magrittr")

res.dist <- get_dist(t(B), stand = TRUE)
jpeg("Halo_KO_eu_dist.jpg", width = 3500, height = 3500)
fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()




library(RColorBrewer)
library(limma)
library(Glimma)
col.group <- as.factor(Group)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Dark2")
col.group <- as.character(col.group)

#plotMDS(B, labels=Group, col=col.group) Generate the web page for MDS data
glMDSPlot(B, labels=paste0(BT,Site,GT,Gender,ID), groups=Group, launch=FALSE)

library(factoextra)
library(NbClust)
#kmeans, pam, clara and hcut (for hierarchical clustering).
c=t(B)
fviz_nbclust(c, kmeans, method = c("silhouette", "wss", "gap_stat"))
#https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
  # Elbow method
  jpeg("kmeans_elbow.jpg", width = 350, height = 350)
  fviz_nbclust(c, kmeans, method = "wss") +
      geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")
  dev.off()

  # Silhouette method
  jpeg("kmeans_Silhouette.jpg", width = 350, height = 350)
  fviz_nbclust(c, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")
  dev.off()

  # Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  jpeg("kmeans_Gap_Statistic.jpg", width = 350, height = 350)
  set.seed(123)
  fviz_nbclust(c, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
    labs(subtitle = "Gap statistic method")
  dev.off()
  jpeg("BbClust_euclidean_2_15.jpg", width = 350, height = 350)
  NbClust(data = c, diss = NULL, distance = "euclidean",
          min.nc = 2, max.nc = 15, method = NULL)
  dev.off()

print(paste("genename",sep=","))

for (i in 1:d[1]){
    genename=A[i,1]
    gene<-as.numeric(B[i,])
    mydata=data.frame(gene,Group,Group1, Group0, GT, BT, Gender,ID, Site, GSSID )
    KP_Group=kruskal.test(gene ~ Group, data = mydata)$p.value
    KP_Group1=kruskal.test(gene ~ Group1, data = mydata)$p.value
    KP_Group0=kruskal.test(gene ~ Group0, data = mydata)$p.value
    mydataNO=mydata[BT!=49&BT!=28.9,]
    KPNO_Group=kruskal.test(gene ~ Group, data = mydataNO)$p.value
    KPNO_Group1=kruskal.test(gene ~ Group1, data = mydataNO)$p.value
    KPNO_Group0=kruskal.test(gene ~ Group0, data = mydataNO)$p.value

    PB_HN=my.wilcox.p.value(B_N,B_H)
    PG_HN=my.wilcox.p.value(G_N,G_H)
    PI_HN=my.wilcox.p.value(I_N,I_H)
    PP_HN=my.wilcox.p.value(P_N,P_H)
    P_HN=my.wilcox.p.value(R_N,R_H)
    FB_HN=M_B_H/M_B_N
    if(FB_HN<1){FB_HN=-1/FB_HN}
    FG_HN=M_G_H/M_G_N
    if(FG_HN<1){FG_HN=-1/FG_HN}
    FI_HN=M_I_H/M_I_N
    if(FI_HN<1){FI_HN=-1/FI_HN}
    FP_HN=M_P_H/M_P_N
    if(FP_HN<1){FP_HN=-1/FP_HN}
    F_HN=M_R_H/M_R_N
    if(F_HN<1){F_HN=-1/F_HN}
    CorS_All<-cor( gene, RScore,method="spearman")
    PCorS_All<-cor.test( gene, RScore,method="spearman", use="pairwise.complete.obs")$p.value
    CorP_All<-cor( gene,RScore,method="pearson")
    PCorP_All<-cor.test( gene, RScore,method="pearson", use="pairwise.complete.obs")$p.value

    CorS_B<-cor( B, RScore[BodySite=="Buttock"],method="spearman")
    PCorS_B<-cor.test( B, RScore[BodySite=="Buttock"],method="spearman", use="pairwise.complete.obs")$p.value
    CorP_B<-cor( B,RScore[BodySite=="Buttock"],method="pearson")
    PCorP_B<-cor.test( B, RScore[BodySite=="Buttock"],method="pearson", use="pairwise.complete.obs")$p.value

    CorS_G<-cor( G, RScore[BodySite=="Genital"],method="spearman")
    PCorS_G<-cor.test( G, RScore[BodySite=="Genital"],method="spearman", use="pairwise.complete.obs")$p.value
    CorP_G<-cor( G,RScore[BodySite=="Genital"],method="pearson")
    PCorP_G<-cor.test( G, RScore[BodySite=="Genital"],method="pearson", use="pairwise.complete.obs")$p.value

    CorS_I<-cor( I, RScore[BodySite=="Interiginous"],method="spearman")
    PCorS_I<-cor.test( I, RScore[BodySite=="Interiginous"],method="spearman", use="pairwise.complete.obs")$p.value
    CorP_I<-cor( I,RScore[BodySite=="Interiginous"],method="pearson")
    PCorP_I<-cor.test( I, RScore[BodySite=="Interiginous"],method="pearson", use="pairwise.complete.obs")$p.value

    CorS_P<-cor( P, RScore[BodySite=="Perianal"],method="spearman")
    PCorS_P<-cor.test( P, RScore[BodySite=="Perianal"],method="spearman", use="pairwise.complete.obs")$p.value
    CorP_P<-cor( P,RScore[BodySite=="Perianal"],method="pearson")
    PCorP_P<-cor.test( P, RScore[BodySite=="Perianal"],method="pearson", use="pairwise.complete.obs")$p.value

    maxCor=max(CorS_All, CorP_All, CorS_B, CorP_B, CorS_G, CorP_G, CorS_I, CorP_I, CorS_P, CorP_P)
    minCor=min(CorS_All, CorP_All, CorS_B, CorP_B, CorS_G, CorP_G, CorS_I, CorP_I, CorS_P, CorP_P)
    minPCor=min(PCorS_All, PCorP_All, PCorS_B, PCorP_B, PCorS_G, PCorP_G, PCorS_I, PCorP_I, PCorS_P, PCorP_P)
    minPWilcox=min(PB_HN, PG_HN, PI_HN, PP_HN, P_HN)
    mydata=data.frame(gene,log(gene), BodySite, Rash, RScore,SiteRash )
    KP_BodySite=kruskal.test(gene ~ BodySite, data = mydata)$p.value
    KP_Rash=kruskal.test(gene ~ Rash, data = mydata)$p.value
    KP_RScore=kruskal.test(gene ~ RScore, data = mydata)$p.value
    KP_SiteRash=kruskal.test(gene ~ SiteRash, data = mydata)$p.value

    Bdata=mydata[BodySite=="Buttock",]
    KP_Rash_B=kruskal.test(gene ~ Rash, data = Bdata)$p.value
    KP_RScore_B=kruskal.test(gene ~ RScore, data = Bdata)$p.value

    Gdata=mydata[BodySite=="Genital",]
    KP_Rash_G=kruskal.test(gene ~ Rash, data = Gdata)$p.value
    KP_RScore_G=kruskal.test(gene ~ RScore, data = Gdata)$p.value

    Idata=mydata[BodySite=="Interiginous",]
    KP_Rash_I=kruskal.test(gene ~ Rash, data = Idata)$p.value
    KP_RScore_I=kruskal.test(gene ~ RScore, data = Idata)$p.value

    Pdata=mydata[BodySite=="Perianal",]
    KP_Rash_P=kruskal.test(gene ~ Rash, data = Pdata)$p.value
    KP_RScore_P=kruskal.test(gene ~ RScore, data = Pdata)$p.value

    minPk=min(KP_BodySite, KP_Rash, KP_RScore, KP_SiteRash, KP_Rash_B, KP_RScore_B,  KP_Rash_G, KP_RScore_G,  KP_Rash_I, KP_RScore_I,  KP_Rash_P, KP_RScore_P)

    print(paste(genename,genus,minPWilcox,minPk,minPCor,maxCor, minCor,maxAve, KP_BodySite, KP_Rash, KP_RScore, KP_SiteRash, KP_Rash_B, KP_RScore_B,  KP_Rash_G, KP_RScore_G,  KP_Rash_I, KP_RScore_I,  KP_Rash_P, KP_RScore_P, PCorS_All, PCorP_All, PCorS_B, PCorP_B, PCorS_G, PCorP_G, PCorS_I, PCorP_I, PCorS_P, PCorP_P,
                PB_HN, PG_HN, PI_HN, PP_HN, P_HN, FB_HN, FG_HN, FI_HN, FP_HN, F_HN,CorS_All, CorP_All, CorS_B, CorP_B, CorS_G, CorP_G, CorS_I, CorP_I, CorS_P, CorP_P, M_R_N,M_R_L, M_R_H, M_B, M_G, M_I, M_P, M_B_N, M_B_L, M_B_H,  M_G_N, M_G_L, M_G_H,  M_I_N, M_I_L, M_I_H,  M_P_N, M_P_L, M_P_H,sep=","))

}
