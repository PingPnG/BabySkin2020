args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################
#filename="/mnt/G6_2D/project/BabySkin/org/GSS2656_data"
filename="selected_data"
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
B=log10(A[1:d[1], 5:d[2]])
rownames(B)=A[,1]
Cname=colnames(A)[5:d[2]]
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
Group2=Group0
Group2[Group1=="A20s"|Group1=="A60s"]="Adult"
Gall=unique(c(Group, Group1, Group2))
# [1] "A20s"  "A60s"  "FT2"   "FT5"   "LPT2"  "LPT5"  "VPT2"  "VPT5"  "FT"   
#[10] "LPT"   "VPT"   "Adult" "Baby" 


print(paste( "genename","genesymbol", "genename", "ezID",  "mean", paste(Gall, collapse=","),"sd", paste(Gall, collapse=","),"minP", "minAllP","minBabyP", "kruskal-pvalue", "KP_Group", "KP_Group1", "KP_Group0", "KPNO_Group", "KPNO_Group1", "KPNO_Group0",  "KPbaby_Group", "KPbaby_Group1", "KPbabyNO_Group", "KPbabyNO_Group1", sep=","))
   
for (i in 1:d[1]){
    genename=A[i,1]
    genesymbol=A[i,2]
    longname=A[i,3]
    ezid=A[i,4]
    gene<-as.numeric(B[i,])
    gene0<-as.numeric(A[i,5:dim(A)[2]])
    mydata=data.frame(gene,Group,Group1, Group0, GT, BT, Gender,ID, Site, GSSID )
    KP_Group=kruskal.test(gene ~ Group, data = mydata)$p.value
    KP_Group1=kruskal.test(gene ~ Group1, data = mydata)$p.value
    KP_Group0=kruskal.test(gene ~ Group0, data = mydata)$p.value
    mydataNO=mydata[BT!=49&BT!=28.9,]
    KPNO_Group=kruskal.test(gene ~ Group, data = mydataNO)$p.value
    KPNO_Group1=kruskal.test(gene ~ Group1, data = mydataNO)$p.value
    KPNO_Group0=kruskal.test(gene ~ Group0, data = mydataNO)$p.value
    mybaby=mydata[Group0=="Baby",]
    mybabyNO=mybaby[BT!=49&BT!=28.9,]
    KPbaby_Group=kruskal.test(gene ~ Group, data = mybaby)$p.value
    KPbaby_Group1=kruskal.test(gene ~ Group1, data = mybaby)$p.value
    KPbabyNO_Group=kruskal.test(gene ~ Group, data = mybabyNO)$p.value
    KPbabyNO_Group1=kruskal.test(gene ~ Group1, data = mybabyNO)$p.value
 
    Amean0=aggregate(x = gene0,by = list(Group),FUN=mean)
    Amean1=aggregate(x = gene0,by = list(Group1),FUN=mean)
    Amean2=aggregate(x = gene0,by = list(Group2),FUN=mean)
    Asd0=aggregate(x = gene0,by = list(Group),FUN=sd)
    Asd1=aggregate(x = gene0,by = list(Group1),FUN=sd)
    Asd2=aggregate(x = gene0,by = list(Group2),FUN=sd)
    
   x1=rbind(Amean0, Amean1, Amean2)
   y1=rbind(Asd0, Asd1, Asd2)
   z1=cbind(x1,y1)
   z1=unique(z1)
   colnames(z1)=c("Group", "aveSignal", "Group2", "sd")
   z1$Group2[startsWith(z1$Group2,"A")]="Adult" 
   z1$Group2[startsWith(z1$Group2,"F")]="FullTermBaby"
   z1$Group2[startsWith(z1$Group2,"V")]="VeryPrematureBaby"
   z1$Group2[startsWith(z1$Group2,"L")]="LatePrematureBaby"
   meanstr=paste(z1$aveSignal, collapse=",")
   sdstr=paste(z1$sd, collapse=",")
  
   #http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization 
#   library(ggplot2)
#   jpeg(paste0(genesymbol,".",genename,".jpg"))
   
#   p <- ggplot(z1, aes(x=Group, y=aveSignal,group=Group2, color=Group2))+geom_pointrange(aes( ymin=aveSignal-sd, ymax=aveSignal+sd)) +theme_classic()+labs(title=paste0(genename," | ", genesymbol, " | ", longname), x="Group of Samples", y = paste0(genename," ", genesymbol," Signal"))
#   print(p)
   
#   dev.off()

   #z1$Group
   #[1] "A20s"  "A60s"  "FT2"   "FT5"   "LPT2"  "LPT5"  "VPT2"  "VPT5"  "FT"   
   #[10] "LPT"   "VPT"   "Adult" "Baby" 
   # have to make sure that z1$Group== Gall

   minP=min( KP_Group, KP_Group1, KP_Group0, KPNO_Group, KPNO_Group1, KPNO_Group0,  KPbaby_Group, KPbaby_Group1, KPbabyNO_Group, KPbabyNO_Group1)
   minbabyP=min( KPbaby_Group, KPbaby_Group1, KPbabyNO_Group, KPbabyNO_Group1)
   
   print(paste(genename,genesymbol, longname, ezid, "mean",meanstr ,"sd", sdstr,  "minP", minP,minbabyP, "kruskal-pvalue",  KP_Group, KP_Group1, KP_Group0, KPNO_Group, KPNO_Group1, KPNO_Group0,  KPbaby_Group, KPbaby_Group1, KPbabyNO_Group, KPbabyNO_Group1, sep=","))
   
} 
