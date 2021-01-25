##########################
#filename="metaphlan.estimated_number.ann.xls.rearranged"
#filename="genefamily.community.relab10.ann.xls.rearranged"
#filename="baby_antimicrobial_protein.xls"
#filename="Baby.signal.used.filtered.xls"
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#GSS<-args[2]

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
B=A[1:d[1], 4:d[2]]
#ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
#C=B+ZZ 
C=B
Cname=colnames(A)[4:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
BabyAdult=rep(NA, Clen)
BabyCat=rep(NA,Clen)
id=rep(NA,Clen)
FTPTCat=rep(NA,Clen)
gender=rep(NA,Clen)
site=rep(NA, Clen)
chipid=rep(NA,Clen)
for(mm in  1:Clen ){
    BabyAdult[mm]=splitname[[mm]][1]
    BabyCat[mm]=splitname[[mm]][2]
    FTPTCat[mm]=splitname[[mm]][3]
    id[mm]=splitname[[mm]][6]
    site[mm]=splitname[[mm]][5]
    gender[mm]=splitname[[mm]][4]
    chipid[mm]=splitname[[mm]][7]
}

      
    print(paste("probename","genesymbol", "genename", "MA20","MA60","MBaby","MFT","MVeryPT","MFT2","MPT2","MPT5",
    "Pwilcox_A20vA60","Pwilcox_BabyvA20","Pwilcox_VeryPTvFT",
    "Pwilcox_PT2vFT2","Pwilcox_PT5vFT2","Pwilcox_PT2vPT5",
    "Pt2_A20vA60","Pt2_BabyvA20","Pt2_VeryPTvFT",
    "Pt2_PT2vFT2","Pt2_PT5vFT2","Pt2_PT2vPT5",
    "truefc_A20vA60","truefc_BabyvA20","truefc_VeryPTvFT",
    "truefc_PT2vFT2","truefc_PT5vFT2","truefc_PT2vPT5",
    "minPwilcox", "minPt", "maxAve","maxAbS_FC",
    sep=","))
 
for (i in 1:d[1]){
    probename=A[i,1]
    genesymbol=A[i,2]
    genename=A[i,3]
    gene<-as.numeric(C[i,])
    A20<-gene[BabyAdult=="X20s"]
    A60<-gene[BabyAdult=="X60s"]
    Baby<-gene[BabyAdult=="baby"]
    FT<-gene[BabyCat=="FT"]
    VeryPT<-gene[BabyCat=="veryPT"]
    FT2<-gene[FTPTCat=="FT2"]
    PT2<-gene[FTPTCat=="PT2"]
    PT5<-gene[FTPTCat=="PT5"]
   
    MA20=mean(A20)
    MA60=mean(A60)
    MBaby=mean(Baby)
    MFT=mean(FT)
    MVeryPT=mean(VeryPT)
    MFT2=mean(FT2)
    MPT2=mean(PT2)
    MPT5=mean(PT5)
    
    Pwilcox_A20vA60<-my.wilcox.p.value(A20, A60, na.rm=TRUE)
    Pt2_A20vA60<-my.t.test.p.value(log(A20), log(A60), na.rm=TRUE)
    FC_A20vA60<-MA20/MA60
    truefc_A20vA60=truefc(FC_A20vA60)

    Pwilcox_BabyvA20<-my.wilcox.p.value(Baby, A20, na.rm=TRUE)
    Pt2_BabyvA20<-my.t.test.p.value(log(Baby), log(A20), na.rm=TRUE)
    FC_BabyvA20<-MBaby/MA20
    truefc_BabyvA20=truefc(FC_BabyvA20)
   
    Pwilcox_VeryPTvFT<-my.wilcox.p.value(VeryPT, FT, na.rm=TRUE)
    Pt2_VeryPTvFT<-my.t.test.p.value(log(VeryPT), log(FT), na.rm=TRUE)
    FC_VeryPTvFT<-MVeryPT/MFT
    truefc_VeryPTvFT=truefc(FC_VeryPTvFT)
    
    Pwilcox_PT2vFT2<-my.wilcox.p.value(PT2, FT2, na.rm=TRUE)
    Pt2_PT2vFT2<-my.t.test.p.value(log(PT2), log(FT2), na.rm=TRUE)
    FC_PT2vFT2<-MPT2/MFT2
    truefc_PT2vFT2=truefc(FC_PT2vFT2)

    Pwilcox_PT5vFT2<-my.wilcox.p.value(PT5, FT2, na.rm=TRUE)
    Pt2_PT5vFT2<-my.t.test.p.value(log(PT5), log(FT2), na.rm=TRUE)
    FC_PT5vFT2<-MPT5/MFT2
    truefc_PT5vFT2=truefc(FC_PT5vFT2)

    Pwilcox_PT2vPT5<-my.wilcox.p.value(PT2, PT5, na.rm=TRUE)
    Pt2_PT2vPT5<-my.t.test.p.value(log(PT2), log(PT5), na.rm=TRUE)
    FC_PT2vPT5<-MPT2/MPT5
    truefc_PT2vPT5=truefc(FC_PT2vPT5)

    maxAve=max(MA20,MA60,MBaby,MFT,MVeryPT,MFT2,MPT2,MPT5)
    minPwilcox=min( Pwilcox_A20vA60,Pwilcox_BabyvA20,Pwilcox_VeryPTvFT,
    Pwilcox_PT2vFT2,Pwilcox_PT5vFT2,Pwilcox_PT2vPT5)
    
    minPt=min(Pt2_A20vA60,Pt2_BabyvA20,Pt2_VeryPTvFT,
    Pt2_PT2vFT2,Pt2_PT5vFT2,Pt2_PT2vPT5)
    
    maxFC=max(truefc_A20vA60,truefc_BabyvA20,truefc_VeryPTvFT,
              truefc_PT2vFT2,truefc_PT5vFT2,truefc_PT2vPT5)
    minFC=min(truefc_A20vA60,truefc_BabyvA20,truefc_VeryPTvFT,
              truefc_PT2vFT2,truefc_PT5vFT2,truefc_PT2vPT5) 
    maxAbS_FC=max(abs(maxFC), abs(minFC))

    print(paste(probename,genesymbol, genename, MA20,MA60,MBaby,MFT,MVeryPT,MFT2,MPT2,MPT5, 
    Pwilcox_A20vA60,Pwilcox_BabyvA20,Pwilcox_VeryPTvFT,
    Pwilcox_PT2vFT2,Pwilcox_PT5vFT2,Pwilcox_PT2vPT5, 
    Pt2_A20vA60,Pt2_BabyvA20,Pt2_VeryPTvFT,
    Pt2_PT2vFT2,Pt2_PT5vFT2,Pt2_PT2vPT5,
    truefc_A20vA60,truefc_BabyvA20,truefc_VeryPTvFT,
    truefc_PT2vFT2,truefc_PT5vFT2,truefc_PT2vPT5,
    minPwilcox, minPt, maxAve,maxAbS_FC,
    sep=","))
   
} 
