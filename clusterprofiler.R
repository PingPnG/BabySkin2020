rm(list = ls())
#########################################
#Ping Hu
##########################################
library(clusterProfiler)
library(DOSE)
######Given a list tell you GO and KEGG matching U133 entrenz gene id
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#filename <- args[1]
#COLN <-as.numeric(args[2])
#chip<-args[3]
filename="20_60.stat.xls"
COLN=4 ###entrenzID column
chip="U219"

#rm(args)
#filename="agree" ###remove long gene name, there are issure with it
A<-read.table(filename, sep="\t", header=TRUE)
mygene<-as.character(unique(A[,COLN]))

#if(chip =="U133"){
#    B<-read.table("HG-U133_Plus_2.na36.Entrez", sep="\t", header=TRUE)
#    allgene<-as.character(unique(B[,2]))
#}else if (chip == "U219"){
    B<-read.table("HG-U219.na36.annot.csv.EZ", sep="\t", header=TRUE)
    allgene<-as.character(unique(B[,2])) 
#}

#5.3 GO over-representation test
#Gene ID can be mapped to gene Symbol by using paramter readable=TRUE or setReadable function.
ego <- enrichGO(gene= mygene, universe= allgene, 'org.Hs.eg.db',ont= "BP",pAdjustMethod = "BH",pvalueCutoff=0.001,qvalueCutoff=0.05, readable= TRUE)
write.table(summary(ego), file=paste0(filename, ".ego"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
#cnetplot(ego, categorySize="geneNum")

kk <- enrichKEGG(gene= mygene,'hsa',pvalueCutoff = 0.05)
write.table(summary(kk), file=paste0(filename, ".kk"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
#C5: GO gene sets
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt <- enricher(mygene, TERM2GENE=c5)
write.table(summary(egmt), file=paste0(filename, ".C5egmt"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
#wikipthway
wpgmtfile <- system.file("extdata", "wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
library(tidyr)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

ewp <- enricher(mygene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
ewp <- setReadable(ewp, org.Hs.eg.db, keytype = "ENTREZID")
write.table(summary(ewp), file=paste0(filename, ".ewp"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
