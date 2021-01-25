library(clusterProfiler)
library(DOSE)
######Given a list tell you GO and KEGG matching U133 entrenz gene id
args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
COLN <-as.numeric(args[2])
chip<-args[3]

rm(args)
#filename="agree" ###remove long gene name, there are issure with it
A<-read.table(filename, sep="\t", header=TRUE)
 mygene<-as.character(unique(A[,COLN]))

if(chip =="U133"){
    B<-read.table("/home/ping/db/affy/NA36/HG-U133_Plus_2.na36.annot.csv.EZ", sep="\t", header=TRUE)
    allgene<-as.character(unique(B[,2]))
}else if (chip == "U219"){
    B<-read.table("/home/ping/db/affy/NA36/HG-U219.na36.annot.csv.EZ", sep="\t", header=TRUE)
    allgene<-as.character(unique(B[,2])) 
}
#5.2 GO classification "MF", "BP", and "CC" subontologies.
ggo <- groupGO(gene= mygene, organism="human",ont ="MF",level=3,readable = TRUE)
write.table(summary(ggo), file=paste0(filename, ".ggo"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
#5.3 GO over-representation test
#Gene ID can be mapped to gene Symbol by using paramter readable=TRUE or setReadable function.
ego <- enrichGO(gene= mygene, universe= allgene,organism="human",ont= "BP",pAdjustMethod = "BH",pvalueCutoff=0.001,qvalueCutoff=0.05, readable= TRUE)
write.table(summary(ego), file=paste0(filename, ".ego"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
cnetplot(ego, categorySize="geneNum")
####
kk <- enrichKEGG(gene= mygene,organism='human',pvalueCutoff = 0.05)
write.table(summary(kk), file=paste0(filename, ".kk"), sep="\t", row.name=TRUE, col.names = TRUE, eol = "\n", na = "NA")
