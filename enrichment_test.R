library(bc3net)
inlist <- strsplit(readLines("/media/ping/_media_G6D/db/ko/kegg_tool_2017/pathway_ko.list"), '\t')
kopathways <- lapply(inlist, tail, n = -1)
names(kopathways) <- lapply(inlist, head, n = 1)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

x=read.table(filename, header=TRUE, sep='\t')
colnames(x)=c("KOID", "Pvalue", "FC")
present=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED"),]$KOID
sig=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED") & (x$Pvalue <=0.05),]$KOID
ko.result=enrichment(sig, present, kopathways, verbose=TRUE)
write.csv(ko.result,paste0(filename, ".enrichment.csv"))

presentU=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED")&(x$FC>0),]$KOID
sigU=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED") & (x$Pvalue <=0.05)&(x$FC >0),]$KOID
koU.result=enrichment(sigU, presentU, kopathways, verbose=TRUE)
write.csv(koU.result,paste0(filename, ".up.enrichment.csv"))


presentD=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED")&(x$FC<0),]$KOID
sigD=x[(x$KOID != "UNMAPPED") & (x$KOID != "UNGROUPED") & (x$Pvalue <=0.05)&(x$FC<0),]$KOID
koD.result=enrichment(sigD, presentD, kopathways, verbose=TRUE)
write.csv(koD.result,paste0(filename, ".down.enrichment.csv"))