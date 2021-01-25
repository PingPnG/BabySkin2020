#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
#http://yulab-smu.top/clusterProfiler-book/chapter12.html#upset-plot

#https://f1000research.com/articles/5-1384/v1
rm(list=ls())
filename="GSS2656.stat.xls"
A<-read.csv(filename, sep="\t", header=TRUE)
d=dim(A)
###original 34442 gene
library("hgu219.db")
colnames(A)[1]='probe_id'
egids2=hgu219ENTREZID[A$probe_id]
annots=toTable(egids2)
x<- merge(A,annots,by="probe_id")
library(dplyr)
####only use the single gene with highest average expression
#x<- total %>% group_by(gene_id) %>% filter(AverageSignal == max(AverageSignal)) 
#dim(x)
#### reduce into 14018 gene
#colnames(A)
#[1] "probe"                 "syn"                   "GeneName"              "EZID"                 
#[5] "KP"                    "wilcox.Adult_Baby"     "wilcox.20_60"          "wilcox.FT_LPT"        
#[9] "wilcox.FT3_PT3"        "wilcox.FT3_PT5"        "wilcox.PT3_PT5"        "truefc.Adult_Baby"    
#[13] "truefc.20_60"          "truefc.FT_LPT"         "truefc.FT3_PT3"        "truefc.FT3_PT5"       
#[17] "truefc.PT3_PT5"                        "fdr.KP"               
#[29] "fdr.wilcox.Adult_Baby" "fdr.wilcox.20_60"    
###significant gene list
genelist<-x$wilcox.FT3_PT5[x$wilcox.FT3_PT5 <=0.05]
names(genelist) <-x$EZID[x$wilcox.FT3_PT5 <=0.05]
genelist <- sort(genelist, decreasing = TRUE)
library("DOSE")
data(geneList, package="DOSE")

foldchanges=x$truefc.FT_LPT[x$wilcox.FT3_PT5 <=0.05]
names(foldchanges)=x$EZID[x$wilcox.FT3_PT5 <=0.05]
gene<-names(foldchanges)
geneU<-names(foldchanges [foldchanges >0])
geneD<-names(foldchanges [foldchanges <0])
allgenes<-x$EZID
library(magrittr)
library(clusterProfiler)
library(tidyr)
#gene.df <- bitr(gene, fromType = "ENTREZID",
#                toType = c("ENSEMBL", "SYMBOL"),
#                OrgDb = org.Hs.eg.db)
####this type of method only have ration and number
#ggo_bp <- groupGO(gene     = gene,
#                  OrgDb    = org.Hs.eg.db,
#                  ont      = "BP",
#                  level    = 3,
#                  readable = TRUE)
#barplot(ggo_bp, showCategory = 100)
ego_bp <- enrichGO(gene          = geneU,
                   universe      = allgenes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
dotplot(ego_bp,showCategory = 20 )
barplot(ego_bp, ,showCategory = 50)
write.table(summary(ego_bp), file="FT_LPT.ego_bp", sep="\t", row.name=TRUE
            , col.names = TRUE, eol = "\n", na = "NA")


ego3_bp <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
dotplot(ego3_bp, ,showCategory = 20)
#'quantreg''RcppArmadillo' 'statmod'
#installation of both Clang and GFortran solved my problem.
#https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/

#cnetplot(ego, categorySize="geneNum")
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
head(keggres)
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways  
pathview(gene.data=foldchanges, pathway.id="hsa04520", species="hsa")
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
write.table(gobpres, file="FT3_PT5.gobp.txt",sep="\t")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")
library(STRINGdb)
###use gene symbol to match
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )
enrichment <- string_db$get_enrichment( hits )

ego_cc <- enrichGO(gene          = gene,
                universe      = allgenes,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego_bp <- enrichGO(gene          = gene,
                   universe      = allgenes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego_mf <- enrichGO(gene          = gene,
                   universe      = allgenes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
#head(ego)
ego3_cc <- gseGO(geneList     = genelist,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego3_bp <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego3_mf <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk2 <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
dotplot(kk2)
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa')
mkk2 <- gseMKEGG(geneList = genelist,
                 organism = 'hsa')
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt <- enricher(gene, TERM2GENE=c5)
library(meshes)
x <- enrichMeSH(gene, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', category = 'C')
head(x)

dotplot(ego3_bp)
dotplot(kk2)
barplot(ego_bp, showCategory = 20)
library(ggplot2)
dotplot(ego_bp, showCategory=30) + ggtitle("dotplot for GO BP")

edox <- setReadable(ego_bp, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=genelist)
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(edox, categorySize="pvalue", foldChange=genelist)
cnetplot(edox, foldChange=genelist, circular = TRUE, colorEdge = TRUE)
#library(cowplot)
#cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
#Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
#                     Viewport has zero dimension(s)
cnetplot(edox, node_label="category") 
cnetplot(edox, node_label="gene")
cnetplot(edox, node_label="all") 
cnetplot(edox, node_label="none")
heatplot(edox)
heatplot(edox, foldChange=genelist)
emapplot(edox)
emapplot(edox, pie_scale=1.5)
emapplot(edox,layout="kk")
emapplot(edox, pie_scale=1.5,layout="kk") 
ridgeplot(kk2)
gseaplot(kk2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
gseaplot(kk2, geneSetID = 1, by = "preranked", title = kk2$Description[1])
gseaplot2(kk2, geneSetID = 1, title = kk2$Description[1])


#gseaplot2(kk2, geneSetID = 1, title = kk2$Description[1])
gsearank(kk2, 1, title = kk2[1, "Description"])

#terms <- kk2$Description[1:3]
#pmcplot(terms, 2010:2017)
goplot(ego_bp)

library("pathview")
hsa04110 <- pathview(gene.data  = genelist,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(genelist)), cpd=1))
install.packages("UpSetR")