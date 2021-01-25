#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
#http://yulab-smu.top/clusterProfiler-book/chapter12.html#upset-plot


rm(list=ls())
setwd("/Users/ping/Desktop/BabySkin/Results_10_6/Primary_Results/")
data <- read.csv(file = 'FTvsLPT.all.summary.csv')
library("hgu219.db")

colnames(data)[1]='probe_id'
egids2=hgu219ENTREZID[data$probe_id]

annots=toTable(egids2)
total <- merge(data,annots,by="probe_id")
library(dplyr)
####only use the single gene with highest average expression
x<- total %>% group_by(gene_id) %>% filter(AveExpr == max(AveExpr)) 

remove(data)
###significant gene list
Sig<-x[x$P.Value <=0.05,]
foldchanges=Sig$logFC
names(foldchanges)=Sig$gene_id
#'quantreg''RcppArmadillo' 'statmod'
#installation of both Clang and GFortran solved my problem.
#https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
library(magrittr)
library(clusterProfiler)
library(tidyr)
genelist<-total$logFC
names(genelist) <-total$gene_id
genelist <- sort(genelist, decreasing = TRUE)
library("DOSE")
data(geneList, package="DOSE")
head(genelist)
gene<-names(foldchanges)
allgenes <-total$gene_id

gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
ggo_cc <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
ggo_bp <- groupGO(gene     = gene,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  level    = 3,
                  readable = TRUE)
ggo_mf <- groupGO(gene     = gene,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "MF",
                  level    = 3,
                  readable = TRUE)
#head(ggo)
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