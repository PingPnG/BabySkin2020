## ----setup, echo=FALSE, results="hide", warning=FALSE-------------------------
knitr::opts_chunk$set(prompt = TRUE, comment = NA)
suppressPackageStartupMessages({library(RnaSeqGeneEdgeRQL)})

## ----sra----------------------------------------------------------------------
library(RnaSeqGeneEdgeRQL)
targetsFile <- system.file("extdata", "targets.txt",
                           package="RnaSeqGeneEdgeRQL")
targets <- read.delim(targetsFile, stringsAsFactors=FALSE)
targets

## ----group--------------------------------------------------------------------
group <- paste(targets$CellType, targets$Status, sep=".")
group <- factor(group)
table(group)

## ----checkdownload, echo=FALSE, results="hide", message=FALSE-----------------
if( !file.exists("GSE60450_Lactation-GenewiseCounts.txt.gz") ) {
FileURL <- paste(
  "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
  "format=file",
  "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
  sep="&")
download.file(FileURL, method="libcurl", "GSE60450_Lactation-GenewiseCounts.txt.gz")
}

## ----download, eval=FALSE-----------------------------------------------------
#  FileURL <- paste(
#    "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
#    "format=file",
#    "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
#    sep="&")
#  download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")

## ----readcounts---------------------------------------------------------------
GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",
                             row.names="EntrezGeneID")
colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)
dim(GenewiseCounts)
head(GenewiseCounts)

## ----DGEList, message=FALSE---------------------------------------------------
library(edgeR)
y <- DGEList(GenewiseCounts[,-1], group=group,
             genes=GenewiseCounts[,1,drop=FALSE])
options(digits=3)
y$samples

## ----symbols, message=FALSE---------------------------------------------------
library(org.Mm.eg.db)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENTREZID", column="SYMBOL")
head(y$genes)

## ----dropNAsymbols------------------------------------------------------------
y <- y[!is.na(y$genes$Symbol), ]
dim(y)

## ----design-------------------------------------------------------------------
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

## ----keep---------------------------------------------------------------------
keep <- filterByExpr(y, design)
table(keep)

## ----filter-------------------------------------------------------------------
y <- y[keep, , keep.lib.sizes=FALSE]

## ----aveLogCPM, fig.cap="Histogram of average log2 CPM."----------------------
AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

## ----norm---------------------------------------------------------------------
y <- calcNormFactors(y)
y$samples

## ----mdsplot, fig.width=7, fig.height=7, fig.cap="MDS plot showing distances between expression profiles. This provides a type of unsupervised clustering of the samples. In this case, the first dimension separates samples by cell type while the second dimension corresponds to mouse status."----
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

## ----mdplot1, fig.width=7, fig.height=7, fig.cap="Mean-difference plot of log2-expression in sample 1 versus the average log2-expression across all other samples. Each point represents a gene, and the red line indicates a log-ratio of zero. The majority of points cluster around the red line."----
plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

## ----mdplot11, fig.width=7, fig.height=7, fig.cap="Mean-difference plot of log2-expression in sample 11 versus the average log2-expression across all other samples. The plot shows a number of genes that are both highly expressed and highly up-regulated."----
plotMD(y, column=11)
abline(h=0, col="red", lty=2, lwd=2)

## ----estimateDisp-------------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)

## ----plotBCV, fig.width=7, fig.height=7, fig.cap="Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene. The plot shows the square-root estimates of the common, trended and tagwise NB dispersions."----
plotBCV(y)

## ----glmQLFit-----------------------------------------------------------------
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

## ----QLDisp, fig.width=7, fig.height=7, fig.cap="A plot of the quarter-root QL dispersion against the average abundance of each gene. Estimates are shown for the raw (before EB moderation), trended and squeezed (after EB moderation) dispersions. Note that the QL dispersions and trend shown here are relative to the NB dispersion trend show in the previous figure produced using `plotBCV()`."----
plotQLDisp(fit)

## ----df.prior-----------------------------------------------------------------
summary(fit$df.prior)

## ----B.PvsL-------------------------------------------------------------------
B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)

## ----glmQLFTest---------------------------------------------------------------
res <- glmQLFTest(fit, contrast=B.LvsP)

## ----topTags------------------------------------------------------------------
topTags(res)

## ----decideTests--------------------------------------------------------------
is.de <- decideTestsDGE(res)
summary(is.de)

## ----plotMDfit, fig.width=7, fig.height=7, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Significantly up and down DE genes are highlighted in red and blue, respectively."----
plotMD(res, status=is.de)

## ----treat--------------------------------------------------------------------
tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
topTags(tr)

## ----treatdecideTests---------------------------------------------------------
is.de <- decideTestsDGE(tr)
summary(is.de)

## ----plotMDtreat, fig.width=7, fig.height=7, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Genes with fold-changes significantly greater than 1.5 are highlighted."----
plotMD(tr, status=is.de)

## ----cpm----------------------------------------------------------------------
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")

## ----order--------------------------------------------------------------------
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:30],]

## ----heatmap, message=FALSE, fig.width=8, fig.height=12, fig.cap="Heat map across all the samples using the top 30 most DE genes between the basal lactating group and the basal pregnancy group."----
coolmap(logCPM, margins=c(7,7), lhei=c(1,6), lwid=c(1,3))

## ----makeContrasts------------------------------------------------------------
con <- makeContrasts(
     L.PvsL = L.pregnant - L.lactating,
     L.VvsL = L.virgin - L.lactating,
     L.VvsP = L.virgin - L.pregnant, levels=design)

## ----anovaQLFtest-------------------------------------------------------------
res <- glmQLFTest(fit, contrast=con)
topTags(res)

## ----complicatedContrasts-----------------------------------------------------
con <- makeContrasts(
     (L.lactating-L.pregnant)-(B.lactating-B.pregnant), 
     levels=design)

## ----complicatedQLTest--------------------------------------------------------
res <- glmQLFTest(fit, contrast=con)
topTags(res)

## ----goana--------------------------------------------------------------------
go <- goana(tr, species="Mm")
topGO(go, n=15)

## ----kegga--------------------------------------------------------------------
keg <- kegga(tr, species="Mm")
topKEGG(keg, n=15, truncate=34)

## ----select, message=FALSE----------------------------------------------------
library(GO.db)
cyt.go <- c("GO:0032465", "GO:0000281")
term <- select(GO.db, keys=cyt.go, columns="TERM")
term

## ----GO2ALLEGS----------------------------------------------------------------
Rkeys(org.Mm.egGO2ALLEGS) <- cyt.go
cyt.go.genes <- as.list(org.Mm.egGO2ALLEGS)

## ----fry----------------------------------------------------------------------
B.VvsL <- makeContrasts(B.virgin-B.lactating, levels=design)
fry(y, index=cyt.go.genes, design=design, contrast=B.VvsL)

## ----barcode, fig.width=8, fig.height=5, fig.cap="Barcode plot showing enrichment of the GO term GO:0032465 in the basal virgin group compared to the basal lactating group. X-axis shows logFC for B.virgin vs B.lactating. Black bars represent genes annotated with the GO term. The worm shows relative enrichment."----
res <- glmQLFTest(fit, contrast=B.VvsL)
index <- rownames(fit) %in% cyt.go.genes[[1]]
barcodeplot(res$table$logFC, index=index, labels=c("B.lactating","B.virgin"), 
            main=cyt.go[1])

## ----camera-------------------------------------------------------------------
Mm.c2 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds"))

## ----id2indices.Mm.c2---------------------------------------------------------
idx <- ids2indices(Mm.c2,id=row.names(y))

## ----cam.BVvsLV---------------------------------------------------------------
BvsL.v <- makeContrasts(B.virgin - L.virgin, levels=design)
cam <- camera(y, idx, design, contrast=BvsL.v, inter.gene.cor=0.01)
options(digits=2)
head(cam,14)

## ----barcode2, fig.width=8, fig.height=6.4, fig.cap="Barcode plot showing strong enrichment of mammary stem cell signature in the stem cell vs luminal cell comparison. Red bars show up signature genes, blue bars show down genes. The worms show relative enrichment."----
res <- glmQLFTest(fit, contrast=BvsL.v)
barcodeplot(res$table$logFC,
            index=idx[["LIM_MAMMARY_STEM_CELL_UP"]],
            index2=idx[["LIM_MAMMARY_STEM_CELL_DN"]],
            labels=c("L.virgin","B.virgin"),
            main="LIM_MAMMARY_STEM_CELL",
            alpha=1)

## ----info---------------------------------------------------------------------
sessionInfo()

