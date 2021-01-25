## ----setup, message=FALSE, echo = FALSE-------------------------------------------------
options(digits=3)
options(width=90)

## ----setup2, eval=TRUE, message=FALSE---------------------------------------------------
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

## ----downloadData, eval=TRUE------------------------------------------------------------
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

## ----import1----------------------------------------------------------------------------
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

## ----import2----------------------------------------------------------------------------
x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)

## ----annotatesamples--------------------------------------------------------------------
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

## ----annotategenes, message=FALSE-------------------------------------------------------
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

## ----removedups-------------------------------------------------------------------------
genes <- genes[!duplicated(genes$ENTREZID),]

## ----assigngeneanno---------------------------------------------------------------------
x$genes <- genes
x

## ----cpm--------------------------------------------------------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

## ----lcpm-------------------------------------------------------------------------------
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

## ----zeroes-----------------------------------------------------------------------------
table(rowSums(x$counts==0)==9)

## ----filter-----------------------------------------------------------------------------
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## ----filterplot1, fig.height=4, fig.width=8, fig.cap="The density of log-CPM values for raw pre-filtered data (A) and post-filtered data (B) are shown for each sample. Dotted vertical lines mark the log-CPM threshold (equivalent to a CPM value of about 0.2) used in the filtering step."----
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

## ----normalize--------------------------------------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## ----normalizemodifieddata--------------------------------------------------------------
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

## ----plotmodifieddata, fig.height=4, fig.width=8, fig.cap="Example data: Boxplots of log-CPM values showing expression distributions for unnormalised data (A) and normalised data (B) for each sample in the modified dataset where the counts in samples 1 and 2 have been scaled to 5% and 500% of their original values respectively."----
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

## ----MDS1, fig.height=4, fig.width=8, fig.cap="MDS plots of log-CPM values over dimensions 1 and 2 with samples coloured and labeled by sample groups (A) and over dimensions 3 and 4 with samples coloured and labeled by sequencing lane (B). Distances on the plot correspond to the leading fold-change, which is the average (root-mean-square) log2-fold-change for the 500 genes most divergent between each pair of samples by default."----
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

## ----GlimmaMDSplot----------------------------------------------------------------------
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)

## ----design-----------------------------------------------------------------------------
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts--------------------------------------------------------------------------
contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP, 
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML, 
   levels = colnames(design))
contr.matrix

## ----voom, fig.height=4, fig.width=8, fig.cap="Means (x-axis) and variances (y-axis) of each gene are plotted to show the dependence between the two before `voom` is applied to the data (left panel) and how the trend is removed after `voom` precision weights are applied to the data (right panel). The plot on the left is created within the `voom` function which extracts residual variances from fitting linear models to log-CPM transformed data. Variances are then rescaled to quarter-root variances (or square-root of standard deviations) and plotted against the average log2 count for each gene. The plot on the right is created using `plotSA` which plots log2 residual standard deviations against mean log-CPM values. In both plots, each black dot represents a gene. On the left plot, the red curve shows the estimated mean-variance trend used to compute the voom weights. On the right plot, the average log2 residual standard deviation estimated by the empirical Bayes algorithm is marked by a horizontal blue line. "----
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## ----decidetests------------------------------------------------------------------------
summary(decideTests(efit))

## ----treat------------------------------------------------------------------------------
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

## ----venn, fig.height=6, fig.width=6, fig.cap="Venn diagram showing the number of genes DE in the comparison between basal versus LP only (left), basal versus ML only (right), and the number of genes that are DE in both comparisons (center). The number of genes that are not DE in either comparison are marked in the bottom-right."----
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

## ----toptables--------------------------------------------------------------------------
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

## ----MDplot, fig.keep='none'------------------------------------------------------------
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

## ----GlimmaMDplot-----------------------------------------------------------------------
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)

## ----heatmap, fig.height=8, fig.width=5, fig.cap="Heatmap of log-CPM values for top 100 genes DE in basal versus LP. Expression across each gene (or row) have been scaled so that mean expression is zero and standard deviation is one. Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes have been reordered by the method of hierarchical clustering. A dendrogram is shown for the sample clustering.", message=FALSE----
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
   labRow=v$genes$SYMBOL[i], labCol=group, 
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")

## ----camera-----------------------------------------------------------------------------
load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
idx <- ids2indices(Mm.c2,id=rownames(v))
cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
head(cam.BasalvsLP,5)
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
head(cam.BasalvsML,5)
cam.LPvsML <- camera(v,idx,design,contrast=contr.matrix[,3])
head(cam.LPvsML,5)

## ----barcodeplot, fig.height=6, fig.width=6, fig.cap="Barcode plot of `LIM_MAMMARY_LUMINAL_MATURE_UP` (red bars, top of plot) and `LIM_MAMMARY_LUMINAL_MATURE_DN` (blue bars, bottom of plot) gene sets in the LP versus ML contrast. For each set, an enrichment line that shows the relative enrichment of the vertical bars in each part of the plot is displayed. The experiment of Lim *et al.* [@Lim:BreastCancerRes:2010] is very similar to the current one, with the same sorting strategy used to obtain the different cell populations, except that microarrays were used instead of RNA-seq to profile gene expression. Note that the inverse correlation (the up gene set is down and the down gene set is up) is a result of the way the contrast has been set up (LP versus ML) -- if reversed, the directionality would agree."----
barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")

## ----softwareinfo-----------------------------------------------------------------------
sessionInfo()

