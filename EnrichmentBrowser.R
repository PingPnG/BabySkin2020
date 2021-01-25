## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex()
knitr::opts_chunk$set(message = FALSE)

## ----setup, echo = FALSE---------------------------------------------------
library(EnrichmentBrowser)
library(ALL)
library(airway)

## ----readSE----------------------------------------------------------------
library(EnrichmentBrowser)
data.dir <- system.file("extdata", package = "EnrichmentBrowser")
exprs.file <- file.path(data.dir, "exprs.tab")
cdat.file <- file.path(data.dir, "colData.tab")
rdat.file <- file.path(data.dir, "rowData.tab")
se <- readSE(exprs.file, cdat.file, rdat.file)

## ----help, eval=FALSE------------------------------------------------------
#  ?readSE
#  ?SummarizedExperiment

## ----sexp2eset-------------------------------------------------------------
eset <- as(se, "ExpressionSet")

## ----eset2sexp-------------------------------------------------------------
se <- as(eset, "SummarizedExperiment")

## ----load-ALL--------------------------------------------------------------
library(ALL)
data(ALL)

## ----subset-ALL------------------------------------------------------------
ind.bs <- grep("^B", ALL$BT)
ind.mut <- which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
sset <- intersect(ind.bs, ind.mut)
all.eset <- ALL[, sset]

## ----show-ALL--------------------------------------------------------------
dim(all.eset)
exprs(all.eset)[1:4,1:4]

## ----probe2gene------------------------------------------------------------
allSE <- probe2gene(all.eset)
head(rownames(allSE))

## ----show-probe2gene-------------------------------------------------------
rowData(se)

## ----load-airway-----------------------------------------------------------
library(airway)
data(airway)

## ----preproc-airway--------------------------------------------------------
airSE <- airway[grep("^ENSG", rownames(airway)),]
airSE <- airSE[rowSums(assay(airSE)) > 4,]
dim(airSE)
assay(airSE)[1:4,1:4]

## ----norm-ma---------------------------------------------------------------
allSE <- normalize(allSE, norm.method = "quantile")

## ----plot-norm, fig.width=12, fig.height=6---------------------------------
par(mfrow=c(1,2))
boxplot(assay(allSE, "raw"))
boxplot(assay(allSE, "norm"))

## ----norm-rseq-------------------------------------------------------------
airSE <- normalize(airSE, norm.method = "quantile")

## ----sample-groups-ALL-----------------------------------------------------
allSE$GROUP <- ifelse(allSE$mol.biol == "BCR/ABL", 1, 0)
table(allSE$GROUP)

## ----sample-groups-airway--------------------------------------------------
airSE$GROUP <- ifelse(airway$dex == "trt", 1, 0)
table(airSE$GROUP)

## ----sample-blocks---------------------------------------------------------
airSE$BLOCK <- airway$cell
table(airSE$BLOCK)

## ----DE-ana-ALL------------------------------------------------------------
allSE <- deAna(allSE, padj.method = "BH")
rowData(allSE)

## ----plot-DE, fig.width=12, fig.height=6-----------------------------------
par(mfrow = c(1,2))
pdistr(rowData(allSE)$PVAL)
volcano(rowData(allSE)$FC, rowData(allSE)$ADJ.PVAL)

## ----DE-exmpl--------------------------------------------------------------
ind.min <- which.min(rowData(allSE)$ADJ.PVAL)
rowData(allSE)[ind.min,]

## ----DE-ana-airway---------------------------------------------------------
airSE <- deAna(airSE, de.method = "edgeR")
rowData(airSE)

## ----idmap-idtypes---------------------------------------------------------
idTypes("hsa")

## ----idmap-airway----------------------------------------------------------
head(rownames(airSE))
airSE <- idMap(airSE, org = "hsa", from = "ENSEMBL", to = "ENTREZID")
head(rownames(airSE))

## ----get-kegg-gs, eval=FALSE-----------------------------------------------
#  kegg.gs <- getGenesets(org = "hsa", db = "kegg")

## ----get-go-gs, eval=FALSE-------------------------------------------------
#  go.gs <- getGenesets(org = "hsa", db = "go", onto = "BP", mode = "GO.db")

## ----parseGMT--------------------------------------------------------------
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- getGenesets(gmt.file)
length(hsa.gs)
hsa.gs[1:2]

## ----sbeaMethods-----------------------------------------------------------
sbeaMethods()

## ----sbea------------------------------------------------------------------
sbea.res <- sbea(method = "ora", se = allSE, gs = hsa.gs, perm = 0, alpha = 0.1)
gsRanking(sbea.res)

## ----vst-------------------------------------------------------------------
airSE <- normalize(airSE, norm.method = "vst")

## ----gsea-rseq, eval=FALSE-------------------------------------------------
#  air.res <- sbea(method = "gsea", se = airSE, gs = hsa.gs)
#  gsRanking(sbea.res)

## ----fullrank--------------------------------------------------------------
gsRanking(sbea.res, signif.only = FALSE)

## ----eaBrowse, eval=FALSE--------------------------------------------------
#  eaBrowse(sbea.res)

## ----compile-grn-----------------------------------------------------------
hsa.grn <- compileGRN(org="hsa", db="kegg")
head(hsa.grn)

## ----nbeaMethods-----------------------------------------------------------
nbeaMethods()

## ----nbea------------------------------------------------------------------
nbea.res <- nbea(method="ggea", se=allSE, gs=hsa.gs, grn=hsa.grn)
gsRanking(nbea.res)

## ----ggea-graph, fig.width=12, fig.height=6--------------------------------
par(mfrow=c(1,2))
ggeaGraph(
    gs=hsa.gs[["hsa05217_Basal_cell_carcinoma"]], 
    grn=hsa.grn, se=allSE)
ggeaGraphLegend()

## ----dummySBEA-------------------------------------------------------------
dummySBEA <- function(se, gs)
{
        sig.ps <- sample(seq(0, 0.05, length = 1000), 5)
        insig.ps <- sample(seq(0.1, 1, length = 1000), length(gs) - 5)
        ps <- sample(c(sig.ps, insig.ps), length(gs))
        names(ps) <- names(gs)
        return(ps)
}

## ----sbea2-----------------------------------------------------------------
sbea.res2 <- sbea(method = dummySBEA, se = allSE, gs = hsa.gs)
gsRanking(sbea.res2)

## ----combine---------------------------------------------------------------
res.list <- list(sbea.res, nbea.res)
comb.res <- combResults(res.list)

## ----browse-comb, eval=FALSE-----------------------------------------------
#  eaBrowse(comb.res, graph.view=hsa.grn, nr.show=5)

## ----all-in-one, eval=FALSE------------------------------------------------
#  ebrowser(   meth=c("ora", "ggea"),
#          exprs=exprs.file, cdat=cdat.file, rdat=rdat.file,
#          org="hsa", gs=hsa.gs, grn=hsa.grn, comb=TRUE, nr.show=5)

## ----config-set------------------------------------------------------------
configEBrowser(key="OUTDIR.DEFAULT", value="/my/out/dir")

## ----config-get------------------------------------------------------------
configEBrowser("OUTDIR.DEFAULT")

## ----config-man, eval=FALSE------------------------------------------------
#  ?configEBrowser

## ----deTbl-----------------------------------------------------------------
deTable <-
     matrix(c(28, 142, 501, 12000),
            nrow = 2,
            dimnames = list(c("DE", "Not.DE"),
                            c("In.gene.set", "Not.in.gene.set")))
deTable

## ----fisher----------------------------------------------------------------
fisher.test(deTable, alternative = "greater")

