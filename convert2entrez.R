#pcaMethods, STATegRa, ropls, biosigner, omicade4 fdrtool BioMark biosigner caret klaR KODAMA pathway PCA 
#MultiVarSel OmicsPLS RankProd IntLim(omic integration) DiffCor RedeR BioNeyStat pwOmics mseapca tmod
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6835268/table/metabolites-09-00200-t006/?report=objectonly
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6835268/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6835268/table/metabolites-09-00200-t008/?report=objectonly
#inmput KNN based imputation for microarray MetaboAnalystR MetaDBparse	CompoundDb bioDB RaMP graphite
#LipidMatch 500000 lipid 60 lipid type omu kegg compound and keggrest extract  gene mzmatch hmdb  and kegg
#this is to convert a affymetrix probe set to entrez id and return a matrix for gage/pathview application.
#exp: expression profile with the first row as the probe set ids. log2 transformaiton will be applied to the returned matrix
#chip:gene chip whose annotation is provided by bioconductor annotation packages:
# entering hgu133plus2 for hgu133plus2.db; hgu133a for hgu133a.db; entering hgu219 for hgu219.db
convert2entrez<-function (exp, chip){
  # read expression file to a frame, from GoldenEagle
  if (grepl("csv$", exp, ignore.case =T)){
    f<-read.table(exp, sep=",", header=T, strip.white=T)
  }else{
    f<-read.table(exp, sep="\t", header=T, strip.white=T)
  }
  
  # convert to a matrix
  m=as.matrix(f[,2:ncol(f)])
  # set row names
  row.names(m)=f[,1]
  # convert to log2 valuesls()
   m=log2(m)
  #load libraries
  require('AnnotationDbi')
  require("BiocGenerics")
  require("parallel")
  ### convert affy probe sets to Entrez IDs
  # load the chip annotaton database such as hgu133plus2.db and hgu219.db
  cmd=paste('library(',chip,'.db)',sep='')  
  eval(parse(text=cmd))
  # get affyids from user's matrix
  affyid=rownames(m)
  # get entrez ids and store them into egids2
  cmd=paste('egids2=', chip,'ENTREZID[affyid]',sep='')
  #execute command and get the entrez ids
  eval(parse(text=cmd))
  #store mapping info into annots which is workable by R directly
  annots=toTable(egids2)
  # remove records which don't have any entrez ids mapped
  m=m[annots$probe_id,]
  ### select porbe set which has the maximal IQR if multiple probe sets map to the same gene
  iqrs=apply(m,1,IQR)
  sel.rn=tapply(1:nrow(annots),annots$gene_id, function(x){x[which.max(iqrs[x])]})
  #get matrix with the maximum iqrs
  m_egid=m[sel.rn,]
  # switch to/assign entrez id to row names
  rownames(m_egid)=names(sel.rn)
  return (m_egid)
}
