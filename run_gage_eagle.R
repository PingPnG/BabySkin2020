#this function generates gage results against go.bp, go.cc, go.mf, kegg based on users choices
#and create a ouput file for each gene dataset
#gss_signal: a tab or comma delimited gene expression file located on GoldenEagle which only contains chips passed QC and data are normalized 
#gss_meta: a tab or comma delimited meta file for gss_signal which contains treatment, QC, chip name information
#idx_chip_meta: position of column for chip names in the meta file, typically, it is 2. index starting from 1
#ctl_name: the name for control for naming the output file
#trt_name: the name for treatment for naming the output file
#ctl_meta: a named list for parameters to define control chip-ids in the gss meta file from GoldenEagle
#trt_meta: a named list for parameters to define treatment chip-ids in the gss meta file in GoldenEagle
#chip:gene chip whose annotation is provided by bioconductor annotation packages:
# 1. hgu133plus2 for hgu133plus2.db
# 2. hgu133a for hgu133a.db
# 3. hgu219 for hgu219.db
#mode:unpaired or 1ongroup

### note: meta files are very different regarding to columns to be used for selection of chip names used in gss_signal file, so, modify some codes
###       in the function to make sure your final index is right. Current paramers to be modified are: extra_value, bad, sufix

run_gage_eagle<-function (gss_signal, gss_meta, idx_chip_meta, ctl_name, trt_name, ctl_meta, trt_meta, chip, mode){
  # check QC_PASS == 1 is passed in ctl_meta and trt_meta list
  if (!('QC_PASS' %in% names(ctl_meta))){
    stop('QC_PASS must be included in ctl_meta list') 
  }else if (ctl_meta[['QC_PASS']]!=1){
    stop('QC_PASS value in ctl_meta must be 1')
  }
  if (!('QC_PASS' %in% names(trt_meta))){
    stop ('QC_PASS must be included in trt_meta list')
  }else if(trt_meta[['QC_PASS']]!=1){
    stop('QC_PASS in trt_meta must be 1')
  }
  
  #read in the input expression file [exp] and covert it to a matrix for gage
  if(!exists("convert2entrez")) source("~/R_script/convert2entrez.R") 
  # read in gss expression data
  exp_m<-convert2entrez(gss_signal, chip)
  print ("converting expression to matrix...done")
  # determin the separater in gss_meta file from its header line
  f <- file(gss_meta, 'r')
  header = readLines(f, n=1)
  close(f)
  
  # read in gss_meta file
  if (grepl(',', header, fixed=TRUE)){
    meta<-read.table(gss_meta, sep=",", header=T, strip.white=T,fileEncoding="latin1",skipNul = T)
  }else{
    meta<-read.table(gss_meta, sep="\t", header=T, strip.white=T,fileEncoding="latin1",skipNul = T)
  }
  #print (head(meta))
  
  #select headers for control group in gss_signal file, those headers are stored in the second column of gss_meta
  #note: QC_PASS must be 1 for those passed QC chips
  ctl_which=paste (paste0('meta$',names(ctl_meta)), unlist(ctl_meta, use.names = F), sep = '=="', collapse = '" & ')
  # add last "
  ctl_which = paste0(ctl_which,'"')
  ctl_query = paste('meta[which(', ctl_which, '), idx_chip_meta]')
  ctl_headers<-eval(parse(text = ctl_query))
  #print(ctl_headers)
  
  #select headers for treatment group in gss_signal file, those headers are stored in the second column of gss_meta
  trt_which = paste(paste0('meta$', names(trt_meta)), unlist(trt_meta, use.names = F), sep='=="', collapse = '" &')
  # add the last "
  trt_which = paste0(trt_which, '"')
  trt_query = paste('meta[which(', trt_which, '), idx_chip_meta]')
  trt_headers<-eval(parse(text= trt_query))
  #print(trt_headers)
 
   #get chip names from gss_signal
  headers<-colnames(exp_m)
  
  ##### note: chip names from mega file sometimes don't match the column names in gss_signal file, so, code below need to be 
  #####       adjusted for each study
  #### sufix is the extra string used in gss_signal file which needs to be added to headers stored in gss_meta file
  sufix<-""
  ctl<-paste(ctl_headers,sufix,sep="")
  trt<-paste(trt_headers,sufix,sep="")
  #print(trt)
  idx_ctl=which(headers %in% ctl)
  idx_trt=which(headers %in% trt)
  print("treatment group chip index:")
  print (idx_trt)
  print("control group chip index:")
  print(idx_ctl)
  
   # for gage analysis
  BiocManager::install("gage")
  library(gage)
  ### build most recent or current human GO datasets
  if(!exists('go.hs')){
    go.hs <<- go.gsets(species='human', id.type = 'eg')
  }
  if (!exists("go.mf")){
    go.mf <<- go.hs$go.sets[go.hs$go.subs$MF]
  }
  if (!exists("go.cc")){
    go.cc <<- go.hs$go.sets[go.hs$go.subs$CC]

  }
  if (!exists("go.bp")){
    go.bp <<- go.hs$go.sets[go.hs$go.subs$BP]
  }
  #load KEGG metabolic and signaling pathway sets from gage directly
  if(!exists("kegg.gs")){
    kegg.gs <<- data("kegg.gs")
  }
  
  #assign gene sets to score
  gs=c("go.bp", "go.cc", "go.mf","kegg.gs")
  #gs=c("dej.gs")
  #find the GSS number
  gss<-substr(gss_signal,regexpr("GSS\\d*_signal",gss_signal)[[1]],regexpr("_signal",gss_signal)[[1]]-1)
  #run gage against gene sets
  for (g in gs){
    rst=gage(exp_m,gsets=get(g),ref=idx_ctl, samp=idx_trt, compare=mode)
    outF = paste(trt_name,"vs",ctl_name,"." ,gss, ".",g,".",mode,".txt",sep="")
    write.table(rbind(rst$greater,rst$less), file=outF,sep="\t") 
    #write.table(rst$stats, file=paste(key_trt,".",exp, ".",g,".",mode,".", stats.txt",sep=""),sep="\t")  
  }
  
}
