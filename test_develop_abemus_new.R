library( devtools )
devtools::install_github("cibiobcg/abemus", build_vignettes = T,ref = 'develop')

library(data.table)
library(parallel)
library(Matrix)
library(abemus)

outdir <-  "/BCGLAB/ncasiraghi/abemus_test_develop_branch"
sample.info.file <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/sample_info_toy.txt"
targetbed <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/TST170_DNA_target.nochr.merged.nochromY.bed"
pacbamfolder_bychrom <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/bychrom/"

import_sif <- function(main_sif){
  df  <-  read.delim(main_sif,as.is=TRUE,stringsAsFactors = FALSE,header = FALSE)
  colnames(df) <- c("patient","plasma","plasma.bam","germline","germline.bam")
  # remove NAs and keep only unique germline samples
  df_ctrl <- df[which(!is.na(df$germline.bam)),,drop=FALSE]
  df_ctrl <- unique(df_ctrl[,c(4,5)])
  # remove NAs and keep only case samples having matched germline samples
  df_cases <- df[which(!is.na(df$plasma.bam)),,drop=FALSE]
  return(list(df_ctrl=df_ctrl,df_cases=df_cases))
}

sif <- import_sif(sample.info.file)

bed <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)

get_germlineset <- function(sifgerm,pacbamfolder_bychrom,chrom){
  germlineset = c()
  chrom <- gsub(chrom,pattern = 'chr',replacement = '')
  for(id in seq_len(nrow(sifgerm))){
    thisSample=sifgerm[id,]
    name = gsub(basename(thisSample$germline.bam),pattern = ".bam",replacement = "")
    thisPileup.file = list.files(file.path(pacbamfolder_bychrom,name,"pileup"),full.names = TRUE,pattern = paste0("_chr",chrom,"\\.pileup$"))
    germlineset = c(germlineset,thisPileup.file)
  }
  germlineset <- unique(germlineset)
  return(germlineset)
}

getInfo <- function(chrom,sif,pacbamfolder_bychrom,select){
  germlineset <- get_germlineset(sifgerm = sif$df_ctrl,
                                 pacbamfolder_bychrom = pacbamfolder_bychrom,
                                 chrom = chrom)
  df <- fread(germlineset[1],
              sep = '\t',
              stringsAsFactors = FALSE,
              skip = 1,select = select,
              verbose = FALSE,
              data.table=FALSE,
              nThread = 4)
  return(as.vector(unlist(df)))
}

pileup2mat <- function(chrom,sif,pacbamfolder_bychrom,select,sparse=FALSE){
  germlineset <- get_germlineset(sifgerm = sif$df_ctrl,
                                 pacbamfolder_bychrom = pacbamfolder_bychrom,
                                 chrom = chrom)
  df <- lapply(germlineset,fread,sep = '\t',stringsAsFactors = FALSE,skip = 1,select = select,verbose = FALSE,data.table=FALSE, nThread = 4)
  if(sparse){
    return(Matrix(as.matrix(do.call(cbind,df)), sparse = TRUE))
  } else {
    return(as.matrix(do.call(cbind,df)))
  }
}

ref <- lapply(bed$chromosomes,getInfo,sif,pacbamfolder_bychrom,select=3)
pos <- lapply(bed$chromosomes,getInfo,sif,pacbamfolder_bychrom,select=2)

mat_vaf <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=8,sparse=TRUE)
mat_cov <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=9)

dict <- 4:7
names(dict) <- c('A','C','G','T')
mat_cov_base <- list()
for(base in names(dict)){
  mat_cov_base[[base]] <- lapply(bed$chromosomes,pileup2mat,sif,pacbamfolder_bychrom,select=as.integer(dict[base]),sparse=TRUE)
}

# [ compute pbem ]
max.vaf.pbem <- af_max_to_compute_pbem <- 0.2
cov.min.pbem <- coverage_min_to_compute_pbem <- 10

chrom_pbem <- function(i,mat_cov,mat_cov_base,ref){

  x <- rep(NA,length(ref[[i]]))

  mat_cov_chr <- mat_cov[[i]]

  filtout <- unique(rbind(which(mat_cov_chr < cov.min.pbem,arr.ind = TRUE),
                          which(mat_vaf[[i]] > max.vaf.pbem,arr.ind = TRUE)))

  mat_cov_chr[filtout] <- NA

  den <- rowSums(mat_cov_chr,na.rm = TRUE)

  # hist(den)

  ext <- function(lst,n){
    sapply(lst,FUN = '[',n)
  }

  get_num <- function(alt,idx,filtout){
    alt[filtout,] <- NA
    rowSums(alt[idx,],na.rm = TRUE)
  }

  for(base in c('A','C','G','T')){

    alt_list <- ext(mat_cov_base[setdiff(c('A','C','G','T'),base)],n = i)

    idx <- which(ref[[i]] == base)

    num <- sapply(alt_list, get_num, idx, filtout)

    x[idx] <- rowSums(num)/den[idx]

  }

  return(x)
}

pbem_by_chrom <- lapply(seq_len(length(ref)),chrom_pbem,mat_cov,mat_cov_base,ref)

# [ compute af threshold ] all chroms together
# af_max_to_compute_thresholds = 0.2
# coverage_min_to_compute_thresholds = 10
# coverage_binning = 50,
# probs = seq(0.9,1,0.0001)){

spec <- seq(0.9,1,0.0001)
max.vaf <- 0.2
min.cov <- 10
bin.cov <- 50

# roll over chromosomes [!]

chrom_cov_vaf <- function(i,mat_cov,mat_vaf){

  cov <- as.vector(mat_cov[[i]])
  cov[which(cov < min.cov)] <- NA

  vaf <- as.vector(mat_vaf[[i]])
  vaf[which(vaf > max.vaf)] <- NA

  return(list(cov,vaf))
}

cov_vaf_by_chrom <- lapply(seq_len(length(ref)),chrom_cov_vaf,mat_cov,mat_vaf)

# vaf threshold not cov based

ext <- function(lst,n){
  sapply(lst,FUN = '[',n)
}

vaf <- as.vector(unlist(ext(cov_vaf_by_chrom,2)))

vaf.th <- data.frame(spec=spec,
                     th=as.numeric(quantile(vaf,probs = spec,na.rm = TRUE)),
                     stringsAsFactors = FALSE)

# vaf threshold cov based

cov <- as.vector(unlist(ext(cov_vaf_by_chrom,1)))

covbin <- cut(cov,
              breaks=seq(min.cov,max(cov,na.rm = TRUE),by=bin.cov),
              include.lowest=TRUE)

barplot(table(covbin),las=2,cex.names = 0.4)

bin_vafth <- function(bin,vaf,covbin,spec){

  vaf.thbin <- data.frame(bin=bin,
                          spec=spec,
                          th=as.numeric(quantile(vaf[which(covbin == bin)],probs = spec,na.rm = TRUE)),
                          stringsAsFactors = FALSE)

  return(vaf.thbin)
}

vafth_by_bin <- do.call(rbind,lapply(levels(covbin),bin_vafth,vaf,covbin,spec))

# hm <- ggplot(vafth_by_bin, aes(x = bin, y = as.factor(spec), fill = th)) + geom_tile()












