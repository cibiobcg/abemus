library( devtools )
devtools::install_github("cibiobcg/abemus", build_vignettes = F)
library(abemus)

library(data.table)
library(parallel)
library(Matrix)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

sample.info.file <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/sif.txt"
pacbam <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/pacbam"

chroms <- c(1:22,'X','Y')

read.sif <- function(file){
  sif <- list(ctrl=NA)
  tsv <- read.delim(file,header = FALSE,sep = '\t',stringsAsFactors = FALSE)
  if(ncol(tsv) == 1){
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- sif$ctrl
  } else{
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- gsub(unique(basename(tsv[,2])),pattern = '\\.bam$',replacement = '')
  }
  return(sif)
}

get_loci <- function(chr,pacbam,select='ref',nThread=4){
  p <- list.files(pacbam,pattern = paste0('_chr',chr,'.pileup$'),full.names = TRUE)[1]
  if (is.na(p)) {
    return(NA)
  } else{
    df <- fread(input = p,
                sep = '\t',
                stringsAsFactors = FALSE,
                header = TRUE,
                select = select,
                verbose = FALSE,
                data.table = FALSE,
                na.strings = '',
                nThread = nThread)
    return(as.vector(unlist(df)))
  }
}

pileup2mat <- function(chr,sif,pacbam,select,sparse=FALSE,nThread=4){
  id <- paste0(sif$ctrl,paste0('_chr',chr),'.pileup')
  ff <- file.path(pacbam,id)
  ff <- ff[which(file.exists(ff))]
  if(length(ff) == 0){
    return(NA)
  } else{
    df <- lapply(ff,
                 fread,
                 sep = '\t',
                 stringsAsFactors = FALSE,
                 header = TRUE,
                 select = select,
                 verbose = FALSE,
                 data.table = FALSE,
                 na.strings = '',
                 nThread = nThread)
    if(sparse){
      return(Matrix(as.matrix(do.call(cbind,df)), sparse = TRUE))
    } else {
      return(as.matrix(do.call(cbind,df)))
    }
  }
}

# workflow

sif <- read.sif(sample.info.file)

ref <- lapply(chroms,get_loci,pacbam,select = 'ref')
pos <- lapply(chroms,get_loci,pacbam,select = 'pos')

mat_vaf <- lapply(chroms,pileup2mat,sif,pacbam,select='af',sparse=TRUE)

mat_cov <- lapply(chroms,pileup2mat,sif,pacbam,select='cov',sparse=FALSE)

mat_cov_base <- list()
for(base in c('A','C','G','T')){
  mat_cov_base[[base]] <- lapply(chroms,pileup2mat,sif,pacbam,select=base,sparse=TRUE)
}

# [ compute pbem ]
max.vaf.pbem = 0.2
cov.min.pbem = 10

chrom_pbem <- function(i,mat_cov,mat_cov_base,ref){

  if(all(is.na(ref[[i]]))){
    return(NA)
  }

  x <- rep(NA,length(ref[[i]]))

  mat_cov_chr <- mat_cov[[i]]

  filtout <- unique(rbind(which(mat_cov_chr <= cov.min.pbem,arr.ind = TRUE),
                          which(mat_vaf[[i]] >= max.vaf.pbem,arr.ind = TRUE)))

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

# [ compute af thresholds ]
spec <- seq(0.9,1,0.0001)
max.vaf <- 0.2
min.cov <- 10
bin.cov <- 50

chrom_cov_vaf <- function(i,mat_cov,mat_vaf){

  if(all(is.na(mat_cov[[i]]))){
    return(list(NA,NA))
  }

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

vafth <- data.frame(spec=spec,
                    th=as.numeric(quantile(vaf,probs = spec,na.rm = TRUE)),
                    stringsAsFactors = FALSE)

# vaf threshold cov based

cov <- as.vector(unlist(ext(cov_vaf_by_chrom,1)))

covbin <- cut(cov,
              breaks=seq(min.cov,max(cov,na.rm = TRUE),by=bin.cov),
              include.lowest=TRUE)

# barplot(table(covbin),las=2,cex.names = 0.4)

bin_vafth <- function(bin,vaf,covbin,spec){

  vaf.thbin <- data.frame(bin=bin,
                          spec=spec,
                          th=as.numeric(quantile(vaf[which(covbin == bin)],probs = spec,na.rm = TRUE)),
                          stringsAsFactors = FALSE)

  return(vaf.thbin)
}

vafth_by_bin <- do.call(rbind,lapply(levels(covbin),bin_vafth,vaf,covbin,spec))

# improve [!]
# vafth_by_bin$spec.bin <- cut(vafth_by_bin$spec,breaks = seq(0.9,1,0.01),include.lowest = TRUE)
# for(sb in levels(vafth_by_bin$spec.bin)){
#   df <- vafth_by_bin[which(vafth_by_bin$spec.bin == sb),]
#   hm <- ggplot(df, aes(x = bin, y = as.factor(spec), fill = th)) +
#     geom_tile() +
#     facet_wrap(~spec.bin,scales = "free") +
#     scale_fill_distiller(palette = 'RdYlBu')
#   print(hm)
# }

# [ call snvs ]
# replicas = 1000
# replicas.in.parallel = 1
# coeffvar.threshold = 0.01
# mincovgerm = 10
# maxafgerm = 0.2

vaf.th = vafth_by_bin
vaf.th = vafth
spec = 0.995

min.cov = 10
min.alt = 1

callsnvs_tmp <- function(chr,case,pacbam,vaf.th,pos,pbem_by_chrom,det.spec=0.995,min.cov=10,min.alt=1,nThread=4){
  id <- paste0(case,paste0('_chr',chr),'.pabs')
  ff <- file.path(pacbam,id)
  if(!file.exists(ff)){
    return(NA)
  } else{
    out <- list()
    df <- fread(input = ff,
                sep = '\t',
                stringsAsFactors = FALSE,
                header = TRUE,
                verbose = FALSE,
                data.table = FALSE,
                na.strings = '',
                nThread = nThread) %>% filter(cov >= min.cov)

    # filter based on cov
    snvs <- do.call(rbind,lapply(seq_len(nrow(df)),CheckAltReads,df)) %>% filter(cov.alt >= min.alt)
    out$tab1 <- snvs

    # filter based on vaf
    if(ncol(vaf.th)==2){
      snvs <- snvs %>%
        mutate(af.th = vaf.th$th[which(vafth$spec == det.spec)]) %>%
        filter(af >= af.th)
    }

    if(ncol(vaf.th)==3){

      vaf.th.spec <- vaf.th %>% filter(spec == det.spec)
      bin.cov <- diff(parse_number(str_split(vaf.th.spec$bin[1],pattern = ',',simplify = TRUE)))

      snvs$bin <- cut(snvs$cov,
                      breaks=seq(min.cov,max(snvs$cov,na.rm = TRUE),by=bin.cov),
                      include.lowest=TRUE)

      snvs <- left_join(x = snvs,y = vaf.th.spec,by = 'bin') %>%
        select(-spec) %>%
        rename(af.th = th, cov.bin = bin) %>%
        filter(af >= af.th)
    }
    out$tab2 <- snvs

    # filter based on pbem
    j <- which(chr==c(1:22,'X','Y'))
    snvs$pbem <- pbem_by_chrom[[j]][which(pos[[j]] %in% snvs$pos)]

    snvs <- snvs %>% drop_na(pbem)

    snvs$pbem.th <- sapply(seq_len(nrow(snvs)), function(k) bombanel_tab_cov_pbem[min(which(bombanel_covs>=snvs$cov[k])),min(which(bombanel_afs>=snvs$pbem[k]))])

    snvs <- snvs %>% filter(af >= pbem.th)

    out$tab3 <- snvs

    return(out)
  }
}







