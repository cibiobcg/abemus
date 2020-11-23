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
library(gridExtra)
library(ggplot2)

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
rsid <- lapply(chroms,get_loci,pacbam,select = 'rsid') %>% lapply(function(x) which(!is.na(x)))

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
spec <- seq(0.9,1,0.001)
max.vaf <- 0.2
min.cov <- 10
bin.cov <- 50

chrom_cov_vaf <- function(i,mat_cov,mat_vaf,rsid){

  if(all(is.na(mat_cov[[i]]))){
    return(list(NA,NA))
  }

  cov <- as.vector(mat_cov[[i]])
  cov[which(cov < min.cov)] <- NA

  vaf <- as.vector(mat_vaf[[i]])
  vaf[which(vaf > max.vaf)] <- NA
  vaf[rsid[[i]]] <- NA

  return(list(cov,vaf))
}

cov_vaf_by_chrom <- lapply(seq_len(length(ref)),chrom_cov_vaf,mat_cov,mat_vaf,rsid)

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

bin_vafth <- function(bin,vaf,covbin,spec=0.995,replica=10){

  w <- vaf[which(covbin == bin)]

  vaf.thbin <- data.frame(bin=bin,
                     spec=spec,
                     th=as.numeric(quantile(w,probs = spec,na.rm = TRUE)),
                     keep=NA,
                     run=NA,
                     stringsAsFactors = FALSE)

  if( length(w) >= 100 ){

    nn <- length(w) - round(length(w) * seq(0.01,0.99,0.01))

    for(sz in nn){

      for(run in 1:replica){
        this <- data.frame(bin=bin,
                           spec=spec,
                           th=as.numeric(quantile(sample(w,size = sz, replace = FALSE),probs = spec,na.rm = TRUE)),
                           keep=sz,
                           run=run,
                           stringsAsFactors = FALSE)

        vaf.thbin <- rbind(vaf.thbin,this)
      }

    }

    vaf.thbin_full <- vaf.thbin %>%
      filter(is.na(keep)) %>%
      select(-keep,-run)

    vaf.thbin_subs <- vaf.thbin %>%
      filter(!is.na(keep)) %>%
      group_by(spec,keep) %>%
      summarise(cvar=sd(th,na.rm=TRUE)/mean(th,na.rm=TRUE),median=median(th,na.rm = TRUE)) %>%
      arrange(desc(keep))

    vaf.thbin <- full_join(x = vaf.thbin_full,y = vaf.thbin_subs,by='spec') %>%
      mutate(delta = abs(median - th))

    # tested only with one spec and one bin
    p <- ggplot(vaf.thbin, aes(keep, delta)) + geom_bar(stat = 'identity') +
      ggtitle(bin) +
      geom_hline(yintercept = 0.001,linetype="dashed", color = "red")
    print(p)

    return(vaf.thbin)
  }
}

vafth_by_bin <- do.call(rbind,lapply(levels(covbin),bin_vafth,vaf,covbin,spec))

barplot(vafth_by_bin$th,names.arg = vafth_by_bin$bin,las=2)

# [ adjust vafth_by_bin ]

get_afth <- function(num,next.bin.AFgtz,next.bin.AFetz,vaf.gtz){
  afth <- quantile.zaf(x = sample(vaf.gtz, size = next.bin.AFgtz, replace = FALSE),
                       probs = det.spec,
                       nz = next.bin.AFetz)
  return(data.frame(spec=det.spec,
                    th=as.numeric(afth),
                    stringsAsFactors = FALSE))
}

par(mfrow=c(2,1))
barplot(rbind(m$vaf.gtz,m$vaf.etz),beside = F)
barplot(m$th)


adjust_vafth_by_bin <- function(det.spec,vafth_by_bin,vaf,covbin,coeffvar.th=0.001,replicas=100,replicas.in.parallel=2){

  m <- vafth_by_bin %>%
    filter(spec == det.spec) %>%
    arrange(desc(vaf.gtz),desc(vaf.etz))

  stop <- nrow(m)-1
  last.stable.card <- 0
  tab <- c()

  for(i in 1:stop){
    current.bin <- m$bin[i]
    current.bin.card <- m$card[i]

    vaf.gtz <- vaf[which(covbin == current.bin)]
    vaf.gtz <- vaf.gtz[which(vaf.gtz > 0)]

    if(current.bin.card >= last.stable.card){
      for(j in (i+1):(stop+1)){
        next.bin <- m$bin[j]
        next.bin.card <- m$card[j]
        next.bin.AFetz <- m$vaf.etz[j]
        next.bin.AFgtz <- m$vaf.gtz[j]
        out <- mclapply(1:replicas,get_afth, next.bin.AFgtz, next.bin.AFetz, vaf.gtz, mc.cores = replicas.in.parallel)
        ReplicasTable <- do.call(rbind,out)
        coeffvar <- as.numeric(sd(ReplicasTable$th,na.rm=TRUE)/mean(ReplicasTable$th,na.rm=TRUE))
        if(is.na(coeffvar)){
          coeffvar <- 0
        }
        x = data.frame(current.bin=current.bin,
                       next.bin=next.bin,
                       coeffvar=coeffvar,
                       median.vafth=median(ReplicasTable$th,na.rm = TRUE),
                       last.stable.card=last.stable.card,
                       stringsAsFactors = FALSE)
        tab=rbind(tab,x)
        if(coeffvar > coeffvar.th){
          last.stable.card = max(last.stable.card,next.bin.card)
          break
        }
      }
    }
  }

  # correct the original threhsold table
  mcorr <- m
  cc <- tab$last.stable.card[nrow(tab)]

  last.afth.used <- 1
  start <- 2
  end <- length(minaf_cov_corrected)-1

  for(i in start:end){
    bin.name <- mcorr$bin[i]
    bin.card <- mcorr$card[which(mcorr$bin == bin.name)]
    if(!identical(bin.card,numeric(0))){
      if(bin.card < cc & last.afth.used == 1){
        mcorr$th[i] <- last.afth.used
      }
      if(bin.card >= cc){
        last.afth.used <- mcorr$th[i]
      }
      if(bin.card < cc & last.afth.used != 1){
        mcorr$th[i] <- last.afth.used
      }
    }
  }

  mcorr$th[which(is.na(mcorr$th))] <- last.afth.used
  mcorr$th[which(mcorr$th==1)] <- NA

  vafth_by_bin_corr <- full_join(suffix = c('.old','.corr'),x = m, y = mcorr, by=c('bin','spec','vaf.gtz','vaf.etz','card')) %>%
    select(bin,spec,th.old,th.corr,vaf.gtz,vaf.etz,card) %>%
    rename(th = th.old) %>%
    mutate(changed = if_else(th != th.corr,TRUE,FALSE))

  rownames(vafth_by_bin_corr) <- vafth_by_bin_corr$bin
  vafth_by_bin_corr <- vafth_by_bin_corr[vafth_by_bin$bin,]

  barplot(vafth_by_bin_corr$th)
  barplot(vafth_by_bin_corr$th.corr)


  # # correct the AF threhsold in bin where there is Inf as limit
  # N = length(minaf_cov_corrected)
  # minaf_cov_corrected[N] <- minaf_cov_corrected[N-1]
  #
}


# [ call snvs ]
# replicas = 1000
# replicas.in.parallel = 1
# coeffvar.threshold = 0.01
# mincovgerm = 10
# maxafgerm = 0.2

vaf.th = vafth_by_bin %>% select(-vaf.gtz,-vaf.etz)
# vaf.th = vafth
spec = 0.995
min.cov = 10
min.alt = 1

callsnvs_tmp <- function(chr,case,pacbam,vaf.th,pos,pbem_by_chrom,det.spec=0.995,min.cov=10,min.alt=1,nThread=4){
  id <- paste0(case,paste0('_chr',chr),'.pabs')
  ff <- file.path(pacbam,id)
  if(!file.exists(ff)){
    return(NA)
  } else{
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

    # filter based on pbem
    j <- which(chr==c(1:22,'X','Y'))
    snvs$pbem <- pbem_by_chrom[[j]][which(pos[[j]] %in% snvs$pos)]

    snvs <- snvs %>% drop_na(pbem)

    snvs$pbem.th <- sapply(seq_len(nrow(snvs)), function(k) bombanel_tab_cov_pbem[min(which(bombanel_covs>=snvs$cov[k])),min(which(bombanel_afs>=snvs$pbem[k]))])

    snvs <- snvs %>%
      filter(af >= pbem.th) %>%
      select(chr,pos,ref,alt,rsid,cov,cov.bin,cov.alt,strandbias,af,af.th,pbem,pbem.th) %>%
      rename(vaf = af, vaf.th = af.th)

    return(snvs)
  }
}

snvs <-list()
for(case in sif$case){
  message(case)
  snvs[[case]] <- do.call(rbind,mclapply(chroms,callsnvs_tmp,case,pacbam,vaf.th,pos,pbem_by_chrom,mc.cores = 4))
}



