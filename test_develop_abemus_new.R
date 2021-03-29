library( devtools )
devtools::install_github("cibiobcg/abemus", build_vignettes = FALSE, ref = "develop")

library(abemus)
library(data.table)
library(parallel)
library(Matrix)
library(tidyverse)
library(gridExtra)

# sample.info.file <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/sif.txt"
# pacbam <- "/BCGLAB/ncasiraghi/abemus_test_develop_branch/pacbam"

sample.info.file <- "/Users/ncasiraghi/Documents/abemus_datasets/test_dataset_Snv0.01/test_sif_new.tsv"
pacbam <- "/Users/ncasiraghi/Documents/abemus_datasets/test_dataset_Snv0.01/pacbam_data_new"

chroms <- c(1:22,'X','Y')

read.sif <- function(file){
  sif <- list(ctrl=NA)
  tsv <- read.delim(file,header = FALSE,sep = '\t',stringsAsFactors = FALSE)
  if(ncol(tsv) == 1){
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- sif$ctrl
  } else{
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- gsub(basename(tsv[,2]),pattern = '\\.bam$',replacement = '')
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
  id <- paste0(unique(na.omit(sif$ctrl)),paste0('_chr',chr),'.pileup')
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

vaf <- as.vector(unlist(lapply(cov_vaf_by_chrom,`[`,2)))

vafth <- data.frame(spec=spec,
                    th=as.numeric(quantile(vaf,probs = spec,na.rm = TRUE)),
                    stringsAsFactors = FALSE)

# vaf threshold cov based

cov <- as.vector(unlist(lapply(cov_vaf_by_chrom,`[`,1)))

covbin <- cut(cov,
              breaks=seq(min.cov,max(cov,na.rm = TRUE),by=bin.cov),
              include.lowest=TRUE)

bin_vafth <- function(bin,vaf,covbin,spec,replica=10){

  w <- vaf[which(covbin == bin)]

  vaf.thbin <- data.frame(bin=bin,
                          spec=spec,
                          th=as.numeric(quantile(w,probs = spec,na.rm = TRUE)),
                          n=length(w[!is.na(w)]),
                          n_vaf_gtz=length(w[w > 0]),
                          n_vaf_etz=length(w[w == 0]),
                          stringsAsFactors = FALSE)

  if(FALSE){
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

  return(vaf.thbin)
}

vafth_by_bin <- do.call(rbind,lapply(levels(covbin),bin_vafth,vaf,covbin,spec))

tp <- vafth_by_bin %>%
  filter(spec == 0.995) %>%
  mutate(bin = factor(bin,levels = bin)) %>%
  gather(.,class_vaf,n_vaf,n_vaf_gtz:n_vaf_etz) %>%
  mutate(class_vaf = ifelse(class_vaf == "n_vaf_gtz", "VAF > 0", "VAF = 0"))

p1 <- ggplot(tp, aes(x=as.factor(bin), y=n_vaf,fill=class_vaf)) +
  geom_bar(stat="identity") + ylab('VAF count') + xlab('coverage bin') +
  theme(legend.position = 'top',legend.title = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 <- ggplot(tp %>% select(bin,th) %>% distinct(), aes(x=as.factor(bin), y=th)) +
  geom_bar(stat="identity") + ylab('VAF threshold') + xlab('coverage bin') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1)

# [ call snvs ]
# mincovgerm = 10
# maxafgerm = 0.2

vaf.th = vafth_by_bin %>% select(bin,spec,th)
# vaf.th = vafth
det.spec = 0.995
min.cov = 10
min.alt = 1
min.cov.ctrl = 10
max.vaf.ctrl = 0.2

callsnvs_tmp <- function(chr,idx,pacbam,vaf.th,pos,pbem_by_chrom,
                         det.spec=0.995,min.cov=10,min.alt=1,
                         mat_vaf=mat_vaf,mat_cov=mat_cov,min.cov.ctrl=10,max.vaf.ctrl=0.2,nThread=4){
  id <- paste0(sif$case[idx],paste0('_chr',chr),'.pabs')
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

    # filter based on matched ctrl [if available]-------------------------------

    # index for chr
    j <- which(chr==c(1:22,'X','Y'))

    if(!is.na(sif$ctrl[idx])){

      # index for pos
      pos.index <- which(pos[[j]] %in% df$pos)

      # index for ctrl
      ctrl.index <- which(unique(sif$ctrl) == sif$ctrl[idx])

      # add info to snvs table

      snvs <- snvs %>%
        mutate(vaf.ctrl = as.numeric(mat_vaf[[j]][pos.index,ctrl.index]),
               cov.ctrl = as.numeric(mat_cov[[j]][pos.index,ctrl.index])) %>%
        filter(cov.ctrl >= min.cov.ctrl,
               vaf.ctrl <= max.vaf.ctrl)

    }

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

snvs_list <-list()
for(idx in seq_len(length(sif$case))){
  message(sif$case[idx])
  snvs_list[[sif$case[idx]]] <- do.call(rbind,mclapply(chroms,callsnvs_tmp,
                                                       mat_vaf=mat_vaf,
                                                       mat_cov=mat_cov,
                                                       idx=idx,
                                                       pacbam=pacbam,
                                                       vaf.th=vaf.th,
                                                       pos=pos,
                                                       pbem_by_chrom=pbem_by_chrom,
                                                       mc.cores = 4))
}

