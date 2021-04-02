library( devtools )
devtools::install_github("cibiobcg/abemus", build_vignettes = FALSE, ref = "develop")

library(abemus)
library(data.table)
library(parallel)
library(Matrix)
library(tidyverse)
library(gridExtra)

samples.info.file <- "/Users/nicola_casiraghi/abemus_tests/dataset_dave/sif.tsv"
pacbam <- "/Users/nicola_casiraghi/abemus_tests/dataset_dave/pacbam_data_bychrom/"

# workflow

sif <- read.sif(samples.info.file)

ref <- lapply(c(1:22,'X','Y'),getLoci,pacbam,select = 'ref')
pos <- lapply(c(1:22,'X','Y'),getLoci,pacbam,select = 'pos')
rsid<- lapply(c(1:22,'X','Y'),getLoci,pacbam,select = 'rsid') %>% lapply(function(x) which(!is.na(x)))

mat_vaf <- lapply(c(1:22,'X','Y'),pileup2mat,sif,pacbam,select='af',sparse=TRUE)

mat_cov <- lapply(c(1:22,'X','Y'),pileup2mat,sif,pacbam,select='cov',sparse=FALSE)

mat_cov_base <- list()
for(base in c('A','C','G','T')){
  mat_cov_base[[base]] <- lapply(chroms,pileup2mat,sif,pacbam,select=base,sparse=TRUE)
}

# [ compute pbem ]
pbem_by_chrom <- lapply(seq_len(length(ref)),chrom_pbem,mat_cov,mat_cov_base,ref,max.vaf.pbem = 0.2,cov.min.pbem = 10)

# [ compute vaf thresholds ]

cov_vaf_by_chrom <- lapply(seq_len(length(ref)),chrom_cov_vaf,mat_cov,mat_vaf,rsid,max.vaf=0.2,min.cov=10)

# vaf threshold not cov based

vaf <- as.vector(unlist(lapply(cov_vaf_by_chrom,`[`,2)))

spec <- seq(0.9,1,0.001)
vafth <- data.frame(spec=spec,
                    th=as.numeric(quantile(vaf,probs = spec,na.rm = TRUE)),
                    stringsAsFactors = FALSE)

# vaf threshold cov based

cov <- as.vector(unlist(lapply(cov_vaf_by_chrom,`[`,1)))

bin.cov <- 50
covbin <- cut(cov,
              breaks=seq(min.cov,max(cov,na.rm = TRUE),by=bin.cov),
              include.lowest=TRUE)

spec <- seq(0.9,1,0.001)
vafth_by_bin <- do.call(rbind,lapply(levels(covbin),bin_vafth,vaf,covbin,spec))

# tmp plots
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

vaf.th = vafth_by_bin %>% select(bin,spec,th)

snvs_list <-list()
for(idx in seq_len(length(sif$case))){
  message(sif$case[idx])
  snvs_list[[sif$case[idx]]] <- do.call(rbind,mclapply(chroms,callsnvs,
                                                       mat_vaf=mat_vaf,
                                                       mat_cov=mat_cov,
                                                       idx=idx,
                                                       pacbam=pacbam,
                                                       vaf.th=vaf.th,
                                                       pos=pos,
                                                       pbem_by_chrom=pbem_by_chrom,
                                                       mc.cores = 4))
}

