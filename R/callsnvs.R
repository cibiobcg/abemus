callsnvs <- function(chr,idx,pacbam,vaf.th,pos,pbem_by_chrom,
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
