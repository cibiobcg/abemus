filter = function(i,
                  chromosomes,
                  caseout_folder=caseout_folder,
                  plasma.folder,
                  germline.folder,
                  name.patient,
                  name.plasma,
                  name.germline,
                  pbem_dir,
                  out1,
                  out2,
                  out3,
                  mincov,
                  minalt,
                  mincovgerm,
                  maxafgerm,
                  AFbycov,
                  minaf_cov_corrected,
                  xbg,
                  tab_cov_pbem,
                  afs,
                  covs,
                  coverage_binning,
                  njobs = 1){
  chrom <- chromosomes[i]
  chrom <- gsub(chrom,pattern = 'chr',replacement = '')
  # create chromsome sub-folder
  chromdir = file.path(caseout_folder,chrom)
  dir.create(chromdir,showWarnings = TRUE)
  setwd(chromdir)
  # import files
  plasma_snvs = list.files(file.path(plasma.folder,"snvs"),pattern = paste0("_chr",chrom,".pabs"),full.names = TRUE)
  n.rows.plasma_snvs = as.numeric(unlist(strsplit(trimws(x = system(paste("wc -l",plasma_snvs),intern = TRUE),which = "left"),split = " "))[[1]])
  if(n.rows.plasma_snvs == 1){
    return()
  }
  snvs = fread(plasma_snvs,stringsAsFactors = FALSE,showProgress = FALSE,header = FALSE,skip = 1,na.strings = "",colClasses = list(character=3,4,15),verbose = FALSE)
  snvs = unique(snvs)
  snvs = data.frame(snvs)
  names(snvs)=c("chr","pos","ref","alt","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp")
  if(nrow(snvs)==0){
    return()
  }
  # F1) Custom basic filters [ in plasma/tumor ]
  out <-  mclapply(seq(1,nrow(snvs),1),
                   CheckAltReads,
                   snvs=snvs,
                   mc.cores = njobs)
  snvs <- fromListToDF(out)
  snvs <- snvs[which(snvs$af > 0 ),,drop=FALSE]
  snvs <- snvs[which(snvs$cov >= mincov ),,drop=FALSE]
  snvs <- snvs[which(snvs$cov.alt >= minalt ),,drop=FALSE]
  snvs <- unique(snvs)
  if(nrow(snvs)==0){
    return()
  }

  if(is.na(germline.folder)){

    ctrl.pileup <- data.frame(chr=chrom,
                              pos=snvs$pos,
                              ref=snvs$ref,stringsAsFactors = FALSE)
    ctrl.pileup <- cbind(ctrl.pileup,NA,NA,NA,NA,NA,NA)
    names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov")

  } else {
    # print filtered positions and grep these pos only from pileup file of germline sample
    cat(unique(snvs$pos),sep = "\n",file = file.path(chromdir,"postogrep.txt"),append = FALSE)
    controlfolder_pileup <- list.files(file.path(germline.folder,"pileup"),pattern = paste0("_chr",chrom,".pileup"),full.names = TRUE)
    cmd = paste("awk -F'\t' '{if (FILENAME == \"postogrep.txt\") { t[$1] = 1; } else { if (t[$2]) { print }}}' postogrep.txt",controlfolder_pileup,"> filtered.germline.pileup.txt")
    system(cmd)
    ctrl.pileup = fread("filtered.germline.pileup.txt",stringsAsFactors = FALSE,showProgress = TRUE,header = FALSE,na.strings = "",colClasses = list(character=10))
    system("rm postogrep.txt filtered.germline.pileup.txt")
    ctrl.pileup = ctrl.pileup[,1:9]
    ctrl.pileup = unique(ctrl.pileup)
    ctrl.pileup = data.frame(ctrl.pileup)
    names(ctrl.pileup)=c("chr","pos","ref","A","C","G","T","af","cov")
  }

  # F1) Custom basic filters [ in germline ]
  common = merge(x = snvs,y = ctrl.pileup,by = c("chr","pos","ref"),all.x = TRUE,suffixes = c("_case","_control"))

  if(is.na(germline.folder)){
    toremove <- NULL
  } else {
    toremove = which(common$cov_control < mincovgerm | common$af_control > maxafgerm )
  }

  if(length(toremove)>0){
    putsnvs <- common[-toremove,,drop=FALSE]
  } else {
    putsnvs <- common
  }

  if(nrow(putsnvs) > 0){
    # F2) Filters on Variant Allelic Fraction and add pbem [ in plasma/tumor ]
    # import pbem of this chrom
    tabpbem_file = list.files(pbem_dir, pattern = paste0('bperr_',chrom,'.tsv'),full.names = TRUE)
    tabpbem = fread(input = tabpbem_file,stringsAsFactors = FALSE,showProgress = FALSE,header = FALSE,colClasses = list(character=2,character=5),data.table = FALSE)
    colnames(tabpbem) <- c("group","chr","pos","ref","dbsnp","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available",'n_pos_af_lth','n_pos_af_gth','count.A_af_gth','count.C_af_gth','count.G_af_gth','count.T_af_gth',"bperr","tot_reads_supporting_alt")
    # TABLE 1
    chrpmF1 = apply_AF_filters(chrpmF1=putsnvs,
                               AFbycov=AFbycov,
                               mybreaks=define_cov_bins(coverage_binning)[[1]],
                               af.threshold.table=minaf_cov_corrected,
                               minaf=NA,
                               mc.cores=njobs)
    chrpmF1 = chrpmF1[,c("chr","pos","ref","alt","A_case","C_case","G_case","T_case","af_case","cov_case","Ars","Crs","Grs","Trs","rev.ref","fwd.ref","cov.alt","rev.alt","fwd.alt","strandbias","A_control","C_control","G_control","T_control","af_control","cov_control","af_threshold")]
    chrpmF1$group <- paste(chrpmF1$chr,chrpmF1$pos,chrpmF1$ref,sep = ":")
    cpmf1 = merge(x = chrpmF1,y = tabpbem,by = c("group","chr","pos","ref"),all.x = TRUE)
    cpmf1 = cpmf1[,c("group","chr","pos","ref","dbsnp","alt","A_case","C_case","G_case","T_case","af_case","cov_case","Ars","Crs","Grs","Trs","rev.ref","fwd.ref","cov.alt","rev.alt","fwd.alt",
                     "strandbias","A_control","C_control","G_control","T_control","af_control","cov_control","af_threshold","tot_coverage","total.A","total.C","total.G","total.T",
                     "n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth","bperr","tot_reads_supporting_alt")]

    # Add sample/patient IDs
    cpmf1 = add_names(pm = cpmf1,
                      name.patient = name.patient,
                      name.plasma = name.plasma,
                      name.germline = name.germline)

    # compute pbem allele
    Nids = which(cpmf1$alt=='N')
    if(length(Nids)>0){cpmf1 = cpmf1[-Nids,]}
    if(!is.na(germline.folder)){cpmf1 = cpmf1[which(!is.na(cpmf1$cov_control)),,drop=FALSE]}

    if(nrow(cpmf1)==0){
      return()
    }
    out = mclapply(seq(1,nrow(cpmf1),1),
                   compute_pbem_allele,
                   abemus=cpmf1,
                   mc.cores = njobs)
    cpmf1 = fromListToDF(out)

    # add CLASS standard
    if(is.na(germline.folder)){
      cpmf1$CLASS <- NA
    } else {
      cpmf1 <- add_class(pmtab = cpmf1)
    }

    # TABLE 2
    cpmf1$af_threshold[which(is.na(cpmf1$af_threshold))] <- -1
    cpmf2 = cpmf1[which(cpmf1$af_case >= cpmf1$af_threshold),,drop=FALSE]
    if(nrow(cpmf2)==0){
      return()
    }

    # TABLE 3
    # add CLASS background pbem
    if(nrow(cpmf2)>0 & is.na(germline.folder)){
      cpmf3 <- cpmf2
      cpmf3$CLASS.xbg <- NA
    } else {
      cpmf3 <- add_class_xbg(pmtab = cpmf2,xbg = xbg)
    }

    cpmf3$bperr[which(cpmf3$bperr > 0.2)] <- 0.2
    cpmf3$bperr[which(is.na(cpmf3$bperr))] <- 0.2 # assign the highest pbem if it is NA

    covs[which.max(covs)] <- Inf # define as max coverage to avoid NA values

    if(nrow(cpmf3)>0){
      pbem_coverage_filter = sapply(1:nrow(cpmf3), function(k) tab_cov_pbem[min(which(covs>=cpmf3$cov_case[k])),min(which(afs>=cpmf3$bperr[k]))])
      cpmf3$filter.pbem_coverage <- pbem_coverage_filter
      cpmf3$pass.filter.pbem_coverage = 0
      cpmf3$pass.filter.pbem_coverage[which(cpmf3$af_case >= cpmf3$filter.pbem_coverage)] = 1
    }

    # Return chromosome tables
    write.table(cpmf1,file = 'chrpm_f1.tsv',sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
    write.table(cpmf2,file = 'chrpm_f2.tsv',sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
    write.table(cpmf3,file = 'chrpm_f3.tsv',sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
    cat(paste(colnames(cpmf1),collapse='\t'),file = file.path(caseout_folder,out1),sep = '\n')
    cat(paste(colnames(cpmf2),collapse='\t'),file = file.path(caseout_folder,out2),sep = '\n')
    cat(paste(colnames(cpmf3),collapse='\t'),file = file.path(caseout_folder,out3),sep = '\n')
  } else {
    return()
  }
}
