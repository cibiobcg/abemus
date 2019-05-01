#' Compute pbem
#'
#' Compute per-base error model on each targeted position and save AFs by bins of coverage
#' @export
#' @param sample.info.file sample info file listing cases and controls. tab-delimeted file
#' @param targetbp folder with RData for each annotated positions
#' @param pacbamfolder folder with pileups
#' @param outdir output folder for this step analysis
#' @param outdir.bperr.name folder will be created in outdir. default: "BaseErrorModel"
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param bam.with.chr default: FALSE
#' @param af_max_to_compute_thresholds To compute AF thresholds, consider only positions with AF <= af_max_to_compute_thresholds. default 0.2
#' @param coverage_min_to_compute_thresholds To compute AF threshold, consider only positions with coverage >= coverage_min_to_compute_thresholds. default 10
#' @param af_max_to_compute_pbem To compute pbem, consider only positions with AF <= af_max_to_compute_pbem. default: 0.2
#' @param coverage_min_to_compute_pbem To compute pbem, consider only positions with coverage >= coverage_min_to_compute_pbem. default: 10
#' @param n_pos_af_th When compute pbem, count in how many germline samples the position has an AF >= n_pos_af_th. default: 0.2
#' @param mc.cores mc.core param from mclapply. default: 1
#' @param step into how many positions to split the chrom file. default: 5000
#' @return list(bperr, bperr_summary, bperr_tabstat)
compute_pbem <- function(sample.info.file,
                         targetbp,
                         pacbamfolder,
                         outdir,
                         outdir.bperr.name = "BaseErrorModel",
                         coverage_binning = 50,
                         af_max_to_compute_thresholds = 0.2,
                         coverage_min_to_compute_thresholds = 10,
                         af_max_to_compute_pbem = 0.2,
                         coverage_min_to_compute_pbem = 10,
                         n_pos_af_th = 0.2,
                         mc.cores = 1,
                         step = 5000,
                         bam.with.chr = FALSE){
  cat(paste("[",Sys.time(),"]\tReading the sample.info.file","\n"))
  sif <- import_sif(main_sif = sample.info.file)

  cat(paste("[",Sys.time(),"]\tReading chromosomes from bpcovered.tsv","\n"))
  chromosomes = read.delim(file = file.path(targetbp,"bpcovered.tsv"),as.is=T)
  chromosomes = sort(paste0("chr",chromosomes[-nrow(chromosomes),1]))

  cat(paste("[",Sys.time(),"]\tComputation of per-base error model","\n"))

  if(!file.exists(file.path(outdir, outdir.bperr.name))){
    dir.create(file.path(outdir, outdir.bperr.name), showWarnings = TRUE)
  }
  setwd(file.path(outdir, outdir.bperr.name))

  for(chrom in unique(chromosomes)){
    cat(paste("[",Sys.time(),"]\tchromosome:",chrom),"\n")

    tp = list.files(targetbp,pattern = paste0(chrom,"\\.RData$"),full.names = T)
    load(tp,verbose = F)
    chromTargetPositions = as.data.table(chromTargetPositions,keep.rownames = F)

    mytargets = unique(chromTargetPositions)
    targetsFILT <- mytargets
    targetsFILT$randompos <- 1 # deprecated, to be removed

    # Final targets
    targets = unique(targetsFILT)
    targets = targets[with(targets,order(chr,pos)),]
    targets = as.data.frame(targets)

    if(bam.with.chr){
      targets$chr = paste0('chr',targets$chr)
    }

    cat(paste("[",Sys.time(),"]\ttotal positions to check in this chromosome :",nrow(targets)),"\n")
    mclapply(seq(1,nrow(targets),step),pos2bperr,
             targets=targets,
             germlineset=get_germlineset(sifgerm = sif[[1]],pacbamfolder = pacbamfolder,chrom = chrom),
             step=step,
             chrom=chrom,
             covbin=define_cov_bins(coverage_binning)[[1]],
             lev=define_cov_bins(coverage_binning)[[2]],
             af_max_to_compute_thresholds=af_max_to_compute_thresholds,
             coverage_min_to_compute_thresholds=coverage_min_to_compute_thresholds,
             af_max_to_compute_pbem=af_max_to_compute_pbem,
             coverage_min_to_compute_pbem=coverage_min_to_compute_pbem,
             n_pos_af_th=n_pos_af_th,
             mc.cores=mc.cores)

    cmda = paste("cat *_pbem.table.txt >",paste0("bperr_",chrom,".tsv"))
    system(cmda)
    system("rm *_pbem.table.txt")
    cmdb = paste("cat *_afgtz.table.txt >",paste0("afgtz_",chrom,".tsv"))
    system(cmdb)
    system("rm *_afgtz.table.txt")
    cmdc = paste("cat *_afz.table.txt >",paste0("afz_",chrom,".tsv"))
    system(cmdc)
    system("rm *_afz.table.txt")
  }
  cmd.merge = paste("cat afgtz_chr*.tsv > afgtz.tsv")
  system(cmd.merge)
  cmd.merge = paste("cat afz_chr*.tsv > afz.tsv")
  system(cmd.merge)

  # save counter of afz as RData
  afztab = read.delim(file = "afz.tsv",sep="\t",as.is = T,header=F)
  afz = apply(afztab,2,sum)
  names(afz)=define_cov_bins(coverage_binning)[[2]]
  save(afz,file = "afz.RData",compress = T)

  # overall statistics on pbems
  cmd.merge = paste("cat bperr_chr*.tsv > bperr.tsv")
  system(cmd.merge)
  bperr = fread(file.path(outdir, outdir.bperr.name,"bperr.tsv"),stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5))

  # summary stats for pbem across the target
  bperr_summary = summary(bperr$V22)
  names = names(bperr_summary)
  bperr_summary = data.frame(as.numeric(bperr_summary))
  bperr_summary = rbind(bperr_summary,sd(x = bperr$V16,na.rm = T))
  rownames(bperr_summary) = c(names,"std")
  write.table(bperr_summary,file = file.path(outdir, outdir.bperr.name,"bperr_summary.tsv"),row.names = T,col.names = F,quote = F,sep = "\t")

  # compute background pbem
  bperr_subset = bperr[which(bperr$V17 == 0),]
  bgpbem = (sum(as.numeric(bperr_subset$V23)))/(sum(as.numeric(bperr_subset$V10)))
  mean_pbem = mean(as.numeric(bperr_subset$V22),na.rm = T)
  bperr_tabstat = data.frame(background_pbem = bgpbem,
                             mean_pbem = mean_pbem,
                             stringsAsFactors = F)
  write.table(bperr_tabstat,file = file.path(outdir, outdir.bperr.name,"pbem_background.tsv"),row.names = F,col.names = T,quote = F,sep = "\t")

  return(list(bperr=bperr,
              bperr_summary=bperr_summary,
              bperr_tabstat=bperr_tabstat))
}

