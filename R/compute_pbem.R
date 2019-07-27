#' Compute pbem
#'
#' Compute per-base error model on each targeted position and save AFs by bins of coverage
#' @export
#' @param sample.info.file sample info file listing cases and controls. tab-delimeted file
#' @param targetbed targeted regions in BED format.
#' @param pacbamfolder_bychrom folder with pileups
#' @param outdir output folder for this step analysis
#' @param outdir.bperr.name folder will be created in outdir. default: "BaseErrorModel"
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param af_max_to_compute_thresholds To compute AF thresholds, consider only positions with AF <= af_max_to_compute_thresholds. default 0.2
#' @param coverage_min_to_compute_thresholds To compute AF threshold, consider only positions with coverage >= coverage_min_to_compute_thresholds. default 10
#' @param af_max_to_compute_pbem To compute pbem, consider only positions with AF <= af_max_to_compute_pbem. default: 0.2
#' @param coverage_min_to_compute_pbem To compute pbem, consider only positions with coverage >= coverage_min_to_compute_pbem. default: 10
#' @param n_pos_af_th When compute pbem, count in how many germline samples the position has an AF >= n_pos_af_th. default: 0.2
#' @param mc.cores mc.core param from mclapply. default: 1
#' @param step into how many positions to split the chrom file. default: 5000
#' @return list(bperr, bperr_summary, bperr_tabstat)
compute_pbem <- function(sample.info.file,
                         targetbed,
                         pacbamfolder_bychrom,
                         outdir,
                         outdir.bperr.name = "BaseErrorModel",
                         coverage_binning = 50,
                         af_max_to_compute_thresholds = 0.2,
                         coverage_min_to_compute_thresholds = 10,
                         af_max_to_compute_pbem = 0.2,
                         coverage_min_to_compute_pbem = 10,
                         n_pos_af_th = 0.2,
                         mc.cores = 1,
                         step = 5000){
  cat(paste("[",Sys.time(),"]\tReading the sample.info.file","\n"))
  sif <- import_sif(main_sif = sample.info.file)

  cat(paste("[",Sys.time(),"]\tReading chromosomes from 'targetbed'","\n"))
  chromosomes <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)[[1]]

  cat(paste("[",Sys.time(),"]\tComputation of per-base error model","\n"))

  if(!file.exists(file.path(outdir, outdir.bperr.name))){
    dir.create(file.path(outdir, outdir.bperr.name), showWarnings = TRUE)
  }
  setwd(file.path(outdir, outdir.bperr.name))

  for(chrom in unique(chromosomes)){
    cat(paste("[",Sys.time(),"]\tchromosome:",chrom),"\n")

    tp <- bed2positions(targetbed = targetbed,chrom_to_extract = chrom,get_only_chromosomes = FALSE)
    targets <- unique(tp$PosByChrom)

    cat(paste("[",Sys.time(),"]\ttotal positions to check in this chromosome :",nrow(targets)),"\n")
    mclapply(seq(1,nrow(targets),step),pos2bperr,
             targets=targets,
             germlineset=get_germlineset(sifgerm = sif$df_ctrl,pacbamfolder_bychrom = pacbamfolder_bychrom,chrom = chrom),
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
  cmd.merge = paste("cat afgtz_*.tsv > afgtz.tsv")
  system(cmd.merge)
  cmd.merge = paste("cat afz_*.tsv > afz.tsv")
  system(cmd.merge)

  # save counter of afz as RData
  afztab = read.delim(file = "afz.tsv",sep="\t",as.is = T,header=F)
  afz = apply(afztab,2,sum)
  names(afz)=define_cov_bins(coverage_binning)[[2]]
  save(afz,file = "afz.RData",compress = T)

  # overall statistics on pbems
  cmd.merge = paste("cat bperr_*.tsv > bperr.tsv")
  system(cmd.merge)
  pbem_tab <- fread(file.path(outdir, outdir.bperr.name,"bperr.tsv"),stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5),data.table = F)
  header_pbem_tab <- c("group","chr","pos","ref","dbsnp","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth","bperr","tot_reads_supporting_alt")
  colnames( pbem_tab ) <- header_pbem_tab
  rownames( pbem_tab ) <- pbem_tab$group

  # summary stats for pbem across the target
  bperr_summary = summary(pbem_tab$bperr) # bperr
  save(bperr_summary,file = file.path(outdir, outdir.bperr.name,"bperr_summary.RData"))

  # compute background pbem
  bperr_subset <- pbem_tab[which(pbem_tab$n_pos_af_gth == 0),] # select only pos never > n_pos_af_gth
  bgpbem = (sum(as.numeric(bperr_subset$tot_reads_supporting_alt)))/(sum(as.numeric(bperr_subset$tot_coverage)))
  mean_pbem = mean(as.numeric(bperr_subset$bperr),na.rm = T)
  save(bgpbem,mean_pbem,file = file.path(outdir, outdir.bperr.name,"pbem_background.RData"))

  save(pbem_tab,file = file.path(outdir, outdir.bperr.name,"pbem_tab.RData"),compress = T)

  cat(paste("[",Sys.time(),"]\talright.","\n"))
  return(list(pbem_tab=pbem_tab,
              bperr_summary=bperr_summary,
              bgpbem=bgpbem,
              mean_pbem=mean_pbem))
}

