#' Compute the per-base error measure (pbem) for each targeted locus and save allelic (AF) by bins of coverage.
#' @param sample.info.file The sample info file listing CASE and CONTROL samples. The format is simply 5 columns, tab-delimited, and there is no column header.
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @param pacbamfolder_bychrom The folder popluted by outputs by the \code{split_pacbam_bychrom} function.
#' @param outdir The folder where outputs will be saved.
#' @param outdir.bperr.name The subfolder name that will be created in the \code{outdir}. default: "BaseErrorModel"
#' @param coverage_binning Bins of coverage into which divide AFs. default: 50
#' @param af_max_to_compute_thresholds To compute AF thresholds, consider only positions with AF <= \code{af_max_to_compute_thresholds}. default 0.2
#' @param coverage_min_to_compute_thresholds To compute AF threshold, consider only positions with coverage >= \code{coverage_min_to_compute_thresholds}. default 10
#' @param af_max_to_compute_pbem To compute pbem, consider only positions with AF <= \code{af_max_to_compute_pbem}. default: 0.2
#' @param coverage_min_to_compute_pbem To compute pbem, consider only positions with coverage >= \code{coverage_min_to_compute_pbem}. default: 10
#' @param n_pos_af_th When compute pbem, count in how many germline samples the position has an AF >= \code{n_pos_af_th}. default: 0.2
#' @param mc.cores Number of jobs to run in parallel, it is the \code{mc.core} param of the \code{mclapply} function. default: 1
#' @param step Number of positions into split the input chromosome file. default: 5000
#' @return The \code{compute_pbem} will write, in the \code{outdir.bperr.name}, tab-delimeted files reporting the per-base error measure of each targeted locus (saved in \code{bperr.tsv}), the coverage and the AF of loci with AF > 0 (saved in \code{afgtz.tsv}) and the coverage of loci with AF = 0 (saved in \code{afz.tsv}).
#' @return The function also return objects \code{bgpbem}, \code{bperr_summary} and \code{mean_pbem} reporting overall statistics about the per-base error measure.
#' @export
#' @examples
#' sample.info.file <- system.file("extdata", "test_sif_toy.tsv", package = "abemus")
#' targetbed <- system.file("extdata", "regions_toy.bed", package = "abemus")
#' pacbamfolder_bychrom <- system.file("extdata", "pacbam_data_bychrom", package = "abemus")
#' outdir <- tempdir()
#' outpbem <- compute_pbem(sample.info.file,targetbed,pacbamfolder_bychrom,outdir)
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

  old_wd <- getwd()
  on.exit(setwd(old_wd))

  message(paste("[",Sys.time(),"]\tReading the sample.info.file"))
  sif <- import_sif(main_sif = sample.info.file)

  message(paste("[",Sys.time(),"]\tReading chromosomes from 'targetbed'"))
  chromosomes <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)[[1]]

  message(paste("[",Sys.time(),"]\tComputation of per-base error model"))

  if(!file.exists(file.path(outdir, outdir.bperr.name))){
    dir.create(file.path(outdir, outdir.bperr.name), showWarnings = TRUE)
  } else{
    file.remove(list.files(file.path(outdir, outdir.bperr.name), full.names = TRUE))
  }
  setwd(file.path(outdir, outdir.bperr.name))

  for(chrom in unique(chromosomes)){
    message(paste("[",Sys.time(),"]\tchromosome:",chrom))

    tp <- bed2positions(targetbed = targetbed,chrom_to_extract = chrom,get_only_chromosomes = FALSE)
    targets <- unique(tp$PosByChrom)

    message(paste("[",Sys.time(),"]\ttotal positions to check in this chromosome :",nrow(targets)))
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

    # pbem.table
    merge.files <- list.files(path = ".",pattern = "_pbem.table.txt",full.names = TRUE)
    mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
    sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file=paste0("bperr_",chrom,".tsv"),append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

    remove.files <- list.files(path = ".",pattern = "_pbem.table.txt",full.names = TRUE)
    do.call(file.remove,list(remove.files))

    # afgtz.table
    merge.files <- list.files(path = ".",pattern = "_afgtz.table.txt",full.names = TRUE)
    mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
    sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file=paste0("afgtz_",chrom,".tsv"),append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

    remove.files <- list.files(path = ".",pattern = "_afgtz.table.txt",full.names = TRUE)
    do.call(file.remove,list(remove.files))

    # afz.table
    merge.files <- list.files(path = ".",pattern = "_afz.table.txt",full.names = TRUE)
    mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
    sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file=paste0("afz_",chrom,".tsv"),append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

    remove.files <- list.files(path = ".",pattern = "_afz.table.txt",full.names = TRUE)
    do.call(file.remove,list(remove.files))

  }

  # create afgtz.tsv
  merge.files <- list.files(path = ".",pattern = "afgtz_",full.names = TRUE)
  mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
  sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file="afgtz.tsv",append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

  # create afz.tsv
  merge.files <- list.files(path = ".",pattern = "afz_",full.names = TRUE)
  mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
  sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file="afz.tsv",append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

  # save counter of afz as RData
  afztab = read.delim(file = "afz.tsv",sep="\t",as.is = TRUE,header=FALSE)
  afz = apply(afztab,2,sum)
  names(afz)=define_cov_bins(coverage_binning)[[2]]
  save(afz,file = "afz.RData",compress = T)

  # overall statistics on pbems
  merge.files <- list.files(path = ".",pattern = "bperr_",full.names = TRUE)
  mrgd <- lapply(merge.files,fread,data.table=FALSE,stringsAsFactors = FALSE,header = FALSE)
  sapply(seq_len(length(mrgd)),function (x) write.table(mrgd[[x]],file="bperr.tsv",append = TRUE,quote=FALSE,col.names=FALSE,row.names = FALSE,sep = "\t"))

  pbem_tab <- fread(file.path(outdir, outdir.bperr.name,"bperr.tsv"),stringsAsFactors = FALSE,showProgress = FALSE,header = FALSE,colClasses = list(character=2,character=5),data.table = FALSE)
  header_pbem_tab <- c("group","chr","pos","ref","dbsnp","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth","bperr","tot_reads_supporting_alt")
  colnames( pbem_tab ) <- header_pbem_tab
  rownames( pbem_tab ) <- pbem_tab$group

  # summary stats for pbem across the target
  bperr_summary = summary(pbem_tab$bperr) # bperr
  save(bperr_summary,file = file.path(outdir, outdir.bperr.name,"bperr_summary.RData"))

  # compute background pbem
  bperr_subset <- pbem_tab[which(pbem_tab$n_pos_af_gth == 0),] # select only pos never > n_pos_af_gth
  bgpbem = (sum(as.numeric(bperr_subset$tot_reads_supporting_alt)))/(sum(as.numeric(bperr_subset$tot_coverage)))
  mean_pbem = mean(as.numeric(bperr_subset$bperr),na.rm = TRUE)
  save(bgpbem,mean_pbem,file = file.path(outdir, outdir.bperr.name,"pbem_background.RData"))

  save(pbem_tab,file = file.path(outdir, outdir.bperr.name,"pbem_tab.RData"),compress = TRUE)

  message(paste("[",Sys.time(),"]\talright."))
  return(list(pbem_tab=pbem_tab,
              bperr_summary=bperr_summary,
              bgpbem=bgpbem,
              mean_pbem=mean_pbem))
}

