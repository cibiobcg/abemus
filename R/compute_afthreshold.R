#' Compute Allelic Fraction thresholds
#' @export
#' @param outdir output folder for this step analysis
#' @param pbem_dir folder with pbem data. default: file.path(outdir, "BaseErrorModel")
#' @param outdir.afth.name folder will be created in outdir. default: "Controls"
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param probs quatiles to compute
#' @examples
#' outdir <- tempdir()
#' pbem_dir <- system.file("extdata", "BaseErrorModel", package = "abemus")
#' outafth <- compute_afthreshold(outdir = outdir,pbem_dir = pbem_dir)
#' @return allelic fraction thresholds
compute_afthreshold <- function(outdir,
                                pbem_dir = file.path(outdir,"BaseErrorModel"),
                                outdir.afth.name = "Controls",
                                coverage_binning = 50,
                                probs = seq(0.9,1,0.0001)){

  old_wd <- getwd()
  on.exit(setwd(old_wd))

  message(paste("[",Sys.time(),"]\tEstimation of AF thresolds exploiting germline samples only"))

  if(!file.exists(file.path(outdir, outdir.afth.name))){
    dir.create(file.path(outdir, outdir.afth.name), showWarnings = TRUE)
  }
  setwd(file.path(outdir, outdir.afth.name))

  vafcov_file = file.path(pbem_dir,"afgtz.tsv")
  afz_file = file.path(pbem_dir,"afz.RData")

  #  check if data tables exists
  if(file.exists(vafcov_file)){
    message(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ ok ]"))
    vafcov = fread(input = vafcov_file,sep = "\t",header = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
  } else {
    message(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ not found ]"))
    stop()
  }

  if(file.exists(afz_file)){
    message(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ ok ]"))
    load(afz_file)
  } else {
    message(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ not found ]"))
    stop()
  }

  message(paste("[",Sys.time() ,"]\taf.all",dim(vafcov)[1] + sum(afz),"af.gtz",dim(vafcov)[1],"af.etz",sum(afz)))
  message(paste("[",Sys.time(),"]\tcompute AF quantiles +++ not stratified by coverage +++"))
  nz <-  sum(afz)
  th_results <-  quantile.zaf(x = as.numeric(vafcov[,1]),probs = probs,nz = nz)

  message(paste("[",Sys.time(),"]\tcompute AF quantiles +++ stratified by coverage +++"))

  covbin <- define_cov_bins(coverage_binning)[[1]]
  lev <- define_cov_bins(coverage_binning)[[2]]

  a <- cut(x = as.integer(vafcov[,2]),breaks = covbin,include.lowest = TRUE)
  th_results_bin = data.frame(specificity = probs)
  datacount_bin = array(data = 0,dim = length(lev),dimnames = list(lev))
  for(l in lev){
    message(paste("[",Sys.time(),"]\tprocessing coverage bin\t",l))
    idx = which(a==l)
    datacount_bin[l] <- length(idx)
    if(length(idx)>0){
      b = vafcov[idx,1]
      nz = as.integer(afz[l])
      kk = quantile.zaf(x = b,probs = probs,nz = nz)
      th_results_bin[,l] <- kk
    } else {
      th_results_bin[,l] <- NA
    }
  }

  message(paste("[",Sys.time(),"]\tSaving threshold data"))
  save(th_results,th_results_bin,datacount_bin,file = file.path(outdir, outdir.afth.name,paste0("datathreshold.RData")),compress = TRUE)

  message(paste("[",Sys.time(),"]\talright."))
  return(list(th_results=th_results,
              th_results_bin=th_results_bin,
              datacount_bin=datacount_bin))
}
