#' Compute afthresholds
#'
#' @export
#' @param outdir output folder for this step analysis
#' @param pbem_dir folder with pbem data. default: "outdir/BaseErrorModel"
#' @param outdir.afth.name folder will be created in outdir. default: "Controls"
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param probs quatiles to compute
#' @return allelic fraction thresholds
compute_afthreshold <- function(outdir,
                                pbem_dir = file.path(outdir,"BaseErrorModel"),
                                outdir.afth.name = "Controls",
                                coverage_binning = 50,
                                probs = seq(0.9,1,0.0001)){
  cat(paste("[",Sys.time(),"]\tEstimation of AF thresolds exploiting germline samples only","\n"))

  if(!file.exists(file.path(outdir, outdir.afth.name))){
    dir.create(file.path(outdir, outdir.afth.name), showWarnings = TRUE)
  }
  setwd(file.path(outdir, outdir.afth.name))

  vafcov_file = file.path(pbem_dir,"afgtz.tsv")
  afz_file = file.path(pbem_dir,"afz.RData")

  #  check if data tables exists
  if(file.exists(vafcov_file)){
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ ok ]"),"\n")
    #vafcov = read.big.matrix(filename = vafcov_file,header = F,sep = "\t",type = "double")
    vafcov = fread(input = vafcov_file,sep = "\t",header = F,stringsAsFactors = F,data.table = F)
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ not found ]"),"\n")
    quit()
  }

  if(file.exists(afz_file)){
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ ok ]"),"\n")
    load(afz_file)
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ not found ]"),"\n")
    quit()
  }

  cat(paste("[",Sys.time() ,"]\taf.all",dim(vafcov)[1] + sum(afz),"af.gtz",dim(vafcov)[1],"af.etz",sum(afz)),"\n")
  cat(paste("[",Sys.time(),"]\tcompute AF quantiles +++ not stratified by coverage +++"),"\n")
  nz <-  sum(afz)
  th_results <-  quantile.zaf(x = as.numeric(vafcov[,1]),probs = probs,nz = nz)
  #minaf <- as.numeric(th_results[as.character(spec)])

  cat(paste("[",Sys.time(),"]\tcompute AF quantiles +++ stratified by coverage +++"),"\n")

  covbin <- define_cov_bins(coverage_binning)[[1]]
  lev <- define_cov_bins(coverage_binning)[[2]]

  a <- cut(x = as.integer(vafcov[,2]),breaks = covbin,include.lowest = T)
  th_results_bin = data.frame(specificity = probs)
  datacount_bin = array(data = 0,dim = length(lev),dimnames = list(lev))
  for(l in lev){
    cat(paste("[",Sys.time(),"]\tprocessing coverage bin\t",l),"\n")
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
  #minaf_cov = th_results_bin[which(round(th_results_bin$specificity,4)==spec),]

  cat(paste("[",Sys.time(),"]\tSaving threshold data\n"))
  #save(minaf,minaf_cov,th_results,th_results_bin,covbin,lev,datacount_bin,file = file.path(outdir, outdir.afth.name,paste0("datathreshold.RData")),compress = T)
  save(th_results,th_results_bin,datacount_bin,file = file.path(outdir, outdir.afth.name,paste0("datathreshold.RData")),compress = T)
  return(list(th_results,th_results_bin,datacount_bin))
}
