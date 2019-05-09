define_cov_bins <- function(coverage_binning,min.coverage=0,max.coverage=5000){
  covbin = seq(min.coverage,max.coverage,by = coverage_binning)
  covbin[length(covbin)] <- Inf
  lev = levels(cut(1,breaks=covbin,include.lowest=TRUE))
  return(list(covbin,lev))
}
