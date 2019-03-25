#' apply_AF_filters
#'
apply_AF_filters <- function(chrpmF1,
                             AFbycov,
                             af.threshold.table,
                             minaf,
                             mybreaks,
                             mc.cores){
  if (AFbycov == FALSE & !is.numeric(AFbycov)){
    chrpmF1[,'af_threshold'] <- minaf
  } else if (is.numeric(AFbycov)){
    chrpmF1[,'af_threshold'] <- AFbycov
  } else if (AFbycov == TRUE & !is.numeric(AFbycov)){
    thresholds = as.numeric(af.threshold.table[,-1])
    af_filter_by_coverage <- function(y,thresholds,chrpmF1){
      this = chrpmF1[y,,drop=F]
      minaf_covth = thresholds[findInterval(this$cov_case,mybreaks)]
      this[,'af_threshold'] <- minaf_covth
      return(this)
    }
    out = mclapply(seq(1,nrow(chrpmF1),1),af_filter_by_coverage,thresholds=thresholds,chrpmF1=chrpmF1,mc.cores = mc.cores)
    chrpmF1 = fromListToDF(out)
  }
  return(chrpmF1)
}
