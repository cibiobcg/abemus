#' quantile.zaf
#'
#' @param x vafcov[idx,1]
#' @param probs quantiles
#' @param nz zeros
#' @return quantiles
quantile.zaf <- function (x, probs = seq(0, 1, 0.25),nz){
  N <- length(x) + round(nz)
  x <- sort(x)
  qnt <- 1 + round((probs * (N - 1))) # !
  ql <- c()
  for(j in qnt){
    if(j-nz <=0){
      ql = c(ql,0)
    } else {
      ql = c(ql,x[j-nz])
    }
  }
  names(ql) <- probs
  return(ql)
}
