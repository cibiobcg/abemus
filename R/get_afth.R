get_afth <- function(num,
                     idx,
                     next.bin.AFgtz,
                     next.bin.AFetz,
                     vafcov,
                     this.detection.specificity){
  if(length(idx) < next.bin.AFgtz){
    b = vafcov[sample(x = idx,size = next.bin.AFgtz,replace = TRUE),1]
  } else {
    b = vafcov[sample(x = idx,size = next.bin.AFgtz,replace = FALSE),1]
  }
  afth = quantile.zaf(x = b,probs = this.detection.specificity,nz = next.bin.AFetz)
  return(afth)
}
