CheckAltReads <- function(i,snvs){
  this = snvs[i,,drop=F]
  altbase = this$alt
  refbase = this$ref
  cov.ref = this[,refbase]
  rev.ref = this[,paste0(refbase,"rs")]
  fwd.ref = cov.ref - rev.ref
  this$rev.ref = rev.ref
  this$fwd.ref = fwd.ref
  if(altbase=="N"){
    altbases = setdiff(c("A","C","G","T"),refbase)
    cov.alt = as.numeric(max(this[,altbases]))
    this$cov.alt <- cov.alt
    this$af <- round(cov.alt/sum(cov.alt,this[,refbase]),4)
    this$rev.alt = NA
    this$fwd.alt = NA
    this$strandbias = NA
  } else {
    cov.alt = this[,altbase]
    this$cov.alt <- cov.alt
    rev.alt = this[,paste0(altbase,"rs")]
    fwd.alt = cov.alt - rev.alt
    this$rev.alt = rev.alt
    this$fwd.alt = fwd.alt
    this$strandbias = strandbias(fwd.ref=fwd.ref,fwd.alt=fwd.alt,rev.ref=rev.ref,rev.alt=rev.alt)
  }
  return(this)
}
