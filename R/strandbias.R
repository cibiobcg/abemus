#' strandbias
#'
strandbias = function(fwd.ref,fwd.alt,rev.ref,rev.alt){
  if(fwd.alt+rev.alt==0 | fwd.ref+rev.ref==0){
    return(NA)
  } else {
    #sb = abs((fwd.alt/(fwd.ref+fwd.alt))-(rev.alt/(rev.ref+rev.alt)))/((fwd.alt+rev.alt)/sum(fwd.ref,fwd.alt,rev.ref,rev.alt))
    sb = abs((fwd.alt/(fwd.ref+fwd.alt))-(rev.alt/(rev.ref+rev.alt)))
    sb = round(sb,5)
    return(sb)
  }
}
