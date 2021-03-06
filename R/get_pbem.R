get_pbem = function(id,tgf){
  this = tgf[id,,drop=FALSE]
  alts = setdiff(c("total.A","total.C","total.G","total.T"),paste0("total.",this$ref))
  rd.alts = as.numeric(sum(this[alts],na.rm = TRUE))
  this$bperr = sum(rd.alts)/this$tot_coverage
  this$tot_reads_supporting_alt = rd.alts
  return(this)
}
