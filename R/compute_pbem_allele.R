compute_pbem_allele <- function(id,abemus){
  locus = abemus[id,,drop=F]
  locus$pbem_allele = locus[,paste0('total.',locus$alt)]/locus$tot_coverage
  # check if same alt allele is observed in case and germline (when af germline > 0)
  alt_case = as.numeric(locus[,paste0(locus$alt,'_case')])
  alt_ctrl = as.numeric(locus[,paste0(locus$alt,'_control')])
  if(alt_ctrl>0){
    locus$same_allele = 1
  } else {
    locus$same_allele = 0
  }
  return(locus)
}
