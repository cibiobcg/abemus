#' add_class
#'
add_class = function(pmtab){
  pmtab$CLASS = NA
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr==0 & pmtab$pbem_allele==0)] = 1
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr>0 & pmtab$pbem_allele==0)] = 2
  pmtab$CLASS[which(pmtab$af_control==0 & pmtab$bperr>0 & pmtab$pbem_allele>0)] = 3
  pmtab$CLASS[which(pmtab$af_control>0 & pmtab$bperr>0 & pmtab$pbem_allele>=0 & pmtab$same_allele == 0)] = 4
  pmtab$CLASS[which(pmtab$af_control>0 & pmtab$bperr>0 & pmtab$pbem_allele>0 & pmtab$same_allele == 1)] = 5
  return(pmtab)
}
