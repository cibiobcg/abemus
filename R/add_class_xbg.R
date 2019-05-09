add_class_xbg = function(pmtab,xbg){
  if(nrow(pmtab)>0){
    pmtab$CLASS.xbg = NA
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr<=xbg & pmtab$pbem_allele<=xbg)] = 1
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele<=xbg)] = 2
    pmtab$CLASS.xbg[which(pmtab$af_control<=xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg)] = 3
    pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>=xbg & pmtab$same_allele == 0)] = 4
    pmtab$CLASS.xbg[which(pmtab$af_control>xbg & pmtab$bperr>xbg & pmtab$pbem_allele>xbg & pmtab$same_allele == 1)] = 5
    return(pmtab)
  } else {
    return(pmtab)
  }
}
