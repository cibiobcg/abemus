chrom_cov_vaf <- function(i,mat_cov,mat_vaf,rsid,max.vaf,min.cov){

  if(all(is.na(mat_cov[[i]]))){
    return(list(NA,NA))
  }

  cov <- as.vector(mat_cov[[i]])
  cov[which(cov < min.cov)] <- NA

  vaf <- as.vector(mat_vaf[[i]])
  vaf[which(vaf > max.vaf)] <- NA
  vaf[rsid[[i]]] <- NA

  return(list(cov,vaf))
}
