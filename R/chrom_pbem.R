chrom_pbem <- function(i,mat_cov,mat_cov_base,ref,max.vaf.pbem = 0.2,cov.min.pbem = 10){

  if(all(is.na(ref[[i]]))){
    return(NA)
  }

  x <- rep(NA,length(ref[[i]]))

  mat_cov_chr <- mat_cov[[i]]

  filtout <- unique(rbind(which(mat_cov_chr <= cov.min.pbem,arr.ind = TRUE),
                          which(mat_vaf[[i]] >= max.vaf.pbem,arr.ind = TRUE)))

  mat_cov_chr[filtout] <- NA

  den <- rowSums(mat_cov_chr,na.rm = TRUE)

  ext <- function(lst,n){
    sapply(lst,FUN = '[',n)
  }

  get_num <- function(alt,idx,filtout){
    alt[filtout,] <- NA
    rowSums(alt[idx,],na.rm = TRUE)
  }

  for(base in c('A','C','G','T')){

    alt_list <- ext(mat_cov_base[setdiff(c('A','C','G','T'),base)],n = i)

    idx <- which(ref[[i]] == base)

    num <- sapply(alt_list, get_num, idx, filtout)

    x[idx] <- rowSums(num)/den[idx]

  }

  return(x)
}
