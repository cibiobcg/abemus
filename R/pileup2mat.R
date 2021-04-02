pileup2mat <- function(chr,sif,pacbam,select,sparse=FALSE,nThread=4){
  id <- paste0(unique(na.omit(sif$ctrl)),paste0('_chr',chr),'.pileup')
  ff <- file.path(pacbam,id)
  ff <- ff[which(file.exists(ff))]
  if(length(ff) == 0){
    return(NA)
  } else{
    df <- lapply(ff,
                 fread,
                 sep = '\t',
                 stringsAsFactors = FALSE,
                 header = TRUE,
                 select = select,
                 verbose = FALSE,
                 data.table = FALSE,
                 na.strings = '',
                 nThread = nThread)
    if(sparse){
      return(Matrix(as.matrix(do.call(cbind,df)), sparse = TRUE))
    } else {
      return(as.matrix(do.call(cbind,df)))
    }
  }
}
