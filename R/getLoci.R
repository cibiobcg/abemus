getLoci <- function(chr,pacbam,select='ref',nThread=4){
  p <- list.files(pacbam,pattern = paste0('_chr',chr,'.pileup$'),full.names = TRUE)[1]
  if (is.na(p)) {
    return(NA)
  } else{
    df <- fread(input = p,
                sep = '\t',
                stringsAsFactors = FALSE,
                header = TRUE,
                select = select,
                verbose = FALSE,
                data.table = FALSE,
                na.strings = '',
                nThread = nThread)
    return(as.vector(unlist(df)))
  }
}
