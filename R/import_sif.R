#' Read and check the sample info file
#'
#' @param main_sif A tab delimeted file without header where columns are: c("patient","plasma","plasma.bam","germline","germline.bam")
#' @return A list containing 2 data frames:
#' [[1]] Table with unique germline samples (use for GSE distribution and pbem computation steps)
#' [[2]] Table with only case samples having a matched control sample (use for call snvs step)
import_sif <- function(main_sif){
  df  <-  read.delim(main_sif,as.is=T,stringsAsFactors = F,header = F)
  colnames(df) <- c("patient","plasma","plasma.bam","germline","germline.bam")
  # remove NAs and keep only unique germline samples
  df_ctrl = df[which(!is.na(df$germline.bam)),]
  df_ctrl = unique(df_ctrl[,c(1,4,5)])
  # remove NAs and keep only case samples having matched germline samples
  nas = which(is.na(df$plasma.bam) | is.na(df$germline.bam))
  if(length(nas)>0){df_cases = df[-nas,]}
  return(list(df_ctrl,df_cases))
}
