# Read and check the sample info file
#
# @param main_sif A tab delimeted file without header where columns are: c("patient","plasma","path/to/plasma.bam","germline","path/to/germline.bam")
# @return A list containing 2 data frames:
# data.frame with unique germline samples (use for GSE distribution and pbem computation steps)
# data.frame only case samples having or not a matched control sample (use for call snvs step)
import_sif <- function(main_sif){
  df  <-  read.delim(main_sif,as.is=T,stringsAsFactors = F,header = F)
  colnames(df) <- c("patient","plasma","plasma.bam","germline","germline.bam")
  # remove NAs and keep only unique germline samples
  df_ctrl <- df[which(!is.na(df$germline.bam)),,drop=FALSE]
  df_ctrl <- unique(df_ctrl[,c(4,5)])
  # remove NAs and keep only case samples having matched germline samples
  df_cases <- df[which(!is.na(df$plasma.bam)),,drop=FALSE]
  return(list(df_ctrl=df_ctrl,df_cases=df_cases))
}
