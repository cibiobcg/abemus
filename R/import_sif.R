import_sif <- function(main_sif){
  df  <-  read.delim(main_sif,as.is=TRUE,stringsAsFactors = FALSE,header = FALSE)
  colnames(df) <- c("patient","plasma","plasma.bam","germline","germline.bam")
  # remove NAs and keep only unique germline samples
  df_ctrl <- df[which(!is.na(df$germline.bam)),,drop=FALSE]
  df_ctrl <- unique(df_ctrl[,c(4,5)])
  # remove NAs and keep only case samples having matched germline samples
  df_cases <- df[which(!is.na(df$plasma.bam)),,drop=FALSE]
  return(list(df_ctrl=df_ctrl,df_cases=df_cases))
}
