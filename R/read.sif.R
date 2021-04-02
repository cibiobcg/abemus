read.sif <- function(file){
  sif <- list(ctrl=NA)
  tsv <- read.delim(file,header = FALSE,sep = '\t',stringsAsFactors = FALSE)
  if(ncol(tsv) == 1){
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- sif$ctrl
  } else{
    sif$case <- gsub(unique(basename(tsv[,1])),pattern = '\\.bam$',replacement = '')
    sif$ctrl <- gsub(basename(tsv[,2]),pattern = '\\.bam$',replacement = '')
  }
  return(sif)
}
