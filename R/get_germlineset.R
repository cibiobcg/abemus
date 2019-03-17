#' Get list of pileups for germline samples
#'
#' @param sif
#' @param pacbamfolder
#' @param chrom
#' @return A character string with the path to germline pileups
get_germlineset <- function(sif,pacbamfolder,chrom){
  germlineset = c()
  for(id in 1:nrow(sif)){
    thisSample=sif[id,]
    name = gsub(basename(thisSample$germline.bam),pattern = ".bam",replacement = "")
    thisPileup.file = list.files(file.path(pacbamfolder,name,"pileup"),full.names = T,pattern = paste0(chrom,"\\.pileup$"))
    germlineset = c(germlineset,thisPileup.file)
  }
  germlineset = paste(unique(germlineset),collapse = " ")
  return(germlineset)
}
