#' Get list of pileups for germline samples
#'
#' @param main_sif sample info file
#' @param pacbamfolder folder with pileups
#' @param chrom chromosome
#' @return A character string with the path to germline pileups
get_germlineset <- function(main_sif,pacbamfolder,chrom){
  germlineset = c()
  for(id in 1:nrow(main_sif)){
    thisSample=main_sif[id,]
    name = gsub(basename(thisSample$germline.bam),pattern = ".bam",replacement = "")
    thisPileup.file = list.files(file.path(pacbamfolder,name,"pileup"),full.names = T,pattern = paste0(chrom,"\\.pileup$"))
    germlineset = c(germlineset,thisPileup.file)
  }
  germlineset = paste(unique(germlineset),collapse = " ")
  return(germlineset)
}
