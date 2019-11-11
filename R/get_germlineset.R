get_germlineset <- function(sifgerm,pacbamfolder_bychrom,chrom){
  germlineset = c()
  chrom <- gsub(chrom,pattern = 'chr',replacement = '')
  for(id in 1:nrow(sifgerm)){
    thisSample=sifgerm[id,]
    name = gsub(basename(thisSample$germline.bam),pattern = ".bam",replacement = "")
    thisPileup.file = list.files(file.path(pacbamfolder_bychrom,name,"pileup"),full.names = TRUE,pattern = paste0("_chr",chrom,"\\.pileup$"))
    germlineset = c(germlineset,thisPileup.file)
  }
  #germlineset = paste(unique(germlineset),collapse = " ")
  germlineset <- unique(germlineset)
  return(germlineset)
}
