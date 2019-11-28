#' Compute total number of loci covered in a BED file
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @param Mbp Return count as Mbp. default: TRUE
#' @return Mbp covered in the \code{targetbed}
#' @export
#' @examples
#' targetbed <- system.file("extdata", "regions_toy.bed", package = "abemus")
#' target_size <- get_target_size(targetbed=targetbed, Mbp = TRUE)
get_target_size <- function(targetbed,Mbp = TRUE){
  bed <- fread(input = targetbed,colClasses = list(character=1),data.table = FALSE,stringsAsFactors = FALSE,header = FALSE)
  if(Mbp){
    return( sum(bed$V3-bed$V2)/1e+6 )
  } else {
    return( sum(bed$V3-bed$V2) )
  }
}
