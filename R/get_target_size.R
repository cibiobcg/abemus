#' Compute total number of loci covered in a BED file
#' @param targetbed targeted regions in BED format
#' @param Mbp return count as Mbp. default: TRUE
#' @examples
#' targetbed <- system.file("extdata", "regions_toy.bed", package = "abemus")
#' target_size <- get_target_size(targetbed=targetbed, Mbp = TRUE)
#' @return Mbp covered in the BED file
#' @export
get_target_size <- function(targetbed,Mbp = TRUE){
  bed <- fread(input = targetbed,colClasses = list(character=1),data.table = FALSE,stringsAsFactors = FALSE,header = FALSE)
  if(Mbp){
    return( sum(bed$V3-bed$V2)/1e+6 )
  } else {
    return( sum(bed$V3-bed$V2) )
  }
}
