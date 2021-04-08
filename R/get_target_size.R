get_target_size <- function(targetbed,Mbp = TRUE){
  bed <- fread(input = targetbed,colClasses = list(character=1),data.table = FALSE,stringsAsFactors = FALSE,header = FALSE)
  if(Mbp){
    return( sum(bed$V3-bed$V2)/1e+6 )
  } else {
    return( sum(bed$V3-bed$V2) )
  }
}
