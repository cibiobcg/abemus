#' Unwrap a list of genomic intervals into a list of loci grouped by chromosome.
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @param get_only_chromosomes Set TRUE to return only the list of \code{chromosomes} present in the \code{targetbed}. default: FALSE
#' @param chrom_to_extract The chromosome name to extract as i.e. "chr1" or "1" accordingly with BED file annotation
#' @return \code{PosByChrom}:  A list of data.frames (each row is a locus) one for each chromosome present in the \code{targetbed}.
#' @return \code{chromosomes}: A vector listing chromosomes present in the \code{targetbed}.
#' @examples
#' targetbed <- system.file("extdata", "regions_toy.bed",package = "abemus")
#' chromosomes <- bed2positions(targetbed = targetbed,chrom_to_extract="8")
#' targetbp_list <- bed2positions(targetbed = targetbed,get_only_chromosomes = TRUE)
#' @export
bed2positions <- function(targetbed,
                          chrom_to_extract,
                          get_only_chromosomes = FALSE){
  if(get_only_chromosomes){
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = FALSE,stringsAsFactors = FALSE,header = FALSE)
    chromosomes <- sort(unique(gsub(bed$V1,pattern = 'chr',replacement = '')))

    return(list(chromosomes=chromosomes))

  } else {
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = FALSE,stringsAsFactors = FALSE,header = FALSE)
    bed$V2 <- bed$V2+1
    chromosomes <- sort(unique(gsub(bed$V1,pattern = 'chr',replacement = '')))

    unwrap <- function(x){return(seq.int(from = x[2],to = x[3]))}

    bed_chrom <- bed[grep(gsub(bed$V1,pattern = 'chr',replacement = ''),pattern = as.character(chrom_to_extract)),,drop=FALSE]

    PosByChrom <- data.frame(chr=chrom_to_extract,
                             pos=as.numeric(unlist( apply(bed_chrom, MARGIN=1,FUN = unwrap))),
                             stringsAsFactors = FALSE)

    return(list(chromosomes=chromosomes, PosByChrom=PosByChrom))

  }
}

