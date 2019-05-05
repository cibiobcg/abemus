#' bed2positions
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @param get_only_chromosomes return only the list of chromosomes in the BED file
#' @param chrom_to_extract the chromosome name to extract as i.e. "chr1"
#' @return PosByChrom A list of data.frames, each for each chromosome.
#' @return chromosomes The chromosomes included in the BED file.
#' @export
bed2positions <- function(targetbed,
                          chrom_to_extract,
                          get_only_chromosomes = F){
  if(get_only_chromosomes){
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = F,stringsAsFactors = F,header = F)
    bed$V1 <- gsub(bed$V1,pattern = "chr",replacement = "")
    bed$V1 <- paste0("chr",bed$V1)
    chromosomes <- sort(unique(bed$V1))

    return(list(chromosomes=chromosomes))

  } else {
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = F,stringsAsFactors = F,header = F)
    bed$V1 <- gsub(bed$V1,pattern = "chr",replacement = "")
    bed$V1 <- paste0("chr",bed$V1)
    bed$V2 <- bed$V2+1
    chromosomes <- sort(unique(bed$V1))

    unwrap <- function(x){return(seq.int(from = x[2],to = x[3]))}

    PosByChrom <- data.frame(chr=chr,
                             pos=as.numeric(unlist( apply(bed_chrom, MARGIN=1,FUN = unwrap))),
                             stringsAsFactors = F)
    PosByChrom$group <- paste(PosByChrom$chr,PosByChrom$pos,sep=":")

    return(list(chromosomes=chromosomes, PosByChrom=PosByChrom))

  }
}

