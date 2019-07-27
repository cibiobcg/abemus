#' bed2positions
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @param get_only_chromosomes return only the list of chromosomes in the BED file. default: FALSE
#' @param chrom_to_extract the chromosome name to extract as i.e. "chr1" or "1" accordingly with BED file annotation
#' @return PosByChrom A list of data.frames, each for each chromosome.
#' @return chromosomes The chromosomes included in the BED file.
#' @export
bed2positions <- function(targetbed,
                          chrom_to_extract,
                          get_only_chromosomes = F){
  if(get_only_chromosomes){
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = F,stringsAsFactors = F,header = F)
    chromosomes <- sort(unique(bed$V1))

    return(list(chromosomes=chromosomes))

  } else {
    bed <- fread(input = targetbed,colClasses = list(character=1),data.table = F,stringsAsFactors = F,header = F)
    bed$V2 <- bed$V2+1
    chromosomes <- sort(unique(bed$V1))

    unwrap <- function(x){return(seq.int(from = x[2],to = x[3]))}

    bed_chrom <- bed[grep(bed$V1,pattern = chrom_to_extract),,drop=F]

    PosByChrom <- data.frame(chr=chrom_to_extract,
                             pos=as.numeric(unlist( apply(bed_chrom, MARGIN=1,FUN = unwrap))),
                             stringsAsFactors = F)
    #PosByChrom$group <- paste(PosByChrom$chr,PosByChrom$pos,sep=":")

    return(list(chromosomes=chromosomes, PosByChrom=PosByChrom))

  }
}

