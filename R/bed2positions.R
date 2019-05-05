#' bed2positions
#' @param targetbed Genomic regions in the BED tab-delimited format.
#' @return PosByChrom A list of data.frames, each for each chromosome.
#' @return chromosomes The chromosomes included in the BED file.
#' @export
bed2positions <- function(targetbed){
  bed <- fread(input = targetbed,colClasses = list(character=1),data.table = F,stringsAsFactors = F,header = F)
  bed$V2 <- bed$V2+1
  chromosomes <- gsub(unique(bed$V1),pattern = "chr",replacement = "")
  chromosomes <- sort(paste0("chr",chromosomes))
  PosByChrom <- list()

  unwrap <- function(x){
    return(seq.int(from = x[2],to = x[3]))
  }

  for(chr in chromosomes){
    bed_chrom <- bed[grep(bed$V1,pattern = chr),,drop=F]
    PosByChrom[[chr]] <- data.frame(chr=chr,
                                    pos=as.numeric(unlist( apply(bed_chrom, MARGIN=1,FUN = unwrap))),
                                    stringsAsFactors = F)
  }

  return(list(chromosomes=chromosomes, PosByChrom=PosByChrom))

}
