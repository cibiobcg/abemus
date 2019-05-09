pileup_splitbychr <- function(id,CHR,pileup_file){
  chrom = CHR[id]
  filename = gsub(basename(pileup_file),pattern = ".pileup",replacement = paste0("_chr",gsub(chrom,pattern = "chr",replacement = ""),".pileup"))
  cat(c("chr","pos","ref","A","C","G","T","af","cov","dbsnp\n"),sep = "\t",file = filename)
  awk = paste0("/usr/bin/awk '{if ($1 == \"",chrom,"\") { print $0; }}' ",pileup_file," >> ",filename)
  system(awk)
}
