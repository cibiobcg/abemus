pabs_splitbychr <- function(id,CHR,snvs_file){
  chrom = CHR[id]
  filename = gsub(basename(snvs_file),pattern = ".pabs",replacement = paste0("_chr",gsub(chrom,pattern = "chr",replacement = ""),".pabs"))
  cat(c("chr","pos","ref","alt","A","C","G","T","af","cov","Ars","Crs","Grs","Trs","dbsnp\n"),sep = "\t",file = filename)
  awk = paste0("/usr/bin/awk '{if ($1 == \"",chrom,"\") { print $0; }}' ",snvs_file," >> ",filename)
  system(awk)
}
