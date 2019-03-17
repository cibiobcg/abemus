#' Compute pbem
#'
#' Compute per-base error model on each targeted position and save AFs by bins of coverage
#' @export
#' @param chromosomes vector of unique chromosomes to be analysed, i.e. c("chr7","chr21","chr13")
#' @param sif the output [[1]] from import_sif
#' @param targetbp folder with RData for each annotated positions
#' @param step
#' @param bam.with.chr
compute_pbem <- function(chromosomes,sif,step=5000,bam.with.chr=FALSE){

  for(chrom in chromosomes){
    cat(paste("[",Sys.time(),"]\tchromosome:",chrom),"\n")
    germlineset = get_germlineset(sif)

    tp = list.files(targetbp,pattern = paste0(chrom,"\\.RData$"),full.names = T)
    load(tp,verbose = F)
    chromTargetPositions = as.data.table(chromTargetPositions,keep.rownames = F)

    mytargets = unique(chromTargetPositions)
    targetsFILT <- mytargets
    #targetsFILT$randompos <- 1 # deprecated, to remove

    # Final targets
    targets = unique(targetsFILT)
    targets = targets[with(targets,order(chr,pos)),]
    targets = as.data.frame(targets)

    if(bam.with.chr){
      targets$chr = paste0('chr',targets$chr)
    }

    cat(paste("[",Sys.time(),"]\tTotal positions to check in this chromosome :",nrow(targets)),"\n")
    #step = 5000
    mclapply(seq(1,nrow(targets),step),pos2bperr,
             targets=targets,
             germlineset=germlineset,
             step=step,
             chrom=chrom,
             lev=lev,
             covbin=covbin,
             af_max_to_compute_thresholds=af_max_to_compute_thresholds,
             coverage_min_to_compute_thresholds=coverage_min_to_compute_thresholds,
             af_max_to_compute_pbem=af_max_to_compute_pbem,
             coverage_min_to_compute_pbem=coverage_min_to_compute_pbem,
             n_pos_af_th=n_pos_af_th,
             mc.cores=mc.cores)

    cmda = paste("cat *_pbem.table.txt >",paste0("bperr_",chrom,".tsv"))
    system(cmda)
    system("rm *_pbem.table.txt")
    cmdb = paste("cat *_afgtz.table.txt >",paste0("afgtz_",chrom,".tsv"))
    system(cmdb)
    system("rm *_afgtz.table.txt")
    cmdc = paste("cat *_afz.table.txt >",paste0("afz_",chrom,".tsv"))
    system(cmdc)
    system("rm *_afz.table.txt")

  }
  cmd.merge = paste("cat afgtz_chr*.tsv > afgtz.tsv")
  system(cmd.merge)
  cmd.merge = paste("cat afz_chr*.tsv > afz.tsv")
  system(cmd.merge)

  # save counter of afz as RData
  afztab = read.delim(file = "afz.tsv",sep="\t",as.is = T,header=F)
  afz = apply(afztab,2,sum)
  names(afz)=lev
  save(afz,file = "afz.RData",compress = T)

  # overall statistics on pbems
  cmd.merge = paste("cat bperr_chr*.tsv > bperr.tsv")
  system(cmd.merge)

  bperr = fread(file.path(outdir, "BaseErrorModel","bperr.tsv"),stringsAsFactors = F,showProgress = F,header = F,colClasses = list(character=2,character=5))

  # summary stats for pbem across the target
  bperr_summary = summary(bperr$V22)
  names = names(bperr_summary)
  bperr_summary = data.frame(as.numeric(bperr_summary))
  bperr_summary = rbind(bperr_summary,sd(x = bperr$V16,na.rm = T))
  rownames(bperr_summary) = c(names,"std")
  write.table(bperr_summary,file = file.path(outdir, "BaseErrorModel","bperr_summary.tsv"),row.names = T,col.names = F,quote = F,sep = "\t")

  pbem.nas = which(is.na(bperr$V22))
  if(length(pbem.nas)>0){
    warning("NOT ABLE TO COMPUTE PBEM IN ",length(pbem.nas)," POSITIONS.")
  }

  # compute background pbem
  bperr_subset = bperr[which(bperr$V17 == 0),]
  bgpbem = (sum(as.numeric(bperr_subset$V23)))/(sum(as.numeric(bperr_subset$V10)))
  mean_pbem = mean(as.numeric(bperr_subset$V22),na.rm = T)
  tabstat = data.frame(background_pbem = bgpbem,
                       mean_pbem = mean_pbem,
                       stringsAsFactors = F)
  write.table(tabstat,file = file.path(outdir, "BaseErrorModel","pbem_background.tsv"),row.names = F,col.names = T,quote = F,sep = "\t")
}

