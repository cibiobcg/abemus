#' adjust_af_threshold_table
#'
#' @param controls_dir folder with afth data. default: file.path(outdir,"Controls")
#' @param pbem_dir folder with pbem data. default: file.path(outdir,"BaseErrorModel")
#' @param detection.specificity default: 0.995
#' @param replicas default: 1000
#' @param replicas.in.parallel default: 1
#' @param coverage_binning Bins of coverage into which divide allelic fractions. default: 50
#' @param coeffvar.threshold Consider a bin as stable if the Coeff of variations after n replicas is lower than coeffvar.threshold. default: 0.01
#' @return corrected af threshold table
adjust_af_threshold_table <- function(controls_dir,
                                      pbem_dir,
                                      detection.specificity,
                                      replicas,
                                      replicas.in.parallel,
                                      coverage_binning,
                                      coeffvar.threshold){

  dthdataFile <- file.path(controls_dir,"datathreshold.RData")
  if(file.exists(dthdataFile)){
    cat(paste("[",Sys.time(),"]\tlooking for",dthdataFile,"[ ok ]"),"\n")
    load(file = dthdataFile)
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for",dthdataFile,"[ not found ]"),"\n")
    stop()
  }

  vafcov_file = file.path(pbem_dir,"afgtz.tsv")
  if(file.exists(vafcov_file)){
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ ok ]"),"\n")
    #vafcov = read.big.matrix(filename = vafcov_file,header = F,sep = "\t",type = "double")
    vafcov = fread(input = vafcov_file,sep = "\t",header = F,stringsAsFactors = F,data.table = F)
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs > 0 and coverages:",vafcov_file,"[ not found ]"),"\n")
    stop()
  }

  afz_file = file.path(pbem_dir,"afz.RData")
  if(file.exists(afz_file)){
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ ok ]"),"\n")
    load(afz_file)
  } else {
    cat(paste("[",Sys.time(),"]\tlooking for data.table with AFs = 0 and coverages:",afz_file,"[ not found ]"),"\n")
    stop()
  }

  # Find out most represented bin of coverage
  datacount_bin_complete <-  datacount_bin+afz
  datacount_bin_complete <-  datacount_bin_complete[grep(pattern = "Inf",names(datacount_bin_complete),invert = T,value = T)]
  #datacount_bin_complete <-  sort(datacount_bin_complete,decreasing = T) # sort by all
  datacount_bin_complete <- datacount_bin_complete[names(sort(datacount_bin,decreasing = T))] # sort by AF>0
  datacount_bin_complete <- datacount_bin_complete[which(datacount_bin_complete > 0)]

  stop = length(datacount_bin_complete)-1
  last.stable.card = 0
  tab = c()

  covbin <- define_cov_bins(coverage_binning)[[1]]
  lev <- define_cov_bins(coverage_binning)[[2]]
  a <- cut(x = as.integer(vafcov[,2]),breaks = covbin,include.lowest = T)

  for(i in 1:stop){
    current.bin = names(datacount_bin_complete)[i]
    current.bin.card = as.numeric(datacount_bin_complete[i])
    cat(paste("[",Sys.time(),"]\tevaluating bin:",current.bin),"\n")
    if(current.bin.card >= last.stable.card){
      for(j in (i+1):(stop+1)){
        next.bin = names(datacount_bin_complete)[j]
        next.bin.card = as.numeric(datacount_bin_complete[j])
        next.bin.AFetz = as.numeric(afz[next.bin])
        next.bin.AFgtz = as.numeric(datacount_bin[next.bin])
        ReplicasTable = c()
        idx = which(a==current.bin)
        out = mclapply(1:replicas,
                       get_afth,
                       idx=idx,
                       next.bin.AFgtz=next.bin.AFgtz,
                       next.bin.AFetz=next.bin.AFetz,
                       vafcov=vafcov,
                       this.detection.specificity=detection.specificity,
                       mc.cores = replicas.in.parallel)
        ReplicasTable = matrix(unlist(out),ncol = replicas)
        coeffvar = apply(ReplicasTable,MARGIN = 1,FUN = sd,na.rm=T)/apply(ReplicasTable,MARGIN = 1,FUN = mean,na.rm=T)
        if(is.na(coeffvar)){
          coeffvar <- 0
        }
        x = data.frame(current.bin=current.bin,next.bin=next.bin,coeffvar=as.numeric(coeffvar),median.afth=median(ReplicasTable,na.rm = T),last.stable.card=last.stable.card,stringsAsFactors = F)
        tab=rbind(tab,x)
        if(coeffvar > coeffvar.threshold){
          last.stable.card = max(last.stable.card,next.bin.card)
          break
        }
      }
    }
  }
  name.out = paste0("tab",replicas,"r_",detection.specificity,"_b.RData")
  save(tab,file = file.path(controls_dir, name.out),compress = T)

  # correct the original threhsold table
  cat(paste("[",Sys.time(),"]\tcorrecting original afth table"),"\n")
  minaf_cov_corrected <- th_results_bin[which(round(th_results_bin$specificity,4) == as.numeric(detection.specificity)),]
  C = tab$last.stable.card[nrow(tab)]

  last.afth.used = 1
  start = 2
  end = length(minaf_cov_corrected)-1

  for(i in start:end){
    bin.name = names(minaf_cov_corrected)[i]
    bin.card = as.numeric(datacount_bin_complete[which(names(datacount_bin_complete)==bin.name)])
    if(!identical(bin.card,numeric(0))){
      if(bin.card < C & last.afth.used == 1){
        minaf_cov_corrected[i] <- last.afth.used
      }
      if(bin.card >= C){
        last.afth.used <- minaf_cov_corrected[i]
      }
      if(bin.card < C & last.afth.used != 1){
        minaf_cov_corrected[i] <- last.afth.used
      }
    }
  }

  minaf_cov_corrected[which(is.na(minaf_cov_corrected))] <- last.afth.used

  #minaf_cov_corrected[which(minaf_cov_corrected==1)] <- NA #
  minaf_cov_corrected[which(minaf_cov_corrected==1)] <- as.numeric(minaf_cov_corrected[which(minaf_cov_corrected!=1)[2]])

  # correct the AF threhsold in bin where there is Inf as limit
  N = length(minaf_cov_corrected)
  minaf_cov_corrected[N] <- minaf_cov_corrected[N-1]

  cat(paste("[",Sys.time(),"]\tsaving corrected table"),"\n")
  name.out = paste0("minaf_cov_corrected_",replicas,"r_",detection.specificity,"_b.RData") # _b annotation, to be removed
  save(minaf_cov_corrected,file = file.path(controls_dir,name.out))

  cat(paste("[",Sys.time(),"]\talright.","\n"))
  return(list(info.tab=tab,
              minaf_cov_corrected=minaf_cov_corrected))
}
