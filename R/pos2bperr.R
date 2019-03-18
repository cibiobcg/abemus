#' pos2bperr
#'
#' Compute the per-base error model for each position
#' @param id
#' @param targets
#' @param germlineset
#' @param step
#' @param chrom
#' @param lev
#' @param covbin
#' @param af_max_to_compute_thresholds
#' @param coverage_min_to_compute_thresholds
#' @param af_max_to_compute_pbem
#' @param coverage_min_to_compute_pbem
#' @param n_pos_af_th
#'
pos2bperr = function(id,targets,germlineset,step,chrom,lev,covbin,af_max_to_compute_thresholds,coverage_min_to_compute_thresholds,af_max_to_compute_pbem,coverage_min_to_compute_pbem,n_pos_af_th){
  upto = id+step-1
  if(upto>nrow(targets)){upto <- nrow(targets)}
  this = targets[id:upto,,drop=F]
  outfile = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"postogrep.txt",sep="_")
  filtpileup = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"filtered.pileup.txt",sep="_")
  taboutchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"pbem.table.txt",sep="_")
  mytabafchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"afgtz.table.txt",sep="_")
  afzchrom = paste(this$chr[1],this$pos[1],this$pos[nrow(this)],"afz.table.txt",sep="_")
  afz = array(data = 0,dim = length(lev),dimnames = list(lev))
  cat(unique(this$pos),sep = "\n",file = outfile,append = F)
  cmd = paste0("awk -F'\t' '{if (FILENAME == \"",outfile,"\") { t[$1] = 1; } else { if (t[$2]) { print }}}' ",outfile," ",germlineset," > ",filtpileup)
  system(cmd)
  if(file.info(filtpileup)$size == 0){
    cat(paste("[",Sys.time(),"]\tpositions in ",outfile,"not found in any pileups."),"\n")
  } else {
    completetab_all = read.delim(filtpileup,stringsAsFactors = F,header = F,sep = "\t",na.strings = "")
    names(completetab_all)=c("chr","pos","ref","A","C","G","T","af","RD","dbsnp")
    # exclude annotated and private SNPs [ to compute AF threshold]
    completetab = completetab_all[which(is.na(completetab_all$dbsnp)),,drop=F]
    completetab = completetab[which(completetab$af <= af_max_to_compute_thresholds & completetab$RD >= coverage_min_to_compute_thresholds),,drop=F]
    # Compute pbem also in positions that are SNPs
    completetab_dbsnp = completetab_all[which(completetab_all$af <= af_max_to_compute_pbem & completetab_all$RD >= coverage_min_to_compute_pbem),,drop=F]
    if(nrow(completetab)>0){
      # save allelic fractions by bins of coverage
      mytabaf = completetab[which(completetab$af > 0),,drop=F]
      mytabz =  completetab[which(completetab$af == 0),,drop=F]
      afz = afz + as.array(table(cut(mytabz$RD,breaks = covbin,include.lowest = T)))
      write.table(t(afz),file = afzchrom,sep="\t",col.names = F,row.names = F,quote = F,append = F)
      write.table(mytabaf[,8:9],file = mytabafchrom,append = F,sep = "\t",quote = F,row.names = F,col.names = F)
      # compute pbem
      completetab_dbsnp$group = paste(completetab_dbsnp$chr,completetab_dbsnp$pos,completetab_dbsnp$ref,sep=":")
      dat <- data.table(completetab_dbsnp, key="group")
      ans <- dat[,{list(tot_coverage=sum(RD,na.rm = T),
                        total.A=sum(A,na.rm = T),
                        total.C=sum(C,na.rm = T),
                        total.G=sum(G,na.rm = T),
                        total.T=sum(T,na.rm = T),
                        n_pos_available = length(which(RD > 0)),
                        n_pos_af_lth=length(which(af < n_pos_af_th)),
                        n_pos_af_gth=length(which(af >= n_pos_af_th)),
                        count.A_af_gth=length(A[which(af >= n_pos_af_th & A>0)]),
                        count.C_af_gth=length(C[which(af >= n_pos_af_th & C>0)]),
                        count.G_af_gth=length(G[which(af >= n_pos_af_th & G>0)]),
                        count.T_af_gth=length(T[which(af >= n_pos_af_th & T>0)]))},by="group"]
      this = as.data.table(this)
      this$group = paste(this$chr,this$pos,this$ref,sep=":")
      this_ans = merge(x = this,y = ans,by = "group",all.y = T)
      this_ans = as.data.frame(this_ans)
      tabstats = fromListToDF(mclapply(seq(1,nrow(this_ans),1),perbaserror_stats,tgf=this_ans,mc.cores = 2 ))
      tabstats = tabstats[with(tabstats, order(pos)), ]
      cat(paste("[",Sys.time(),"]\tWriting output for positions in: ",filtpileup),"\n")
      write.table(tabstats,file = taboutchrom,append = F,quote = F,row.names = F,col.names = F,sep="\t")
    }
  }
  system(paste("rm",filtpileup,outfile))
}
