#' get pbem statistics
#'
#' @param dt completetab_dbsnp as data.table
#' @param n_pos_af_th When compute pbem, count in how many germline samples the position has an AF >= n_pos_af_th. default: 0.2
#' @return out data.table
get_stats_pbem <- function(dt,n_pos_af_th){
  ids <- unique(dt$group)
  m <- matrix(nrow = length(ids),ncol = 12,data = 0,dimnames = list(ids))
  for(i in 1:nrow(dt)){
    m[dt$group[i],1] <- m[dt$group[i],1] + dt[i,"RD"]
    m[dt$group[i],2] <- m[dt$group[i],2] + dt[i,"Ade"]
    m[dt$group[i],3] <- m[dt$group[i],3] + dt[i,"Cyt"]
    m[dt$group[i],4] <- m[dt$group[i],4] + dt[i,"Gua"]
    m[dt$group[i],5] <- m[dt$group[i],5] + dt[i,"Thy"]
    if(dt[i,"RD"] > 0){m[dt$group[i],6] <- m[dt$group[i],6] + 1}
    if(dt[i,"af"] < n_pos_af_th){m[dt$group[i],7] <- m[dt$group[i],7] + 1}
    if(dt[i,"af"] >= n_pos_af_th){m[dt$group[i],8] <- m[dt$group[i],8] + 1}
  }

  out <- as.data.frame(m)

  colnames(out) <- c("tot_coverage","total.A","total.C","total.G","total.T","n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth")

  out$group <- rownames(out)
  out$chr <- sapply(strsplit(out$group,":"), `[`, 1)
  out$pos <- sapply(strsplit(out$group,":"), `[`, 2)
  out$ref <- sapply(strsplit(out$group,":"), `[`, 3)
  out$dbsnp <- sapply(strsplit(out$group,":"), `[`, 4)

  out <- out[,c("group","chr","pos","ref","dbsnp","tot_coverage","total.A","total.C","total.G","total.T","n_pos_available","n_pos_af_lth","n_pos_af_gth","count.A_af_gth","count.C_af_gth","count.G_af_gth","count.T_af_gth")]

  row.names(out) <- NULL
  out$group <- paste(out$chr,out$pos,out$ref,sep = ":")

  return(out)
}
