#' apply_scaling_factor
#' @export
#' @param tabindex output from callsnvs()
#' @param R scaling factor to be applied at filter.pbem_coverage value. default = 1
apply_scaling_factor <- function(tabindex,R=1){
  tabindex[,paste0("tabcalls_f3","_R",R)] <- NA
  for(i in 1:nrow(tabindex)){
    a <- fread(input = tabindex$tabcalls_f3[i],stringsAsFactors = F,data.table = F)
    a$filter.pbem_coverage <- a$filter.pbem_coverage * R
    a$pass.filter.pbem_coverage <- 0
    a$pass.filter.pbem_coverage[which(a$af_case >= a$filter.pbem_coverage)] <- 1
    a <- a[which(a$pass.filter.pbem_coverage==1),]
    out.name <- gsub(basename(tabindex$tabcalls_f3[i]),pattern = "pmtab_F3_",replacement = paste0("pmtab_F3_R",R,"_"))
    out.path <- gsub(tabindex$tabcalls_f3[i],pattern = basename(tabindex$tabcalls_f3[i]),replacement = out.name)
    write.table(x = a,file = out.path,quote = F,sep = "\t",row.names = F,col.names = T)
    tabindex[i,paste0("tabcalls_f3","_R",R)] <- out.path
  }
  cat(paste("[",Sys.time(),"]\talright.","\n"))
  return(tabindex)
}
