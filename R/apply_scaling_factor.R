apply_scaling_factor <- function(tabindex,
                                 R = 1,
                                 use.optimal.R = FALSE,
                                 target_size = NA){

  if(use.optimal.R){

    #<check>
    if(!"case_mean_coverage"%in%colnames(tabindex)){
      message(paste("[",Sys.time(),"]\tError. column 'case_mean_coverage' must be present in tabindex when use.optimal.R = TRUE"))
      stop()
    }
    if(is.na(target_size)){
      message(paste("[",Sys.time(),"]\tError. 'target_size' needed when use.optimal.R = TRUE"))
      stop()
    }

    tabindex[,paste0("tabcalls_f3","_optimalR")] <- NA
    tabindex[,paste0("tabcalls_f3","_optimalR_used")] <- NA

    for(i in 1:nrow(tabindex)){
      a <- fread(input = tabindex$tabcalls_f3[i],stringsAsFactors = FALSE,data.table = FALSE)
      # select best R for this case sample

      sizes <- sort(unique(tab_optimal_R$target_Mbp))
      closest_target_size <- sizes[which(abs(sizes-target_size)==min(abs(sizes-target_size)))]

      bombanel <- tab_optimal_R[which(tab_optimal_R$target_Mbp == closest_target_size),]

      covs <- bombanel$mean_coverage
      closest_coverage <- covs[which(abs(covs-tabindex$case_mean_coverage[i])==min(abs(covs-tabindex$case_mean_coverage[i])))]

      optR <- bombanel$scalingfactor[which(bombanel$mean_coverage==closest_coverage)]

      a$filter.pbem_coverage <- a$filter.pbem_coverage * optR
      a$pass.filter.pbem_coverage <- 0
      a$pass.filter.pbem_coverage[which(a$af_case >= a$filter.pbem_coverage)] <- 1
      a <- a[which(a$pass.filter.pbem_coverage==1),]
      out.name <- gsub(basename(tabindex$tabcalls_f3[i]),pattern = "pmtab_F3_",replacement = paste0("pmtab_F3_optimalR_"))
      out.path <- gsub(tabindex$tabcalls_f3[i],pattern = basename(tabindex$tabcalls_f3[i]),replacement = out.name)
      #write.table(x = a,file = out.path,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
      tabindex[i,paste0("tabcalls_f3","_optimalR")] <- out.path
      tabindex[i,paste0("tabcalls_f3","_optimalR_used")] <- optR
    }
    message(paste("[",Sys.time(),"]\talright."))
    return(tabindex)

  } else {

    tabindex[,paste0("tabcalls_f3","_R",R)] <- NA
    for(i in 1:nrow(tabindex)){
      a <- fread(input = tabindex$tabcalls_f3[i],stringsAsFactors = FALSE,data.table = FALSE)
      a$filter.pbem_coverage <- a$filter.pbem_coverage * R
      a$pass.filter.pbem_coverage <- 0
      a$pass.filter.pbem_coverage[which(a$af_case >= a$filter.pbem_coverage)] <- 1
      a <- a[which(a$pass.filter.pbem_coverage==1),]
      out.name <- gsub(basename(tabindex$tabcalls_f3[i]),pattern = "pmtab_F3_",replacement = paste0("pmtab_F3_R",R,"_"))
      out.path <- gsub(tabindex$tabcalls_f3[i],pattern = basename(tabindex$tabcalls_f3[i]),replacement = out.name)
      #write.table(x = a,file = out.path,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
      tabindex[i,paste0("tabcalls_f3","_R",R)] <- out.path
    }
    message(paste("[",Sys.time(),"]\talright."))
    return(tabindex)

  }
}
