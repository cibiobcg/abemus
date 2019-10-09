#' Compute the mean coverage of a case sample
#' @param tabindex output calls$tabsnvs_index from callsnvs() function.
#' @param pacbamfolder_bychrom folder with pileups split by chromosome.
#' @return sample.info.file with 1 additional column reporting the mean coverage for each case having a matched-control
#' @export
get_case_mean_coverage <- function(tabindex,
                                   pacbamfolder_bychrom){

  pacbam_list <- list.files(pacbamfolder_bychrom,full.names = TRUE)

  chrom_mean_cov <- function(cp){
    chrom_cov <- fread(input = cp,sep = "\t",skip = 1,stringsAsFactors = FALSE,select = 9,data.table = FALSE)
    return( mean(as.numeric(chrom_cov[,1])) )
  }

  get_pacbam_case_pileups <- function(pat, pacbam_list){
    pacbam_case_folder <- pacbam_list[grep(pattern = paste("^",pat,"$",sep=""),x = basename(pacbam_list))]
    pacbam_case_pileups <- list.files(file.path(pacbam_case_folder,"pileup"),full.names = TRUE)
    mean_covs <- as.numeric(sapply(pacbam_case_pileups,FUN = chrom_mean_cov))
    return( mean(mean_covs) )
  }

  tabindex$case_mean_coverage <- as.numeric(sapply(tabindex[,2], get_pacbam_case_pileups, pacbam_list))

  return( tabindex )
}
