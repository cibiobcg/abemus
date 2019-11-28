#' Compute the mean coverage of a CASE sample
#' @param tabindex The data.frame output by the \code{callsnvs} function.
#' @param pacbamfolder_bychrom The folder popluted by outputs by the \code{split_pacbam_bychrom} function.
#' @return The \code{sample.info.file} with an additional column reporting the mean coverage for each case (only for those with a matched-CONTROL)
#' @export
#' @examples
#' sample.info.file <- system.file("extdata", "test_sif_toy.tsv", package = "abemus")
#' outdir <- tempdir()
#' targetbed <- system.file("extdata", "regions_toy.bed", package = "abemus")
#' pacbamfolder_bychrom <- system.file("extdata", "pacbam_data_bychrom", package = "abemus")
#' pbem_dir <- system.file("extdata", "BaseErrorModel", package = "abemus")
#' controls_dir <- system.file("extdata", "Controls", package = "abemus")
#' m<-callsnvs(sample.info.file,outdir,targetbed,pbem_dir,controls_dir,pacbamfolder_bychrom,replicas=1)
#' tabindex <- m$tabsnvs_index
#' tabindex <- get_case_mean_coverage(tabindex = tabindex,pacbamfolder_bychrom = pacbamfolder_bychrom)
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
