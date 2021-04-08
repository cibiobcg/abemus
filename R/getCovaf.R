#' Compute the per-base error model
#' @param data The output from \code{importData()} function.
#' @param max.vaf Threshold on VAF. default: 0.2
#' @param min.cov Threshold on coverage. default: 10
#' @return A list with coverage (cov) and VAF (vaf) values split by chromosome.
#' @export
getCovaf <- function(data,max.vaf=0.2,min.cov=10){

  out <- lapply(seq_len(length(data$ref)),
                chrom_cov_vaf,
                mat_cov = data$mat_cov,
                mat_vaf = data$mat_vaf,
                rsid = data$rsid,
                max.vaf = max.vaf,
                min.cov = min.cov)

  return(out)

}
