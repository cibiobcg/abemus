#' Compute the per-base error model
#' @param data Output from \code{importData()} function.
#' @param max.vaf.pbem Threshold on VAF.
#' @param cov.min.pbem Threshold on coverage.
#' @return A list with pbem values split by chromosome.
#' @export
getPbem <- function(data,max.vaf.pbem = 0.2,cov.min.pbem = 10){

  out <- lapply(seq_len(length(data$ref)),
                chrom_pbem,
                mat_cov = data$mat_cov,
                mat_vaf = data$mat_vaf,
                mat_cov_base = data$mat_cov_base,
                ref = data$ref,
                max.vaf.pbem = max.vaf.pbem,
                cov.min.pbem = cov.min.pbem)

  return(out)

}
