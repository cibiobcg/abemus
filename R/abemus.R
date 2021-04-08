#' ABEMUS:  Adaptive per Base Error Model in Ultra-deep Sequencing Data
#'
#' The package abemus is a NGS-based computational method that uses control samples to build global and local sequencing error reference models that are then used to improve the detection of somatic SNVs in cfDNA samples.
#'
#' In step_1 and step_2 only control samples are used to build the per-base error model and compute both coverage-dependent and coverage-independent allelic fraction thresholds.\cr
#' In step_3 tumor samples are inspected and somatic snvs are called based on both custom filtering criteria (i.e. minimum coverage, minimum number of reads supporting an alternative allele) and allelic fraction thresholds.\cr
#' Then in step_4, the per-base error model computed in step 1 is exploited to further refine the final set of putative somatic snvs.
#
#' @docType package
#' @name abemus
#' @import utils
#' @import stats
#' @import parallel
NULL

globalVariables(c("fread",
                  "afz",
                  "bgpbem",
                  "bombanel_afs",
                  "bombanel_covs",
                  "bombanel_tab_cov_pbem",
                  "datacount_bin",
                  "tab_optimal_R",
                  "th_results_bin"))
