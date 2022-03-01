
#' Cluster atac peaks based on atac distributions
#'
#' @param atac_mc a McATAC object
#' @param peak_set a PeakInervals object (or a misha intervals set). if NULL - he peaks from \code{atac_mc} would be used.
#' @param k number of clusters
#'
#' @export
gen_atac_peak_clust <- function(atac_mc, k, peak_set = NULL) {

}


#' Cluster metacells based on atac profiles
#'
#' @param use_prior_annot when TRUE - use the metacell annotation to generate metacell clusters. Clusters would be generated based on a categorical field \code{annot} from the \code{metadata} slot in the McATAC object.
#' @param annot name of the field to use when \code{use_prior_annot} is TRUE.
#'
#' @inheritParams gen_atac_peak_clust
#' @export
gen_atac_mc_clust <- function(atac_mc, k, peak_set = NULL, use_prior_annot = FALSE, annot = "cell_type") {

}
