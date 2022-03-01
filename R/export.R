#' Write a McATAC or ScATAC object to an h5ad file
#'
#'
#' @param object McATAC or ScATAC object
#' @param out_file name of the output file
#'
#' @return None.
#'
#' @export
export_to_h5ad <- function(object, out_file) {

}

#' Generate a ucsc genome browser file for each metacell cluster
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust}
#' @param fn name of the genome browser file (should have a '.ucsc' extension)
#'
#' @export
export_atac_clust_ucsc <- function(mc_atac, mc_atac_clust, fn) {

}

#' Generate a misha track for each atac metacell cluster
#'
#' @description generate a track for each metacell cluster, of the form \code{track_prefix.name}, where names
#' are given at \code{clust_names}
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust}
#' @param clust_names names for each metacell cluster
#' @param track_prefix prefix for generated misha tracks.
#'
#' @export
export_atac_clust_misha <- function(mc_atac, mc_atac_clust, clust_names) {

}
