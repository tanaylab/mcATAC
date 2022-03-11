
#' Cluster atac peaks based on atac distributions
#'
#' @param atac_mc a McATAC object
#' @param peak_set a PeakIntervals object (or a misha intervals set). if NULL - the peaks from \code{atac_mc} would be used.
#' @param k number of clusters
#'
#' @inheritDotParams tglkmeans::TGL_kmeans
#' @return atac_mc with added cluster_k_{k} column to $peaks specifying the cluster for each enhancer
#' @export
gen_atac_peak_clust <- function(atac_mc, k, peak_set = NULL) {
    atac_peak_km <- tglkmeans::TGL_kmeans(as.matrix(atac_mc$mat), k)
    atac_mc$peaks[, glue::glue("cluster_k={k}")] <- atac_peak_km$cluster
    return(atac_mc)
}

#' Cluster metacells based on atac profiles
#'
#' @param use_prior_annot when TRUE - use the metacell annotation to generate metacell clusters. Clusters would be generated based on a categorical field \code{annot} from the \code{metadata} slot in the McATAC object.
#' @param annot name of the field to use when \code{use_prior_annot} is TRUE.
#'
#' @inheritParams gen_atac_peak_clust
#' @inheritDotParams tglkmeans::TGL_kmeans
#' @return a named numeric vector specifying the cluster for each metacell
#' @export
gen_atac_mc_clust <- function(atac_mc, k = NULL, peak_set = NULL, use_prior_annot = TRUE, annot = "cell_type") {
    if (!use_prior_annot) {
        if (!is.null(k)) {
            atac_mc_km <- tglkmeans::TGL_kmeans(as.matrix(atac_mc$mat), k)
            return(setNames(atac_mc_km$cluster, 1:length(atac_mc_km$cluster)))
        }
        else {cli_abort('Must choose k if clustering with use_prior_annot == FALSE')}
    }
    else {
        assert_that(!is.null(atac_mc$metadata), any(grepl(colnames(atac_mc$metadata) == annot)), 
                    msg = 'There is no metadata or the field "{annot}" does not exist in it.')
        return(deframe(atac_mc$metadata %>% select(metacell, !!annot)))
    }
}
