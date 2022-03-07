
#' Cluster atac peaks based on atac distributions
#'
#' @param atac_mc a McATAC object
#' @param peak_set a PeakIntervals object (or a misha intervals set). if NULL - the peaks from \code{atac_mc} would be used.
#' @param k number of clusters
#'
#' @export
gen_atac_peak_clust <- function(atac_mc, k, peak_set = NULL) {
    atac_peak_km = tglkmeans::TGL_kmeans(as.matrix(atac_mc$mat), k)
    atac_mc$peaks[,glue::glue('cluster_k={k}')] = atac_peak_km$cluster
    return(atac_mc)
}

##### Why do we need this function? Manifold structure across metacells should come from RNA
##### -YSh

#' Cluster metacells based on atac profiles
#'
#' @param use_prior_annot when TRUE - use the metacell annotation to generate metacell clusters. Clusters would be generated based on a categorical field \code{annot} from the \code{metadata} slot in the McATAC object.
#' @param annot name of the field to use when \code{use_prior_annot} is TRUE.
#'
#' @inheritParams gen_atac_peak_clust
#' @export
gen_atac_mc_clust <- function(atac_mc, k=NULL, peak_set = NULL, use_prior_annot = TRUE, annot = "cell_type") {
    if (!use_prior_annot) {
        if (!is.null(k)) {
            atac_mc_km = tglkmeans::TGL_kmeans(as.matrix(atac_mc$mat), k)
            return(setNames(atac_mc_km$cluster, 1:length(atac_mc_km$cluster)))
        }
        else {stop('Must choose k if clustering with use_prior_annot == FALSE')}
    }
    else {
        res = match(unlist(atac_mc$metadata[,annot]), sort(unique(unlist(atac_mc$metadata[,annot]))))
        return(setNames(res, 1:length(res)))
    }
}
