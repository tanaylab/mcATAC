
#' Cluster atac peaks based on atac distributions
#'
#' @param atac_mc a McATAC object
#' @param peak_set a PeakIntervals object (or a misha intervals set). if NULL - the peaks from \code{atac_mc} would be used.
#' @param k number of clusters
#' @param peak_set - (optional) a subset of peaks of \code{atac_mc@peaks} on which to cluster

#' @inheritDotParams tglkmeans::TGL_kmeans
#' @return a tglkmeans clustering object

#' @examples
#' \dontrun{
#' my_atac_mc <- gen_atac_peak_clust(my_atac_mc, )
#' }
#'
#' @export
gen_atac_peak_clust <- function(atac_mc, k, peak_set = NULL, ...) {
    if (!is.null(peak_set)) {
        atac_mc <- subset_peaks(atac_mc, peak_set)
    }
    return(tglkmeans::TGL_kmeans(as.matrix(atac_mc@mat), k, ...))
}

#' Cluster metacells based on atac profiles using the k-means algorithm
#'
#' @param atac_mc - an McATAC object
#' @param use_prior_annot (optional) when TRUE - use the metacell annotation to generate metacell clusters. Clusters would be generated based on a categorical field \code{annot} from the \code{metadata} slot in the McATAC object.
#' @param k - (optional, when \code{use_prior_annot == F}) number of clusters to generate
#' @param peak_set - a subset of peaks of \code{atac_mc@peaks} on which to cluster
#' @param annot - name of the field to use when \code{use_prior_annot == T}.
#'
#' @inheritParams gen_atac_peak_clust
#' @inheritDotParams tglkmeans::TGL_kmeans
#' @return a named numeric vector specifying the cluster for each metacell
#' @examples
#' \dontrun{
#' ## Use "default clustering" - the existing annotations
#' mc_clusters <- gen_atac_mc_clust(my_atac_mc, use_prior_annot = T)
#'
#' ## Identify peaks of interest, namely peaks neighboring a set of feature genes, and use only them for clustering
#' nei_peaks_feat_genes <- gintervals.neighbors(my_atac_mc@peaks, tss[tss$name %in% feature_genes, ], maxdist = 5e+5)
#' peaks_of_interest <- nei_peaks_feat_genes[, c("chrom", "start", "end")]
#' mc_clusters <- gen_atac_mc_clust(my_atac_mc, k = 16, peak_set = peaks_of_interest, use_prior_annot = F)
#' }
#' @export
gen_atac_mc_clust <- function(atac_mc, use_prior_annot = TRUE, k = NULL, peak_set = NULL, annot = "cell_type", ...) {
    if (!is.null(peak_set)) {
        atac_mc <- subset_peaks(atac_mc, peak_set)
    }
    if (!use_prior_annot) {
        if (!is.null(k)) {
            atac_mc_km <- tglkmeans::TGL_kmeans(t(as.matrix(atac_mc@mat)), k, ...)
            return(setNames(atac_mc_km$cluster, 1:length(atac_mc_km$cluster)))
        } else {
            cli_abort("Must choose k if clustering with use_prior_annot == FALSE")
        }
    } else {
        assert_that(!is.null(atac_mc@metadata), any(grepl(colnames(atac_mc@metadata) == annot)),
            msg = 'There is no metadata or the field "{annot}" does not exist in it.'
        )
        return(deframe(atac_mc@metadata %>% select(metacell, !!annot)))
    }
}

#' Subset peaks of an McATAC object
#'
#' @param atac_mc - an McATAC object
#' @param peak_set - a subset of peaks of \code{atac_mc@peaks} to keep
#'
#' @examples
#' \dontrun{
#' ## Use "default clustering" - the existing annotations
#' mc_clusters <- gen_atac_mc_clust(my_atac_mc, use_prior_annot = T)
#'
#' ## Identify peaks of interest, namely peaks neighboring a set of feature genes, and subset by them
#' nei_peaks_feat_genes <- gintervals.neighbors(my_atac_mc@peaks, tss[tss$name %in% feature_genes, ], maxdist = 5e+5)
#' peaks_of_interest <- nei_peaks_feat_genes[, c("chrom", "start", "end")]
#' mc_clusters <- subset_peaks(my_atac_mc, peaks_of_interest)
#' }
#' @export
subset_peaks <- function(atac_mc, peak_set) {
    atac_mc@peaks$tmpID <- 1:nrow(atac_mc@peaks)
    pks_filt <- dplyr::semi_join(atac_mc@peaks, peak_set, by = c("chrom", "start", "end"))
    atac_mc@peaks <- atac_mc@peaks[pks_filt$tmpID, ]
    atac_mc@mat <- atac_mc@mat[pks_filt$tmpID, ]
    return(atac_mc)
}
