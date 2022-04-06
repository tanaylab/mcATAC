
#' Cluster atac peaks based on atac distributions
#'
#' @param atac_mc a McATAC object
#' @param k number of clusters; must be specified if \code{clustering_algoritm = 'kmeans'}
#' @param clustering_algoritm (optional) either "kmeans" or "louvain"
#' @param cluster_on (optional; default - fp) which matrix (\code{mat}/\code{fp}/\code{egc})to cluster on
#'
#' @inheritDotParams tglkmeans::TGL_kmeans
#' @return a named numeric vector specifying the cluster for each peak

#' @examples
#' \dontrun{
#' my_atac_mc <- gen_atac_peak_clust(my_atac_mc, k = 16, cluster_on = "mat")
#'
#' dyn_p <- identify_dynamic_peaks(my_atac_mc)
#' my_atac_mc <- gen_atac_peak_clust(my_atac_mc, k = 16, cluster_on = "fp", peak_set = dyn_p)
#' }
#'
#' @export
gen_atac_peak_clust <- function(atac_mc, k = NULL, clustering_algoritm = "kmeans", cluster_on = "fp", ...) {
    louvain_k <- 5
    if (cluster_on %!in% c("fp", "mat", "egc")) {
        cli_abort("{.var cluster_on} must be either 'fp', 'mat' or 'egc'")
    }
    assert_atac_object(atac_mc)
    if (clustering_algoritm == "kmeans" && is.null(k)) {
        cli_abort("Specify {.var k} when clustering with kmeans")
    }
    if (clustering_algoritm == "louvain") {
        cli_li("Clustering using {.val louvain}")
        if (is.null(k)) {
            k <- louvain_k
        }
        mca_knn <- tgs_cor_knn(x = t(atac_mc@fp), y = t(atac_mc@fp), knn = k, spearman = T)
        gknn <- igraph::graph_from_data_frame(mca_knn[, c("col1", "col2")], directed = F)
        louv_cl <- igraph::cluster_louvain(graph = gknn)
        atac_peak_cl <- setNames(louv_cl$membership, rownames(atac_mc@mat))
    } else {
        cli_li("Clustering using {.val kmeans++}. k = {.val {k}}")
        atac_peak_km <- tglkmeans::TGL_kmeans(as.matrix(slot(atac_mc, cluster_on)), k, id_column = FALSE, ...)
        atac_peak_cl <- setNames(atac_peak_km$cluster, rownames(atac_mc@mat))
    }
    return(atac_peak_cl)
}

#' Cluster metacells based on atac profiles using the k-means algorithm
#'
#' @param atac_mc - an McATAC object
#' @param use_prior_annot (optional) when TRUE - use the metacell annotation to generate metacell clusters. Clusters would be generated based on a categorical field \code{annot} from the \code{metadata} slot in the McATAC object.
#' @param k - (optional, when \code{use_prior_annot == F}) number of clusters to generate
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
gen_atac_mc_clust <- function(atac_mc, use_prior_annot = TRUE, k = NULL, annot = "cell_type", ...) {
    assert_atac_object(atac_mc)

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
#' @return the atac_mc object only with the peaks of interest (not saved in the "ignore_..." slots)
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
    assert_atac_object(atac_mc)
    pks_filt <- semi_join(atac_mc@peaks, peak_set, by = c("chrom", "start", "end", "peak_name"))
    atac_mc@peaks <- atac_mc@peaks[atac_mc@peaks$peak_name %in% pks_filt$peak_name, ]
    atac_mc@mat <- atac_mc@mat[rownames(atac_mc@mat) %in% pks_filt$peak_name, ]
    return(atac_mc)
}
