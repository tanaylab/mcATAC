#' Find dynamic peaks in McATAC matrix
#'
#' @description This function identifies "dynamic" peaks, i.e. those that have high expression only in a subset of the cells. They are identified by overdispersion in the coefficient of variation (std.dev./mean) per quantiles.
#' @param mcatac the McATAC object to analyze
#' @param method (optional) either 'bmq' (default) or 'gmm'; 'bmq' (binned-mean quantiles) bins the log-mean of all peaks (averaged across metacells) and
#'   selects all enhancers with a coefficient of variation above some quantile in each bin. More controlled 
#'  'gmm' fits a Gaussian mixture model to the log10(COV) vs. log10(mean) distribution, 
#'   and selects peaks in clusters that show overdispersion in the COV.
#' @param mean_thresh_q (optional) threshold quantile on peaks' mean 
#' @param cov_q_thresh (optional) threshold on minimum COV quantile to consider as dynamic in each bin
#' @param num_bins (optional) number of bins to divide features' means into
#' @param gmm_g (optional) number of groups for 'gmm'
#'
#' @return dynamic_peaks - a PeakIntervals object with peaks identified as dynamic
#' @examples
#' \dontrun{
#' my_intervals <- annotate_intervals(my_intervals, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
#' table(my_intervals$peak_annot)
#' }
#' @export
identify_dynamic_peaks <- function(mcatac, method = 'bmq', mean_thresh_q = 0.1, cov_q_thresh = 0.75, num_bins = 200) {
    mns <- rowMeans(mcatac)
    sds <- apply(mcatac, 1, sd)
    covs <- sds/mns
    mn_thresh <- quantile(mns, mean_thresh_q)
    lmi <- lm(cov ~ mn, data <- as.data.frame(list('mn' <- log10(mns), 'cov' <- log10(covs))))
    pred_df <- as.data.frame(list('pred' <- lmi$coefficients[[1]] + lmi$coefficients[[2]]*log10(mns), 'cov' <- log10(covs)))
    pred_df <- dplyr::mutate(pred_df, diff = cov - pred)
    if (method == 'bmq') {
        lpm <- log10(mns)
        qcut <- cut(lpm, breaks <- seq(min(lpm), max(lpm), l=num_bins))
        feat_select <- unlist(sapply(sort(unique(qcut[which(lpm >= log10(mn_thresh))])),
                             function(nm) rownames(mcatac)[which(covs >= quantile(covs[qcut == nm], cov_q_thresh, na.rm=T) & qcut == nm)]))
        pred_df$is_bmq <- rownames(pred_df) %in% feat_select
        smoothScatter(lpm, log10(covs), cex = 0.05)
        points(lpm[feat_select], log10(covs[feat_select]), cex = 0.05, col = 'red')
        return(mcatac@peaks[mcatac@peaks$peak_name %in% feat_select,])
    }
    else if (method == 'gmm') {
        X <- cbind(lpm, log10(covs))
        fit_gmm <- mclust::Mclust(X, G=gmm_g, model="VVV")
        pred_df$fc <- fit_gmm$classification
        cluster_mean_diff <- tapply(pred_df$diff, pred_df$fc, mean)
        clusters_selected <- names(cluster_mean_diff)[cluster_mean_diff > 0]
        return(mcatac@peaks[mcatac@peaks$peak_name %in% rownames(pred_df)[pred_df$fc %in% clusters_selected],])
    }
}

#' Find overlaps with ENCODE blacklists
#'
#' @description Identify peaks in the data which overlap (or are adjacent to?) regions blacklisted by ENCODE as having universally high DNAse HS or ChIP signal (basically mapping artifacts)
#' @param scatac (optional) an ScATAC object
#' @param mcatac (optional) an McATAC object
#' @param peaks (optional) the intervals set to check
#' @param genome (optional, required if checking peaks directly) the interval set to check
#'
#' @return blacklist_overlaps - a PeakIntervals object with peaks identified as overlapping blacklisted regions
#' @examples
#' \dontrun{
#' my_intervals <- annotate_intervals(my_intervals, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
#' table(my_intervals$peak_annot)
#' }
#' @export
find_blacklist_overlaps = function(scatac = NULL, mcatac = NULL, peaks = NULL, genome = NULL) {
    if (!is.null(scatac)) {
        peaks <- scatac@peaks
        genome <- scatac@genome
    }
    else if (!is.null(mcatac)) {
        peaks <- mcatac@peaks
        genome <- mcatac@genome
    }
    else if (is.null(peaks)) {
        cli_abort('Must specify either scatac, mcatac or peaks')
    }
    else if (!is.null(peaks) && is.null(genome)) {
        cli_abort('Must specify genome if annotating peaks directly')
    }
    misha.ext::gset_genome(genome)
    blacklist_name <- glue::glue('ENCODE_blacklist_{genome}'))
    if (gintervals.exists() {
        blacklist <- gintervals.load(blacklist_name)
    }
    nei_blk_pks = gintervals.neighbors(blacklist, peaks, maxdist = 0, mindist = 0, maxneighbors = 10)
    blacklist_overlaps = peaks[unique(nei_blk_pks$intervalID),]
    return(blacklist_overlaps)
}