#' Filter features by summary statistics
#'
#' @description Remove ATAC peaks with low coverage, or that are too long etc.
#' @param scatac the McATAC object to analyze
#' @param minimal_max_umi (optional) threshold on minimum maximal coverage - i.e. remove all peaks that DON'T have a cell with at least \code{abs_cov_thresh} UMIs
#' @param min_peak_length (optional) remove all peaks that are less than \code{min_cov_thresh} base-pairs long
#' @param max_peak_length (optional) remove all peaks that are more than \code{min_cov_thresh} base-pairs long
#'
#' @return scatac - the filtered ScATAC object
#' @examples
#' \dontrun{
#' my_scatac <- filter_features(my_scatac, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
#' table(my_intervals$peak_annot)
#' }
#' @export
filter_features <- function(scatac, minimal_max_umi = NULL, min_peak_length = NULL, max_peak_length = NULL) {
    if (!is.null(minimal_max_umi)) {
        rmx <- sparseMatrixStats::rowMaxs(scatac@mat)
        low_max_peaks <- scatac@peaks$peak_name[rmx < minimal_max_umi]
    } else {low_max_peaks <- c()}
    if (!is.null(min_peak_length) || !is.null(max_peak_length)) {
        peak_len <- scatac@peaks$end - scatac@peaks$start
        if (!is.null(min_peak_length)) {
            too_short_peaks <- scatac@peaks$peak_name[scatac@peaks$len < min_peak_length]
        } else {too_short_peaks <- c()}
        if (!is.null(max_peak_length)) {
            too_long_peaks <- scatac@peaks$peak_name[scatac@peaks$len > max_peak_length]
        } else {too_long_peaks <- c()}
    }
    peaks_to_remove <- c(low_max_peaks, too_short_peaks, too_long_peaks)
    scatac <- ignore_peaks(scatac, peaks_to_remove)
    return(scatac)
}
    




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
    mcatac@mat <- mcatac@mat + quantile(mns, mean_thresh_q)
    mns <- rowMeans(mcatac@mat)
    sds <- MatrixGenerics::rowSds(mcatac@mat)
    covs <- sds/mns
    lmi <- lm(cov ~ mn, data <- as.data.frame(list('mn' <- log10(mns), 'cov' <- log10(covs))))
    pred_df <- as.data.frame(list('pred' <- lmi$coefficients[[1]] + lmi$coefficients[[2]]*log10(mns), 'cov' <- log10(covs)))
    pred_df <- dplyr::mutate(pred_df, diff = cov - pred)
    if (method == 'bmq') {
        lpm <- log10(mns)
        qcut <- cut(lpm, breaks <- seq(min(lpm), max(lpm), l=num_bins))
        feat_select <- unlist(sapply(sort(unique(qcut[which(lpm >= log10(mn_thresh))])),
                             function(nm) rownames(mcatac@mat)[which(covs >= quantile(covs[qcut == nm], cov_q_thresh, na.rm=T) & qcut == nm)]))
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
        clusters_selected <- as.numeric(names(cluster_mean_diff)[cluster_mean_diff > 0])
        plot(X[,1], X[,2], main = glue::glue('Selected clusters - {clusters_selected}'), cex = 0.05)
        clrs = c('red', 'blue', 'green', 'cyan', 'purple', 'orange')
        purrr::walk(seq_along(clusters_selected), function(cl, i) {
            points(X[fc == cl[[i]],1], X[fc == cl[[i]],2], col = clrs[[i]], cex= 0.05)
        }, cl = clusters_selected)
        return(mcatac@peaks[mcatac@peaks$peak_name %in% rownames(pred_df)[pred_df$fc %in% clusters_selected],])
    }
    else {cli_abort('Method should be NULL, "bmq" or "gmm"')}
}

#' Find overlaps with ENCODE blacklists
#'
#' @description Identify peaks in the data which overlap (or are adjacent to?) regions blacklisted by ENCODE as having universally high DNAse HS or ChIP signal (basically mapping artifacts)
#' See https://doi.org/10.1038/s41598-019-45839-z for more details
#' @param scatac (optional) an ScATAC object
#' @param mcatac (optional) an McATAC object
#' @param peaks (optional) the intervals set to check
#' @param genome (optional, required if checking peaks directly) the interval set to check
#'
#' @return blacklist_overlaps - a PeakIntervals object with peaks identified as overlapping blacklisted regions
#' @examples
#' \dontrun{
#' blacklist_overlaps <- find_blacklist_overlaps(my_intervals, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
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
    peaks$intervalID <- 1:nrow(peaks)
    misha.ext::gset_genome(genome)
    blacklist_name <- glue::glue('ENCODE_blacklist_{genome}')
    if (gintervals.exists(blacklist_name)) {
        blacklist <- gintervals.load(blacklist_name)
    }
    nei_blk_pks <- gintervals.neighbors(blacklist, peaks, maxdist = 0, mindist = 0, maxneighbors = 10)
    blacklist_overlaps <- peaks[unique(nei_blk_pks$intervalID),]
    return(blacklist_overlaps)
}