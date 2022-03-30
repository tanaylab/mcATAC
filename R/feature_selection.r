#' Filter features by summary statistics
#'
#' @description Remove ATAC peaks with low coverage, or that are too long etc.
#' @param scatac the ScATAC object to analyze
#' @param minimal_max_umi (optional) threshold on minimum maximal coverage - i.e. remove all peaks that DON'T have a cell with at least \code{abs_cov_thresh} UMIs
#' @param min_peak_length (optional) remove all peaks that are less than \code{min_cov_thresh} base-pairs long
#' @param max_peak_length (optional) remove all peaks that are more than \code{min_cov_thresh} base-pairs long
#'
#' @return scatac - the filtered ScATAC object
#' @examples
#' \dontrun{
#' ## quantiles of peak lengths before filtering
#' quantile(my_scatac@peaks$end - my_scatac@peaks$start, (0:10) / 10)
#' my_scatac <- filter_features(scatac = my_scatac, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
#'
#' ## quantiles of peak lengths after filtering
#' quantile(my_scatac@peaks$end - my_scatac@peaks$start, (0:10) / 10)
#' }
#' @export
filter_features <- function(scatac, minimal_max_umi = NULL, min_peak_length = NULL, max_peak_length = NULL) {
    if (!is.null(minimal_max_umi)) {
        rmx <- sparseMatrixStats::rowMaxs(scatac@mat)
        low_max_peaks <- scatac@peaks$peak_name[rmx < minimal_max_umi]
    } else {
        low_max_peaks <- c()
    }
    peak_len <- scatac@peaks$end - scatac@peaks$start
    if (!is.null(min_peak_length)) {
        too_short_peaks <- scatac@peaks$peak_name[peak_len < min_peak_length]
    } else {
        too_short_peaks <- c()
    }
    if (!is.null(max_peak_length)) {
        too_long_peaks <- scatac@peaks$peak_name[peak_len > max_peak_length]
    } else {
        too_long_peaks <- c()
    }
    peaks_to_remove <- unique(c(low_max_peaks, too_short_peaks, too_long_peaks))
    if (length(peaks_to_remove) > 0) {
        if (length(low_max_peaks) > 0) {
            cli::cli_li("{.val {length(low_max_peaks)}} features had a maximal UMI count less than {.field {minimal_max_umi}}")
        }
        if (length(too_short_peaks) > 0) {
            cli::cli_li("{.val {length(too_short_peaks)}} features were shorter than {.field {min_peak_length}bp}")
        }
        if (length(too_long_peaks) > 0) {
            cli::cli_li("{.val {length(too_long_peaks)}} features were longer than {.field {max_peak_length}bp}")
        }
        scatac <- atac_ignore_peaks(scatac, peaks_to_remove)
    } else {
        cli_alert_warning("No peaks that violate the criteria were found. Returning original object")
    }
    return(scatac)
}


#' Plot the distribution of the peak length
#'
#' @param atac an ScATAC or McATAC object
#'
#' @return a ggplot object with the distribution of peak length
#'
#' @examples
#' \dontrun{
#' plot_peak_length_distribution(scatac)
#' }
#'
#' @export
plot_peak_length_distribution <- function(atac) {
    gg <- atac@peaks %>%
        mutate(len = end - start) %>%
        ggplot(aes(x = len)) +
        stat_density() +
        scale_x_log10() +
        labs(
            x = "Peak length (bp)",
            y = "Density",
            title = "Peak length distribution",
            subtitle = glue("n = {nrow(atac@peaks)}"),
            caption = glue("object id: {atac@id}")
        )

    return(gg)
}

#' Plot the coverage distribution of the peaks
#'
#' @param atac an ScATAC or McATAC object
#'
#' @return a ggplot object with the distribution of peak coverage
#'
#' @examples
#' \dontrun{
#' plot_peak_coverage_distribution(scatac)
#' }
#' @export
plot_peak_coverage_distribution <- function(atac) {
    gg <- atac@peaks %>%
        mutate(cov = sparseMatrixStats::rowSums2(atac@mat)) %>%
        ggplot(aes(x = cov)) +
        stat_density() +
        scale_x_log10() +
        labs(
            x = "# of UMIs",
            y = "Density",
            title = glue("Peak coverage distribution"),
            subtitle = glue("n = {nrow(atac@peaks)}"),
            caption = glue("object id: {atac@id}")
        )
    return(gg)
}


#' Plot the maximal per-cell coverage distribution of the peaks
#'
#' @param atac an ScATAC or McATAC object
#'
#' @return a ggplot object with the distribution of maximal per-cell coverage for each peak
#'
#' @examples
#' \dontrun{
#' plot_peak_max_cov_distribution(scatac)
#' }
#'
#' @export
plot_peak_max_cov_distribution <- function(atac) {
    gg <- atac@peaks %>%
        mutate(cov = sparseMatrixStats::rowMaxs(atac@mat)) %>%
        count(cov) %>%
        mutate(
            cov = ifelse(cov > 10, ">10", as.character(cov)),
            cov = factor(cov, levels = c(0:10, ">10"))
        ) %>%
        ggplot(aes(x = cov, y = n)) +
        geom_col() +
        labs(
            x = "Maximal per-cell coverage",
            y = "Density",
            title = glue("Maximal per-cell coverage distribution"),
            subtitle = glue("n = {nrow(atac@peaks)}"),
            caption = glue("object id: {atac@id}")
        )
    return(gg)
}

#' Find dynamic peaks in McATAC matrix
#'
#' @description This function identifies "dynamic" peaks, i.e. those that have high expression only in a subset of the cells. They are identified by overdispersion in the coefficient of variation (std.dev./mean) per quantiles.
#' @param mcatac the McATAC object to analyze
#' @param method (optional) either 'bmq' (default) or 'gmm'; 'bmq' (binned-mean quantiles) bins the log-mean of all peaks (averaged across metacells) and
#'   selects all peaks with a coefficient of variation above some quantile in each bin. More controlled
#'  'gmm' fits a Gaussian mixture model to the log10(COV) vs. log10(mean) distribution,
#'   and selects peaks in clusters that show overdispersion in the COV.
#' @param plot plot the peak mean vs coefficient of variation (both in log10 scale). Note that it is highly recommended to look at the
#' scatter plot before proceeding, so set this parameter to FALSE only __after__ you made sure that the scatter looks reasonable.
#' @param mean_thresh_q (optional) threshold quantile on peaks' mean
#' @param cov_q_thresh (optional) threshold on minimum COV quantile to consider as dynamic in each bin
#' @param num_bins (optional) number of bins to divide features' means into
#' @param gmm_g (optional) number of groups for 'gmm'
#'
#' @return a PeakIntervals object with peaks identified as dynamic. If \code{plot = TRUE} the selected points would plotted.
#' @examples
#' \dontrun{
#' dynamic_peaks_by_bmq <- identify_dynamic_peaks(mcatac, method = "bmq", mean_thresh_q = 0.1, cov_q_thresh = 0.6, num_bins = 100)
#' dynamic_peaks_by_gmm <- identify_dynamic_peaks(mcatac, method = "gmm", gmm_g = 3)
#' }
#' @export
identify_dynamic_peaks <- function(mcatac, method = "bmq", plot = TRUE, mean_thresh_q = 0.1, cov_q_thresh = 0.75, num_bins = 200, gmm_g = 4) {
    assert_atac_object(mcatac)
    mns <- Matrix::rowMeans(mcatac@mat, na.rm = TRUE)
    sds <- sparseMatrixStats::rowSds(mcatac@mat, na.rm = TRUE)
    covs <- sds / mns
    lpm <- log10(mns)
    data <- tibble(mn = log10(mns), cov = log10(covs))
    lmi <- lm(cov ~ mn, data)
    pred_df <- tibble(
        pred = lmi$coefficients[[1]] + lmi$coefficients[[2]] * log10(mns),
        cov = log10(covs),
        peak_name = mcatac@peaks$peak_name
    ) %>%
        mutate(diff = cov - pred)
    if (method == "bmq") {
        qcut <- cut(lpm, breaks <- seq(min(lpm), max(lpm), l = num_bins))
        mn_thresh <- quantile(mns, mean_thresh_q)
        feat_select <- unlist(sapply(
            sort(unique(qcut[which(lpm >= log10(mn_thresh))])),
            function(nm) rownames(mcatac@mat)[which(covs >= quantile(covs[qcut == nm], cov_q_thresh, na.rm = T) & qcut == nm)]
        ))
        pred_df$is_bmq <- pred_df$peak_name %in% feat_select
        if (plot) {
            cli_li("Plotting log10(mean) vs log10(sd/mean)")
            smoothScatter(lpm, log10(covs), cex = 0.05, xlab = "log10(mean)", ylab = "log10(sd/mean)")
            points(lpm[feat_select], log10(covs[feat_select]), cex = 0.05, col = "red")
        }
    } else if (method == "gmm") {
        X <- cbind(lpm, log10(covs))
        fit_gmm <- mclust::Mclust(X, G = gmm_g, modelNames = "VVV")
        pred_df$fc <- fit_gmm$classification
        cluster_mean_diff <- tapply(pred_df$diff, pred_df$fc, mean)
        clusters_selected <- as.numeric(names(cluster_mean_diff)[cluster_mean_diff > 0])
        if (plot) {
            cli_li("Plotting log10(mean) vs log10(sd/mean)")
            plot(X[, 1], X[, 2], main = glue::glue("Selected clusters: {paste(clusters_selected, collapse = ',')}"), cex = 0.05, xlab = "log10(mean)", ylab = "log10(sd/mean)")
            clrs <- c("red", "blue", "green", "cyan", "purple", "orange")
            purrr::walk(seq_along(clusters_selected), function(cl, i) {
                points(X[pred_df$fc == cl[[i]] & pred_df$diff > 0, 1], X[pred_df$fc == cl[[i]] & pred_df$diff > 0, 2], col = clrs[[i]], cex = 0.05)
            }, cl = clusters_selected)
        }
        feat_select <- pred_df %>%
            filter(fc %in% clusters_selected, diff > 0) %>%
            pull(peak_name)
    } else {
        cli_abort("Method should be {.code NULL}, {.code bmq} or {.code gmm}")
    }
    res <- mcatac@peaks %>% filter(peak_name %in% feat_select)
    cli_alert_success("Identified {.val {nrow(res)}} dynamic peaks (out of {.val {nrow(mcatac@peaks)}}) using the {.field '{method}'} method.")
    return(res)
}

#' Find overlaps with ENCODE blacklists
#'
#' @description Identify peaks in the data which overlap (or are adjacent to?) regions blacklisted by ENCODE as having universally high DNAse HS or ChIP signal (basically mapping artifacts)
#' See https://doi.org/10.1038/s41598-019-45839-z for more details
#' @param atac (optional) an ScATAC or McATAC object
#' @param peaks (optional) the intervals set to check
#' @param genome (optional, required if checking peaks directly) the genome of the peaks
#' @param max_dist_to_blacklist_region (optional) distance to nearest blacklist region which still qualifies for being blacklisted
#' @param blacklist_name name of the blacklist intervals to use (default: "ENCODE.blacklist")
#'
#' @return blacklist_overlaps - a PeakIntervals object with peaks identified as overlapping blacklisted regions
#' @examples
#' \dontrun{
#' blacklist_overlaps <- find_blacklist_overlaps(atac = my_mcatac, max_dist_to_blacklist_region = 100)
#' blacklist_overlaps <- find_blacklist_overlaps(peaks = my_peak_set, genome = "mm10")
#' }
#' @export
find_blacklist_overlaps <- function(atac = NULL, peaks = NULL, genome = NULL, max_dist_to_blacklist_region = 0, blacklist_name = "ENCODE.blacklist") {
    assert_atac_object(atac)
    if (!is.null(atac) && is.null(peaks)) {
        peaks <- atac@peaks
        genome <- atac@genome
    } else if (!is.null(atac) && !is.null(peaks)) {
        cli_alert_warning("Both ATAC object and peaks specified, which is redundant. Ignoring peaks.")
        peaks <- atac@peaks
        genome <- atac@genome
    } else if (is.null(peaks)) {
        cli_abort("Must specify either scatac, mcatac or peaks")
    } else if (!is.null(peaks) && is.null(genome)) {
        cli_abort("Must specify genome if analyzing peaks directly")
    }

    misha.ext::gset_genome(genome)
    if (!gintervals.exists(blacklist_name)) {
        cli_abort("Blacklist intervals {.field {blacklist_name}} does not exist")
    }
    nei_blk_pks <- misha.ext::gintervals.filter(as.data.frame(peaks), blacklist_name, max_distance = max_dist_to_blacklist_region)

    blacklist_overlaps <- peaks %>% filter(peak_name %in% nei_blk_pks$peak_name)

    return(blacklist_overlaps)
}
