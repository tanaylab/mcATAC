#' Calculate coverage statistics of the peaks
#'
#' @description For each peak, calculate:
#' \enumerate{
#' \item{cov: }{The total number of UMIs in the peak}
#' \item{max_cov: }{The maximal per-cell coverage of the peak}
#' \item{len: }{The length of the peak}
#' \item{cov_density: }{The number of UMIs the number of UMIs per \code{scale} bp.}
#' }
#'
#' @param atac an ScATACPeaks or McATACPeaks object
#' @param scale for the 'cov_density' field, calculate the number of UMIs per \code{scale}
#' bp (default: 100)
#'
#' @return a data frame with peaks, and their length ('len'), their total coverage ('cov'),
#' maximal per-cell coverage ('max_cov') and their coverage density ('cov_density')
#'
#' @examples
#' \dontrun{
#' get_peak_coverage_stats(atac)
#' }
#'
#' @export
get_peak_coverage_stats <- function(atac, scale = 100) {
    assert_atac_object(atac)
    df <- atac@peaks %>%
        mutate(
            len = end - start,
            cov = sparseMatrixStats::rowSums2(atac@mat),
            max_cov = sparseMatrixStats::rowMaxs(atac@mat),
            cov_density = cov / len * scale
        )
    return(df)
}

#' Filter features by summary statistics
#'
#' @description Remove ATAC peaks with low coverage, or that are too long etc.
#' @param atac_sc the ScATACPeaks object to analyze
#' @param minimal_max_umi (optional) threshold on minimum maximal coverage - i.e. remove all peaks that DON'T have a cell with at least \code{abs_cov_thresh} UMIs
#' @param min_peak_length (optional) remove all peaks that are less than \code{min_cov_thresh} base-pairs long
#' @param max_peak_length (optional) remove all peaks that are more than \code{min_cov_thresh} base-pairs long
#' @param max_peak_density (optional) remove all peaks that have more than \code{max_peak_density} UMIs per 100bp
#'
#' @return atac_sc - the filtered ScATACPeaks object
#' @examples
#' \dontrun{
#' ## quantiles of peak lengths before filtering
#' quantile(my_atac_sc@peaks$end - my_atac_sc@peaks$start, (0:10) / 10)
#' my_atac_sc <- filter_features(atac_sc = my_atac_sc, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
#'
#' ## quantiles of peak lengths after filtering
#' quantile(my_atac_sc@peaks$end - my_atac_sc@peaks$start, (0:10) / 10)
#' }
#' @export
filter_features <- function(atac_sc, minimal_max_umi = NULL, min_peak_length = NULL, max_peak_length = NULL, max_peak_density = NULL) {
    assert_atac_object(atac_sc)
    peak_stats <- get_peak_coverage_stats(atac_sc, scale = 100)

    if (!is.null(min_peak_length)) {
        too_short_peaks <- peak_stats$peak_name[peak_stats$len < min_peak_length]
        if (length(too_short_peaks) > 0) {
            peak_stats <- peak_stats %>%
                filter(peak_name %!in% too_short_peaks)
            cli::cli_li("{.val {length(too_short_peaks)}} features were shorter than {.field {min_peak_length}bp}")
        }
    }

    if (!is.null(max_peak_length)) {
        too_long_peaks <- peak_stats$peak_name[peak_stats$len > max_peak_length]
        if (length(too_long_peaks) > 0) {
            peak_stats <- peak_stats %>%
                filter(peak_name %!in% too_long_peaks)
            cli::cli_li("{.val {length(too_long_peaks)}} features were longer than {.field {max_peak_length}bp}")
        }
    }

    if (!is.null(minimal_max_umi)) {
        low_max_peaks <- peak_stats$peak_name[peak_stats$max_cov < minimal_max_umi]
        if (length(low_max_peaks) > 0) {
            peak_stats <- peak_stats %>%
                filter(peak_name %!in% low_max_peaks)
            cli::cli_li("{.val {length(low_max_peaks)}} features had a maximal UMI count less than {.field {minimal_max_umi}}")
        }
    }

    if (!is.null(max_peak_density)) {
        too_dens_peaks <- peak_stats$peak_name[peak_stats$cov_density > max_peak_density]
        if (length(too_dens_peaks) > 0) {
            peak_stats <- peak_stats %>%
                filter(peak_name %!in% too_dens_peaks)
            cli::cli_li("{.val {length(too_dens_peaks)}} features had a peak density of more than {.field {max_peak_density}} UMIs per 100bp")
        }
    }

    peaks_to_remove <- setdiff(atac_sc@peaks$peak_name, peak_stats$peak_name)
    if (length(peaks_to_remove) > 0) {
        atac_sc <- atac_ignore_peaks(atac_sc, peaks_to_remove)
    } else {
        cli_alert_warning("No peaks that violate the criteria were found. Returning original object")
    }
    return(atac_sc)
}


#' Plot the distribution of the peak length
#'
#' @param atac an ScATACPeaks or McATACPeaks object
#'
#' @return a ggplot object with the distribution of peak length
#'
#' @examples
#' \dontrun{
#' plot_peak_length_distribution(atac_sc)
#' }
#'
#' @export
plot_peak_length_distribution <- function(atac) {
    assert_atac_object(atac)
    gg <- get_peak_coverage_stats(atac) %>%
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
#' @param atac an ScATACPeaks or McATACPeaks object
#'
#' @return a ggplot object with the distribution of peak coverage
#'
#' @examples
#' \dontrun{
#' plot_peak_coverage_distribution(atac_sc)
#' }
#' @export
plot_peak_coverage_distribution <- function(atac) {
    assert_atac_object(atac)
    gg <- get_peak_coverage_stats(atac) %>%
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
#' @param atac an ScATACPeaks or McATACPeaks object
#'
#' @return a ggplot object with the distribution of maximal per-cell coverage for each peak
#'
#' @examples
#' \dontrun{
#' plot_peak_max_cov_distribution(atac_sc)
#' }
#'
#' @export
plot_peak_max_cov_distribution <- function(atac) {
    assert_atac_object(atac)
    gg <- get_peak_coverage_stats(atac) %>%
        count(max_cov) %>%
        mutate(
            max_cov = ifelse(max_cov > 10, ">10", as.character(max_cov)),
            max_cov = factor(max_cov, levels = c(as.character(0:10), ">10")),
            p = n / sum(n)
        ) %>%
        ggplot(aes(x = max_cov, y = p)) +
        geom_col() +
        scale_y_continuous(labels = scales::percent) +
        labs(
            x = "Maximal per-cell coverage",
            y = "% of peaks",
            title = glue("Maximal per-cell coverage distribution"),
            subtitle = glue("n = {nrow(atac@peaks)}"),
            caption = glue("object id: {atac@id}")
        )
    return(gg)
}

#' Plot the coverage density of the peaks
#'
#' @description For each peak, plot a scatter of the peak length vs the number
#' of UMIs per \code{scale} bp.
#'
#' @inheritParams get_peak_coverage_stats
#' @inheritParams scattermore::geom_scattermore
#' @inheritDotParams scattermore::geom_scattermore
#'
#' @examples
#' \dontrun{
#' plot_peak_coverage_density(atac)
#' }
#'
#' @export
plot_peak_coverage_density <- function(atac, scale = 100, pointsize = 1.5, ...) {
    assert_atac_object(atac)
    gg <- get_peak_coverage_stats(atac, scale = scale) %>%
        ggplot(aes(x = len, y = cov_density)) +
        scattermore::geom_scattermore(pointsize = 1.5, ...) +
        scale_x_log10() +
        labs(
            x = "Read length (bp)",
            y = glue("Reads per {scale} bp"),
            title = glue("Peak coverage density"),
            subtitle = glue("n = {nrow(atac@peaks)}"),
            caption = glue("object id: {atac@id}")
        )
    return(gg)
}

#' Find dynamic peaks in McATACPeaks matrix
#'
#' @description This function identifies "dynamic" peaks, i.e. those that have high expression only in a subset of the cells. They are identified by overdispersion in the coefficient of variation (std.dev./mean) per quantiles.
#' @param atac_mc the McATACPeaks object to analyze
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
#' dynamic_peaks_by_bmq <- identify_dynamic_peaks(atac_mc, method = "bmq", mean_thresh_q = 0.1, cov_q_thresh = 0.6, num_bins = 100)
#' dynamic_peaks_by_gmm <- identify_dynamic_peaks(atac_mc, method = "gmm", gmm_g = 3)
#' }
#' @export
identify_dynamic_peaks <- function(atac_mc, method = "bmq", plot = TRUE, mean_thresh_q = 0.1, cov_q_thresh = 0.75, num_bins = 200, gmm_g = 4) {
    assert_atac_object(atac_mc)
    mns <- Matrix::rowMeans(atac_mc@mat, na.rm = TRUE)
    sds <- sparseMatrixStats::rowSds(atac_mc@mat, na.rm = TRUE)
    covs <- sds / mns
    lpm <- log10(mns)
    data <- tibble(mn = log10(mns), cov = log10(covs))
    lmi <- lm(cov ~ mn, data)
    pred_df <- tibble(
        pred = lmi$coefficients[[1]] + lmi$coefficients[[2]] * log10(mns),
        cov = log10(covs),
        peak_name = atac_mc@peaks$peak_name
    ) %>%
        mutate(diff = cov - pred)
    if (method == "bmq") {
        qcut <- cut(lpm, breaks <- seq(min(lpm), max(lpm), l = num_bins))
        mn_thresh <- quantile(mns, mean_thresh_q)
        feat_select <- unlist(sapply(
            sort(unique(qcut[which(lpm >= log10(mn_thresh))])),
            function(nm) rownames(atac_mc@mat)[which(covs >= quantile(covs[qcut == nm], cov_q_thresh, na.rm = T) & qcut == nm)]
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
    res <- atac_mc@peaks %>% filter(peak_name %in% feat_select)
    cli_alert_success("Identified {.val {nrow(res)}} dynamic peaks (out of {.val {nrow(atac_mc@peaks)}}) using the {.field '{method}'} method.")
    return(res)
}

#' Find overlaps with ENCODE blacklists
#'
#' @description Identify peaks in the data which overlap (or are adjacent to?) regions blacklisted by ENCODE as having universally high DNAse HS or ChIP signal (basically mapping artifacts)
#' See https://doi.org/10.1038/s41598-019-45839-z for more details
#' @param atac (optional) an ScATACPeaks or McATACPeaks object
#' @param peaks (optional) the intervals set to check
#' @param genome (optional, required if checking peaks directly) the genome of the peaks
#' @param max_dist_to_blacklist_region (optional) distance to nearest blacklist region which still qualifies for being blacklisted
#' @param blacklist_name name of the blacklist intervals to use (default: "ENCODE.blacklist")
#'
#' @return blacklist_overlaps - a PeakIntervals object with peaks identified as overlapping blacklisted regions
#' @examples
#' \dontrun{
#' blacklist_overlaps <- find_blacklist_overlaps(atac = my_atac_mc, max_dist_to_blacklist_region = 100)
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
        cli_abort("Must specify either ScATACPeaks, McATACPeaks or peaks")
    } else if (!is.null(peaks) && is.null(genome)) {
        cli_abort("Must specify genome if analyzing peaks directly")
    }

    misha.ext::gset_genome(genome)
    if (!gintervals.exists(blacklist_name)) {
        cli_abort("Blacklist intervals {.field {blacklist_name}} do not exist")
    }
    nei_blk_pks <- misha.ext::gintervals.filter(as.data.frame(peaks), blacklist_name, max_distance = max_dist_to_blacklist_region)

    blacklist_overlaps <- peaks %>% filter(peak_name %in% nei_blk_pks$peak_name)

    return(blacklist_overlaps)
}
