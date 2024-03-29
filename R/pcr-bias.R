#' Normalize Marginal Track
#'
#' This function normalizes a marginal track by dividing it by a smoothed track.
#' The normalization is performed by taking the geometric mean of window_size windows
#' in the original track and dividing each value by 2 raised to the power of the floored
#' smoothed track value minus 1. The floored smoothed track value is calculated as the
#' floor_percentile of the smoothed track.
#'
#' @param orig_track The original track to be normalized.
#' @param normed_track The name of the normalized track to be created.
#' @param smoothed_track The name of the smoothed track to be created.
#' The original track is normalized by this track. If NULL, a temporary track is created.
#' @param window_size The size of the windows used for smoothing the track.
#' @param iterator The iterator used for smoothing the track.
#' @param floor_percentile The percentile used to calculate the floored smoothed track value.
#' @param overwrite Logical value indicating whether to overwrite an existing normalized track.
#' @return None
#'
#' @examples
#' \dontrun{
#' # mcc_to_marginal runs internally normalize_marginal when normalize = TRUE
#' mcc_to_marginal_track(mc_counts, "pbmc_mc.marginal", normalize = TRUE)
#'
#' # alternatively, you can run normalize_marginal directly
#' # normalize_marginal("pbmc_mc.marginal", "pbmc_mc.marginal_normed", "pbmc_mc.marginal_smoothed")
#' }
#'
#' @export
normalize_marginal <- function(orig_track, normed_track, smoothed_track = NULL, window_size = 20e3, iterator = 1e3, floor_percentile = 0.1, overwrite = FALSE) {
    if (is.null(smoothed_track)) {
        gdir.create("temp", showWarnings = FALSE)
        smoothed_track <- temp_track_name("temp.")
    }

    if (gtrack.exists(normed_track)) {
        if (overwrite) {
            cli_alert_warning("Removing previous track {.val {normed_track}}")
            gtrack.rm(normed_track, force = TRUE)
        } else {
            cli_abort("{normed_track} already exists. Use 'overwrite = TRUE' to overwrite.")
        }
    }

    if (gtrack.exists(smoothed_track)) {
        if (overwrite) {
            cli_alert_warning("Removing previous track {.val {smoothed_track}}")
            gtrack.rm(smoothed_track, force = TRUE)
        } else {
            cli_abort("{smoothed_track} already exists. Use 'overwrite = TRUE' to overwrite.")
        }
    }

    cli_alert_info("Smoothing track {.val {orig_track}} with window size {.val {window_size}} and iterator {.val {iterator}}")
    # Create a track of smoothed ATAC signal (geometric mean of window_size windows)
    gtrack.smooth(
        track = smoothed_track,
        description = glue("Smoothed track of {orig_track} with window size {window_size} and iterator {iterator}"),
        expr = glue("log2(1 + ifelse(is.na({orig_track}), 0, {orig_track}))"),
        winsize = window_size,
        weight_thr = 0,
        smooth_nans = FALSE,
        alg = "MEAN",
        iterator = iterator
    )

    # we set the minimal value to be the floor_percentile of the smoothed track (and not less than 1)
    floor_val <- gquantiles(smoothed_track, floor_percentile)

    # We now normalize the marginal track by the smoothed track.
    # The norm is in log-scale of the value + 1 (geometric mean) hence the exponent of 2 and the -1
    cli_alert_info("Normalizing track {.val {orig_track}} to track {.val {normed_track}}")
    expr <- glue("{orig_track} / 2^({smoo_track_floored} - 1)", smoo_track_floored = glue("pmax(1, 1 + ifelse(is.na({smoothed_track}), 0, {smoothed_track} - {floor_val}))"))
    orig_description <- gtrack.attr.export(orig_track)$description
    description <- glue("{orig_description} (normalized by a geometric mean of {window_size}bp windows and iterator of {iterator})")
    gtrack.create(
        track = normed_track,
        description = description,
        expr = expr,
        iterator = gtrack.info(orig_track)$bin.size
    )
}

#' Normalize McPeaks data using marginal coverage
#'
#' The function calculates the (punctured) marginal coverage in \code{window_size} windows around each peak,
#' and then normalizes the data by dividing the original data by the punctured marginal coverage, while setting
#' the minimum to \code{minimal_quantile} of the punctured marginal coverage.
#'
#' @param mcatac The McPeaks object containing the data to be normalized.
#' @param marginal_track The marginal track
#' @param widnow_size The size of the windows to use around each peak for normalization.
#' @param epsilon The epsilon value added to the egc matrix (default is 1e-5).
#'
#' @return The McPeaks object with normalized egc values.
#'
#'
#' @examples
#' \dontrun{
#' mcatac <- import_from_matrix(mat, peaks, genome = "hg38", class = "McPeaks")
#' mcatac <- normalize_egc(mcatac, "pbmc_mc.marginal")
#' }
#'
#' @export
normalize_egc <- function(mcatac, marginal_track, window_size = 1e4, epsilon = 1e-5, minimal_quantile = 0.1) {
    gvtrack.create("marginal", marginal_track, func = "sum")
    gvtrack.create("marginal_20k", marginal_track, func = "sum")
    gvtrack.iterator("marginal_20k", sshift = -window_size / 2, eshift = window_size / 2)
    peaks_metadata <- misha.ext::gextract.left_join(
        c("marginal", "marginal_20k"),
        intervals = mcatac@peaks,
        iterator = mcatac@peaks
    ) %>%
        mutate(marginal_20k_punc = marginal_20k - marginal) %>%
        mutate(norm_f = marginal_20k_punc / quantile(marginal_20k_punc, minimal_quantile)) %>%
        mutate(norm_f = ifelse(norm_f < 1, 1, norm_f)) %>%
        mutate(norm_f = 1 / norm_f)
    mc_mat_s <- mcatac@mat * peaks_metadata$norm_f
    mc_egc_s <- t(t(mc_mat_s) / colSums(mc_mat_s))
    
    mcatac@egc <- mc_egc_s
    mcatac@fp <- calc_mc_fp(mcatac)

    return(mcatac)
}

