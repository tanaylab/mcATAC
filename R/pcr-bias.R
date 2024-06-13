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


#' Normalize McPeaks data using constitutive peaks
#'
#' The function normalizes the McPeaks data using the values of the mean of constitutive peaks in a specified cell type. The normalized data is then scaled to a specified quantile and capped at 1.
#'
#' @param mcatac The McPeaks object containing the data to be normalized. The object should contain a field named {.field const} in {.field peaks}. Use \code{add_const_peaks} to add this field.
#' @param norm_type The cell type to use as the reference for normalization
#' @param norm_quant The quantile to which the normalized data should be scaled (default is 1).
#'
#' @return The McPeaks object with normalized egc values.
#'
#' @examples
#' \dontrun{
#' mcatac <- import_from_matrix(mat, peaks, genome = "mm10", class = "McPeaks")
#' mcatac <- normalize_egc(mcatac, "gastrulation.marginal")
#' mcatac <- add_const_peaks(mcatac, -16)
#' mcatac <- normalize_const(mcatac, norm_type = "Epiblast")
#' }
#'
#' @export
normalize_const <- function(mcatac, norm_type, norm_quant = 1, scaling_quant = 0.9) {
    if (!has_name(mcatac@peaks, "const")) {
        cli_abort("No field named {.field const} in peaks. Please run {.code add_const_peaks}")
    }

    const_peaks <- mcatac@peaks$peak_name[mcatac@peaks$const]
    egc <- mcatac@egc

    if (!has_cell_type(mcatac)) {
        cli_abort("No cell_type in metadata")
    }

    if (!(norm_type %in% mcatac@metadata$cell_type)) {
        cli_abort("norm_type {.val {norm_type}} is not in {.field metadata$cell_type}")
    }

    norm_type_mcs <- mcatac@metadata$metacell[mcatac@metadata$cell_type == norm_type]

    const_peak_score_all <- colSums(egc[const_peaks, ])

    const_peak_score_norm_type <- mean(const_peak_score_all[norm_type_mcs])
    const_peak_score_norm <- const_peak_score_all / const_peak_score_norm_type

    egc_norm <- t(t(egc) / const_peak_score_norm)
    egc_norm <- egc_norm / quantile(as.vector(egc_norm), norm_quant)
    egc_norm[egc_norm > 1] <- 1

    egc_norm <- egc_norm * max(egc[, norm_type_mcs]) / scaling_quant

    mcatac@egc <- egc_norm

    return(mcatac)
}




#' Calculate constitutive peaks
#'
#' This function calculates the constitutive peaks based on the minimal expression of each peak across different cell types.
#'
#' @param mcatac McPeaks object.
#' @param const_threshold The threshold value for determining constitutive peaks. A peak would be considered constitutive if its minimal expression across all cell types is greater than or equal to this threshold.
#'
#' @return A named logical vector indicating whether each peak is constitutive or not.
#'
#' @examples
#' \dontrun{
#' mcatac <- import_from_matrix(mat, peaks, genome = "mm10", class = "McPeaks")
#' mcatac <- normalize_egc(mcatac, "gastrulation.marginal")
#' const_v <- calc_const_peaks(mcatac, -16)
#' }
#'
#' @export
calc_const_peaks <- function(mcatac, const_threshold) {
    mc_egc <- mcatac@egc
    mc_egc <- as.matrix(log2(mc_egc + 1e-5))

    if (!has_cell_type(mcatac)) {
        cli_abort("No cell_type in metadata")
    }

    # get cell type egc average for each peak
    ct_egc <- t(tgs_matrix_tapply(mc_egc, mcatac@metadata$cell_type, mean, na.rm = TRUE))

    const_v <- tibble(cre = rownames(ct_egc), const = apply(ct_egc, 1, min, na.rm = TRUE) >= const_threshold) %>%
        select(cre, const) %>%
        deframe()

    return(const_v)
}

#' Add constitutive peaks to McPeaks object
#'
#' This function adds a field named {.field const} to the {.field peaks} field of the McPeaks object. The field contains a named logical vector indicating whether each peak is constitutive or not.
#'
#' @return The McPeaks object with the added field.
#'
#' @examples
#' \dontrun{
#' mcatac <- import_from_matrix(mat, peaks, genome = "mm10", class = "McPeaks")
#' mcatac <- normalize_egc(mcatac, "gastrulation.marginal") #'
#' mcatac <- add_const_peaks(mcatac, -16)
#' }
#'
#' @inheritParams calc_const_peaks
#' @export
add_const_peaks <- function(mcatac, const_threshold) {
    mcatac@peaks$const <- calc_const_peaks(mcatac, const_threshold)
    return(mcatac)
}

#' Normalize the mcATAC object to probability values
#'
#' This function normalizes the mcATAC object to probability values by setting values above a threshold to the threshold value,
#' and then scaling the values to the range [0, 1].
#'
#' @param mcatac The mcPeaks object to be normalized.
#' @param prob1_thresh The threshold value for setting values above it to the threshold value. If not provided, the function will calculate the threshold using the "const" field in the peaks.
#' @param const_quantile The quantile value used to calculate the threshold when "prob1_thresh" is not provided.
#'
#' @return The normalized mcATAC object.
#'
#' @examples
#' \dontrun{
#' mcatac <- import_from_matrix(mat, peaks, genome = "mm10", class = "McPeaks")
#' mcatac <- normalize_egc(mcatac, "gastrulation.marginal")
#' mcatac <- add_const_peaks(mcatac, -16)
#' mcatac <- normalize_to_prob(mcatac)
#' }
#'
#' @export
normalize_to_prob <- function(mcatac, prob1_thresh = NULL, const_quantile = 0.8) {
    if (is.null(prob1_thresh)) {
        if (!has_name(mcatac@peaks, "const")) {
            cli_abort("No field named {.field const} in peaks. Please run {.code add_const_peaks}")
        }
        prob1_thresh <- quantile(mcatac@egc[mcatac@peaks$const, ], const_quantile)
        cli::cli_alert("Using {.val {prob1_thresh}} as the p=1 threshold")
    }
    mcatac@egc <- mcatac@egc / prob1_thresh
    mcatac@egc[mcatac@egc > 1] <- 1
    # mcatac@egc <- norm01(mcatac@egc)

    return(mcatac)
}
