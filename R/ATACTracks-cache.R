#' Get the key of ATAC track matrix cache
#'
#' @param mct a McTrack object
#' @param intervals an intervals set
#' @param downsample return a key for the downsampled matrix
#' @param downsample_n the number of counts to downsample to
#'
#' @return a key for the ATAC track matrix cache
#'
#' @noRd
intervs_key <- function(mct, intervals, downsample = FALSE, downsample_n = NULL) {
    key_int <- paste(intervals$chrom, intervals$start, intervals$end, sep = "_")
    key_tracks <- paste(mct@tracks, sep = "_")
    key <- paste0("tracks:", key_tracks, ";intervs:", key_int)
    key <- digest::digest(key, algo = "md5")
    if (downsample) {
        if (is.null(downsample_n)) {
            cli_abort("{.field downsample_n} must be provided in order to cache a downsampled matrix")
        }
        key <- paste0("downsample__{downsample_n}:", key)
    }
    return(key)
}

mct_has_region <- function(mct, intervals, downsample = FALSE, downsample_n = NULL) {
    !is.null(.mcatac_cache__[[intervs_key(mct, intervals, downsample, downsample_n)]])
}

mct_cache_mat <- function(mct, mat, intervals, downsample = FALSE, downsample_n = NULL) {
    .mcatac_cache__[[intervs_key(mct, intervals, downsample, downsample_n)]] <- mat
}

#' Clear the ATAC track matrix cache
#'
#' @return None.
#'
#' @examples
#' clean_cache()
#'
#' @export
clean_cache <- function() {
    rm(list = ls(envir = .mcatac_cache__), envir = .mcatac_cache__)
}

mct_cache_region <- function(mct, intervals, downsample = TRUE, downsample_n = NULL, force = FALSE, seed = NULL) {
    if (force || !mct_has_region(mct, intervals, FALSE)) {
        cli_alert("Extracting region")
        mat <- gextract(mct@tracks, intervals, colnames = mct@metacells, iterator = mct@resolution) %>%
            misha.ext::intervs_to_mat(remove_intervalID = TRUE)
        mat[is.na(mat)] <- 0
        mct_cache_mat(mct, mat, intervals, FALSE)
    }

    if (downsample) {
        if (force || !mct_has_region(mct, intervals, downsample = TRUE, downsample_n = downsample_n)) {
            mat <- .mcatac_cache__[[intervs_key(mct, intervals, FALSE)]]

            ds_mat <- downsamp_region(mat, mct@total_cov, downsample_n, seed = seed)

            removed_metacells <- setdiff(colnames(mat), colnames(ds_mat))
            cli_alert_warning("{.value {length(removed_metacells)}} metacells were removed from the downsampled matrix")

            mct_cache_mat(mct, ds_mat, intervals, TRUE, downsample_n)
        }
    }
}

#' Get an ATAC matrix from a McTracks object
#'
#' @description Retruns a matrix of genomic coordinates over the metacells. The matrix is cached in memory, so subsequent calls to this function will not require the extraction of the reads from the underlying tracks. \cr
#' When \code{downsample} is TRUE, downsampling is done by first setting a coverage goal,
#' then the reads within the region are subsampled relative to the ratio of the metacell coverage and the goal, i.e. \code{N_i * C_i/goal)} where \code{C_i} is the total number of reads in metacell i, \code{N_i} is the total numbers of reads in the region for metacell i, and \code{goal} is the coverage goal. \cr
#' For example, if the goal is 2M reads, a metacell has a total of 1M reads and 50K reads within the region, the reads within the region are subsampled to \code{50K * 1M/2M = 25K} (50%). Metacells with less than the goal are removed.
#'
#'
#' @param mct a McTrack object
#' @param intervlas an intervals set
#' @param downsample return a downsampled matrix. See description.
#' @param downsample_n total coverage goal. See description. Default: lower 5% percentile of the total coverage)
#' @param force force the computation of the matrix. If FALSE, the matrix is
#' retrieved from the cache if it exists.
#' @param seed random seed for the downsampling.
#'
#' @return a matrix where rows are genomic coordinates in the resolution
#' of \code{mct@resolution} and columns are the metacells in \code{mct@metacells}
#'
#' @export
mct_get_mat <- function(mct, intervals, downsample = FALSE, downsample_n = NULL, force = FALSE, seed = NULL) {
    if (downsample && is.null(downsample_n)) {
        downsample_n <- round(quantile(mct@total_cov, 0.05))
        cli_alert_info("Using 5% of the total coverage as the downsample threshold: {.value {downsample_n}}")
    }

    if (force || !mct_has_region(mct, intervals, downsample, downsample_n)) {
        mct_cache_region(mct, intervals, downsample, downsample_n, seed = seed)
    }

    return(.mcatac_cache__[[intervs_key(mct, intervals, downsample, downsample_n)]])
}


#' Downsample a region matrix relative to a total coverage goal
#'
#' See description of \code{mct_get_mat} for more information.
#'
#' @param mat a matrix of genomic coordinates over metacells
#' @param tr_total the total coverage of each metacell (not only at the region)
#' @param n_goal the target coverage. The number of reads within the region is subsampled to N * (n_goal /tr_total) where N is the total number of the metacell in the region.
#' @param seed random seed for the downsampling.
#'
#' @return a subsampled matrix, where metacells with less than n_goal reads are removed.
#'
#'
#' @noRd
downsamp_region <- function(mat, tr_total, n_goal, seed = NULL) {
    if (!is.null(seed)) {
        withr::local_seed(seed)
    }

    ratios <- n_goal / tr_total
    m <- nrow(mat)
    mat_ds <- NULL
    for (i in 1:ncol(mat)) {
        if (ratios[i] <= 1) {
            v <- mat[, i]
            n_i <- floor(sum(v) * ratios[i])
            new_column <- tabulate(sample(rep(1:length(v), times = v), replace = FALSE, size = n_i), nbins = m)
            if (is.null(mat_ds)) {
                mat_ds <- matrix(new_column, ncol = 1)
            } else {
                mat_ds <- cbind(mat_ds, new_column)
            }
        }
    }
    colnames(mat_ds) <- colnames(mat)[which(ratios <= 1)]
    rownames(mat_ds) <- rownames(mat)
    return(mat_ds)
}
