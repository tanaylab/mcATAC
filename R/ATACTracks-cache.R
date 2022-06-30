#' Get the key of ATAC track matrix cache
#'
#' @param mct a McTrack object
#' @param intervals an intervals set
#' @param downsample return a key for the downsampled matrix
#'
#' @return a key for the ATAC track matrix cache
#'
#' @noRd
intervs_key <- function(mct, intervals, downsample = FALSE) {
    key_int <- paste(intervals$chrom, intervals$start, intervals$end, sep = "_")
    key_tracks <- paste(mct@tracks, sep = "_")
    key <- paste0("tracks:", key_tracks, ";intervs:", key_int)
    key <- digest::digest(key, algo = "md5")
    if (downsample) {
        key <- paste0("downsample__:", key)
    }
    return(key)
}

mct_has_region <- function(mct, intervals, downsample = FALSE) {
    !is.null(.mcatac_cache__[[intervs_key(mct, intervals, downsample)]])
}

mct_cache_mat <- function(mct, mat, intervals, downsample = FALSE) {
    .mcatac_cache__[[intervs_key(mct, intervals, downsample)]] <- mat
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

mct_cache_region <- function(mct, intervals, downsample = TRUE, force = FALSE) {
    if (force || !mct_has_region(mct, intervals, FALSE)) {
        cli_alert("Extracting region")
        mat <- gextract(mct@tracks, intervals, colnames = mct@metacells, iterator = mct@resolution) %>%
            misha.ext::intervs_to_mat(remove_intervalID = TRUE)
        mat[is.na(mat)] <- 0
        mct_cache_mat(mct, mat, intervals, FALSE)
    }

    if (downsample) {
        if (force || !mct_has_region(mct, intervals, downsample = TRUE)) {
            mat <- .mcatac_cache__[[intervs_key(mct, intervals, FALSE)]]

            # ds_mat <-
            # mct_cache_mat(ds_mat, intervals, TRUE)
        }
    }
}

#' Get an ATAC matrix from a McTracks object
#'
#' @param mct a McTrack object
#' @param intervlas an intervals set
#' @param downsample return a downsampled matrix. Downsampling is done by
#' TODO
#' @param force force the computation of the matrix. If FALSE, the matrix is
#' retrieved from the cache if it exists.
#'
#' @return a matrix where rows are genomic coordinates in the resolution
#' of \code{mct@resolution} and columns are the metacells in \code{mct@metacells}
#'
#' @export
mct_get_mat <- function(mct, intervals, downsample = FALSE, force = FALSE) {
    if (force || !mct_has_region(mct, intervals, downsample)) {
        mct_cache_region(mct, intervals, downsample)
    }
    return(.mcatac_cache__[[intervs_key(mct, intervals, downsample)]])
}



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
            message("ds at ", i, " n_i = ", n_i, " rat = ", ratios[i])
            a <- tabulate(sample(rep(1:length(v), times = v), replace = F, size = n_i), nbins = m)
            if (is.null(mat_ds)) {
                mat_ds <- matrix(a, ncol = 1)
            } else {
                mat_ds <- cbind(mat_ds, a)
            }
        }
    }
    colnames(mat_ds) <- colnames(mat)[which(ratios <= 1)]
    return(mat_ds)
}
