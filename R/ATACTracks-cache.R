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
    if (downsample && is.null(downsample_n)) {
        cli_abort("{.field downsample_n} must be provided in order to cache a downsampled matrix")
    }

    intervals <- intervals %>%
        mutate(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end)) %>%
        select(chrom, start, end)

    params <- list(
        tracks = mct@tracks,
        intervals = intervals,
        downsample = downsample,
        downsample_n = downsample_n
    )

    hash <- c(
        tracks = TRUE,
        intervals = FALSE,
        downsample = FALSE,
        downsample_n = FALSE
    )

    if (!downsample) {
        params[["downsample_n"]] <- NULL
        params[["downsample"]] <- NULL
        hash <- hash[1:2]
    }

    return(paste0("region;", cache_key(params, hash)))
}

#' Generate mcATAC cache key
#'
#' @description The generated key is of the form: "{name}:{value};" for all elements of the list.
#'
#' @param params a named list of parameters to be used in the key. Values should be coercible to strings.
#' Elements with NULL values will be ignored.
#' @param hash a logical vector the length of \code{params} where names are the parameter names and values indicate if the element should be hashed prior to
#' insertion to the key.
#'
#' @return a key for mcATAC cache
#'
#' @noRd
cache_key <- function(params, hash = NULL) {
    params <- params %>% purrr::discard(is.null)
    if (is.null(hash)) {
        hash <- rep(TRUE, length(params))
    }
    keys <- purrr::pmap_chr(
        list(params, names(params)), function(val, name) {
            key <- paste(as.character(val), collapse = "_")
            if (hash[name]) {
                key <- digest::digest(key, algo = "md5")
            }
            key <- paste0(name, ":", key)
        }
    )

    key <- paste(keys, collapse = ";")
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
#' clear_cache()
#'
#' @export
clear_cache <- function() {
    rm(list = ls(envir = .mcatac_cache__), envir = .mcatac_cache__)
}

#' List the ATAC track matrix cache
#'
#' @return a vector with all the cache keys.
#'
#' @noRd
ls_cache <- function(...) {
    ls(envir = .mcatac_cache__, ...)
}

#' List the cached regions of an McTracks object
#'
#' @param mct a McTracks object
#'
#' @return an intervals set with additional "downsample" and "downsample_n" fields, representing the intervals matrices that are currently cached.
#'
#' @examples
#' \dontrun{
#' mct_ls_cached_regions(mct)
#' }
#'
#' @export
mct_ls_cached_regions <- function(mct) {
    tracks_hash <- paste(as.character(mct@tracks), collapse = "_") %>% digest::digest(algo = "md5")
    prefix <- glue("^region;tracks:{tracks_hash}")
    cache_keys <- ls_cache(pattern = prefix)
    cache <- gsub(prefix, "", cache_keys)
    cache <- gsub("^;", "", cache)

    if (length(cache) == 0) {
        return(data.frame(chrom = character(0), start = numeric(0), end = numeric(0), downsample = logical(0), downsample_n = numeric(0)))
    }

    # This ugly code parses the cache key to extract the intervals
    all_intervals <- stringr::str_split(cache, pattern = ";") %>%
        map(~ map(stringr::str_split(.x, pattern = ":"), function(y) {
            l <- list()
            l[[y[1]]] <- y[2]
            l
        })) %>%
        map(as.data.frame) %>%
        purrr::map_dfr(~.x) %>%
        tidyr::separate(intervals, c("chrom", "start", "end"), sep = "_") %>%
        mutate(start = as.numeric(start), end = as.numeric(end))


    if (!has_name(all_intervals, "downsample")) {
        all_intervals <- all_intervals %>% mutate(downsample = FALSE)
    }
    if (!has_name(all_intervals, "downsample_n")) {
        all_intervals <- all_intervals %>% mutate(downsample_n = NA)
    }

    all_intervals <- all_intervals %>%
        mutate(downsample = as.logical(downsample), downsample_n = as.numeric(downsample_n)) %>%
        tidyr::replace_na(replace = list(downsample = FALSE))

    return(all_intervals)
}

mct_cache_region_ds <- function(mct, intervals, mat, downsample_n, seed) {
    ds_mat <- downsamp_region(mat, mct@total_cov, downsample_n, seed = seed)

    removed_metacells <- setdiff(colnames(mat), colnames(ds_mat))
    cli_alert_warning("{.value {length(removed_metacells)}} metacells were removed from the downsampled matrix")

    mct_cache_mat(mct, ds_mat, intervals, TRUE, downsample_n)
}

mct_cache_region <- function(mct, intervals) {
    intervals <- intervals[1, ]
    cli_alert("Extracting region {.val {intervals$chrom}}:{.val {intervals$start}}-{.val {intervals$end}}")
    mat <- gextract(mct@tracks, intervals, colnames = mct@metacells, iterator = mct@resolution) %>%
        misha.ext::intervs_to_mat(remove_intervalID = TRUE)
    mat[is.na(mat)] <- 0
    mct_cache_mat(mct, mat, intervals, FALSE)
    invisible(mat)
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
#' @param intervals an intervals set with a single interval. Note that if the start or end coordinates are not divisible by the resolution, the region will be extended to the next resolution interval.
#' @param downsample return a downsampled matrix. See description.
#' @param downsample_n total coverage goal. See description. Default: lower 5th percentile of the total coverage)
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
        cli_alert_info("Using 5th percentile of the total coverage as the downsample threshold: {.val {downsample_n}}")
    }

    intervals <- intervals %>%
        mutate(start = floor(start / mct@resolution) * mct@resolution, end = ceiling(end / mct@resolution) * mct@resolution) %>%
        gintervals.force_range()

    if (!force) {
        # if intervals are already cached as a full region - return it
        if (mct_has_region(mct, intervals, downsample = downsample, downsample_n = downsample_n)) {
            mat <- .mcatac_cache__[[intervs_key(mct, intervals, downsample = downsample, downsample_n = downsample_n)]]
            return(mat)
        }

        # check if intervals are within a cached region
        cached_regs <- mct_ls_cached_regions(mct)
        if (downsample) {
            cached_regs <- cached_regs %>% filter(downsample, downsample_n == !!downsample_n)
        }
        nei_regs <- intervals %>%
            gintervals.neighbors1(cached_regs) %>%
            filter(dist == 0) %>%
            filter(start >= cached_regs$start, end <= cached_regs$end)

        # If intervals are within a region - return the sub-matrix
        if (nrow(nei_regs) > 0) {
            region <- nei_regs %>%
                slice(1) %>%
                select(chrom = chrom1, start = start1, end = end1, downsample = downsample, downsample_n = downsample_n)
            full_mat <- .mcatac_cache__[[intervs_key(mct, region %>% select(chrom, start, end), downsample = downsample, downsample_n = downsample_n)]]
            region_l <- region$end - region$start
            sbin <- (intervals$start - region$start) %/% mct@resolution
            ebin <- nrow(full_mat) - (region$end - intervals$end) %/% mct@resolution
            mat <- full_mat[sbin:ebin, ]
            return(mat)
        }
    }

    # Create the matrices and cache them
    mat <- mct_cache_region(mct, intervals)
    if (downsample) {
        mat <- mct_cache_region_ds(mct, intervals, mat, downsample_n, seed = seed)
    }
    return(mat)
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
