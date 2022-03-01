#' Construct a PeakIntervals object from misha intervals
#'
#' @description a PeakIntervals object is simply a misha intervals set which can have, optionally, a field called 'peak_name'
#' which holds the name of the peaks.
#'
#' @param intervals misha intervals (a data frame which has the following columns: 'chrom', 'start' and 'end')
#' @return a PeakIntervals object
#'
#' @export
PeakIntervals <- function(intervals) {
    validate_peaks(intervals)

    class(intervals) <- c("PeakIntervals", class(intervals))
    return(intervals)
}

validate_peaks <- function(peaks) {
    if (!is.data.frame(peaks)) {
        cli_abort("{.field peaks} should be a data frame")
    }

    if (!all(colnames(peaks)[1:3] == c("chrom", "start", "end"))) {
        cli_abort("{.field peaks} does not contain a valid intervals set (the first columns should be 'chrom', 'start' and 'end')")
    }

    distinct_peaks <- peaks %>% distinct(chrom, start, end)
    if (nrow(distinct_peaks) < nrow(peaks)) {
        cli_abort("Some peaks appear more than once.")
    }

    if (has_name(peaks, "peak_name")) {
        dups <- duplicated(peaks$peak_name)
        if (any(dups)) {
            dup_names <- paste(unique(peaks$peak_name[dups]), collapse = ", ")
            cli_abort("The following peak names are duplicated: {dup_names}")
        }
    }
}

#' Return the names of PeakIntervals
#'
#' @param peaks a PeakIntervals object or misha intervals
#'
#' @return the field called 'peak_name' if exists in \code{peaks}, or the coordinates seperated by an underscore.
#'
#' @export
peak_names <- function(peaks) {
    if (has_name(peaks, "peak_name")) {
        return(peaks$peak_name)
    }

    peak_names <- peaks %>%
        tidyr::unite("coord", start, end, sep = "-") %>%
        tidyr::unite("peak_name", chrom, coord, sep = ":") %>%
        pull(peak_name)

    return(peak_names)
}

#' Select peaks features by minimal coverage and threshold max-min fold
#'
#' @param atac a McATAC or ScATAC object
#' @param min_fold minimal max-min fold change
#' @param min_abs minimal total coverage per peak
#' @param only_promoters when TRUE - summarise over features to get promoters
#'
#' @inheritParams misha.ext::get_promoters
#'
#' @export
gen_peak_features <- function(atac, min_fold, min_abs, only_promoters = FALSE, upstream = 500, downstream = 50) {

}
