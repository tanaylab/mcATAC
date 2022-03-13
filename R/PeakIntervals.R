#' Construct a PeakIntervals object from misha intervals
#'
#' @description a PeakIntervals object is simply a misha intervals set which can have, optionally, a field called 'peak_name'
#' which holds the name of the peaks.
#'
#' Note that if you have informative rownames you should make them
#' an explicit column due to the fact that the package might change
#' them occasionally (e.g. when creating h5ad files)
#'
#' @param intervals misha intervals (a data frame which has the following columns: 'chrom', 'start' and 'end')
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10".
#' @return a PeakIntervals object
#'
#' @export
PeakIntervals <- function(intervals, genome = NULL) {
    validate_peaks(intervals)

    if (!is.null(genome)) {
        withr::defer(gsetroot(GROOT))
        gset_genome(genome)
        bad_intervals <- intervals %>%
            filter(chrom %!in% gintervals.all()$chrom)
        if (nrow(bad_intervals) > 0) {
            bad_chroms <- paste(unique(bad_intervals$chrom), collapse = ", ")
            intervals <- intervals %>%
                anti_join(bad_intervals, by = c("chrom", "start", "end"))
            cli_alert_warning("removed {.field {nrow(bad_intervals)}} peak{?s} from the following chromosome(s) which are missing from {.field {genome}}: {.file {bad_chroms}}")
        }
        intervals$end[1] <- gintervals.all()$end[1] + 100
        intervals$end[2] <- gintervals.all()$end[1] + 1e9
        out_intervals <- intervals %>%
            anti_join(
                intervals %>%
                    as.data.frame() %>%
                    gintervals.force_range(),
                by = c("chrom", "start", "end")
            )
        if (nrow(out_intervals) > 0) {
            intervals <- intervals %>%
                anti_join(out_intervals, by = c("chrom", "start", "end"))
            cli_alert_warning("removed {.field {nrow(out_intervals)}} peak{?s} which {?was/were} outside of {.field {genome}} genome.")
        }
    }

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

    assert_that(is.numeric(peaks$start))
    assert_that(is.numeric(peaks$end))

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
