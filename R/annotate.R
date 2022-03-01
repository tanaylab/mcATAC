
#' Annotate ATAC peaks
#'
#' @param atac a McATAC or ScATAC object
#'
#' @return the original \code{atac} object, with the peaks annotated using \code{annotate_intervals}
#'
#' @export
annotate_peaks <- function(atac) {
    if (!has_name("peaks", atac)) {
        cli_abort("{.field atac} doesn't have a field called {.field peaks}. You can annotate intervals directly using the {.code annotate_intervals} function.")
    }

    atac$peaks <- annotate_intervals(atac$peaks)

    return(atac)
}

#' Annotate an intervals set
#'
#' @param intervals the intervals set to annotate
#'
#' @return an intervals set with the following annotations (TODO)
#'
#' @export
annotate_intervals <- function(intervals) {
    # TODO: update peak annotation data frame with neighbor gene, distance, more (exons, enhancers etc.). The annotation might be an optional paramter.
}
