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
        withr::defer(gsetroot(.misha$GROOT))
        gset_genome(genome)
        bad_intervals <- intervals %>%
            filter(chrom %!in% gintervals.all()$chrom)
        if (nrow(bad_intervals) > 0) {
            bad_chroms <- paste(unique(bad_intervals$chrom), collapse = ", ")
            intervals <- intervals %>%
                anti_join(bad_intervals, by = c("chrom", "start", "end"))
            cli_alert_warning("removed {.field {nrow(bad_intervals)}} peak{?s} from the following chromosome(s) which are missing from {.field {genome}}: {.file {bad_chroms}}")
        }

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
#' @param tad_based whether to name peaks based on TADs. If an intervals set named "intervs.global.tad_names" does not
#' exists - the function will return coordinate based names.
#' @param promoters are the peaks promoters? When true, the peak names (which are gene names) are returned.
#'
#' @return the field called 'peak_name' if exists in \code{peaks}, or a tad based peak name if \code{tad_based=TRUE} and the coordinates in
#' the form of "{chrom}:{start}-{end}" otherwise.
#'
#' @examples
#' \dontrun{
#' misha.ext::gset_genome("hg38")
#' head(peak_names(atac_sc@peaks)) # return the "peak_name" field
#' head(peak_names(atac_sc@peaks[, -4])) # return tad based peak names
#' head(peak_names(atac_sc@peaks[, -4], tad_based = FALSE)) # return coordinate based peak names
#' }
#'
#' @export
peak_names <- function(peaks, tad_based = TRUE, promoters = FALSE) {
    if (has_name(peaks, "peak_name")) {
        return(peaks$peak_name)
    } else if (promoters) {
        cli_abort("{.field peaks} does not have a 'peak_name' field. Please use {.code promoters = FALSE} to get the coordinates instead of the gene names.")
    }
    if (tad_based) {
        if (gintervals.exists("intervs.global.tad_names")) {
            peaks <- name_enhancers(peaks)
            peak_names <- peaks$peak_name
            return(peak_names)
        } else {
            cli_alert_warning("There are no TAD names for this assembly (check gintervals.exists('intervs.global.tad_names')). Using the coordinates instead.")
        }
    }
    peak_names <- peaks %>%
        tidyr::unite("coord", start, end, sep = "-") %>%
        tidyr::unite("peak_name", chrom, coord, sep = ":") %>%
        pull(peak_name)

    return(peak_names)
}

#' Get the names of the peaks close to a gene promoter
#'
#' @param peaks a PeakIntervals object
#' @param gene a gene name, e.g. "CD4"
#' @param max_dist_to_promoter_peak how far from \code{gene}'s TSS to search for promoter-proximal peaks. Default: 500
#' @param tss_intervals name of the intervals set containing the TSSs
#'
#' @return names of the peaks close to \code{gene}'s promoter. If no peaks are found, NULL is returned.
#'
#' @examples
#' \dontrun{
#' get_promoter_peaks(mc_atac@peaks, "CD4")
#' }
#' @export
get_promoter_peaks <- function(peaks, gene, max_dist_to_promoter_peak = 5e+2, tss_intervals = "intervs.global.tss") {
    assert_that(methods::is(peaks, "PeakIntervals"))
    if (!gintervals.exists(tss_intervals)) {
        cli_abort("{.val {tss_intervals}} intervals do not exist. You can either define the peak explicitly or import them from ucsc")
    }
    tss <- gintervals.load(tss_intervals)
    tss_gene <- tss[tss$geneSymbol == gene, c("chrom", "start", "end", "geneSymbol")]
    if (nrow(tss_gene) == 0) {
        return(NULL)
    } else if (nrow(tss_gene) > 1) {
        cli_alert("The gene {.val {gene}} has {.val {nrow(tss_gene)}} alternative promoters.")
    }

    nei_peaks_tss <- misha.ext::gintervals.neighbors1(
        tss_gene[, 1:3],
        as.data.frame(peaks),
        mindist = -max_dist_to_promoter_peak,
        maxdist = max_dist_to_promoter_peak,
        maxneighbors = 500
    ) %>%
        filter(!is.na(peak_name))

    peak <- unique(nei_peaks_tss$peak_name)

    return(peak)
}
