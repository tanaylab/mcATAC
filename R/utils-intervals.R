#' Zoom-in for an intervals set.
#'
#' @param intervals An intervals set.
#' @param zoom the zoom-in factor.
#'
#' @return A new intervals set where the coordinates are \code{c + round(l/2)} and \code{c - round(l/2)} where \code{l} is the length of the interval divided by \code{zoom}, and \code{c} is the center of the interval.
#'
#' @examples
#' gdb.init_examples()
#' gintervals(1, 0, 200) %>% gintervals.zoom_in(2)
#'
#' @export
gintervals.zoom_in <- function(intervals, zoom) {
    clean_intervs <- intervals %>% select(chrom, start, end)
    clean_intervs <- clean_intervs %>%
        mutate(
            l = end - start,
            midpoint = start + round(l / 2),
            new_l = round(l / zoom),
            new_start = midpoint - round(new_l / 2),
            new_end = midpoint + round(new_l / 2)
        ) %>%
        select(chrom, start = new_start, end = new_end)

    intervals <- intervals %>%
        mutate(
            start = clean_intervs$start,
            end = clean_intervs$end
        )
    return(intervals)
}

#' Zoom-out for an intervals set.
#'
#' @param intervals An intervals set.
#' @param zoom the zoom-our factor.
#'
#' @return A new intervals set where the coordinates are \code{c + round(l/2)} and \code{c - round(l/2)} where \code{l} is the length of the interval multiplied by \code{zoom}, and \code{c} is the center of the interval. Note that the intervals remain bounded by the chromosome boundaries.
#'
#' @examples
#' gdb.init_examples()
#' gintervals(1, 50, 150) %>% gintervals.zoom_out(2)
#'
#' @export
gintervals.zoom_out <- function(intervals, zoom) {
    clean_intervs <- intervals %>% select(chrom, start, end)
    clean_intervs <- clean_intervs %>%
        mutate(
            l = end - start,
            midpoint = start + round(l / 2),
            new_l = round(l * zoom),
            new_start = midpoint - round(new_l / 2),
            new_end = midpoint + round(new_l / 2)
        ) %>%
        select(chrom, start = new_start, end = new_end) %>%
        gintervals.force_range()

    intervals <- intervals %>%
        mutate(
            start = clean_intervs$start,
            end = clean_intervs$end
        )
    return(intervals)
}

#' Shift the coordinates of an intervals set to the left.
#'
#' @param intervals An intervals set.
#' @param shift number of bp to shift to the left.
#'
#' @examples
#' gdb.init_examples()
#' gintervals(1, 50, 150) %>% gintervals.shift_left(20)
#'
#' @return A new intervals set where the coordinates are \code{start - shift} and \code{end - shift}. Note that the intervals remain bounded by the chromosome boundaries.
#'
#' @export
gintervals.shift_left <- function(intervals, shift) {
    clean_intervs <- intervals %>% select(chrom, start, end)
    clean_intervs <- clean_intervs %>%
        mutate(
            new_start = start - shift,
            new_end = end - shift
        ) %>%
        select(chrom, start = new_start, end = new_end) %>%
        gintervals.force_range()

    intervals <- intervals %>%
        mutate(
            start = clean_intervs$start,
            end = clean_intervs$end
        )
    return(intervals)
}

#' Shift the coordinates of an intervals set to the right.
#'
#' @param intervals An intervals set.
#' @param shift number of bp to shift to the right.
#'
#' @examples
#' gdb.init_examples()
#' gintervals(1, 50, 150) %>% gintervals.shift_right(20)
#'
#' @return A new intervals set where the coordinates are \code{start + shift} and \code{end + shift}. Note that the intervals remain bounded by the chromosome boundaries.
#'
#' @export
gintervals.shift_right <- function(intervals, shift) {
    clean_intervs <- intervals %>% select(chrom, start, end)
    clean_intervs <- clean_intervs %>%
        mutate(
            new_start = start + shift,
            new_end = end + shift
        ) %>%
        select(chrom, start = new_start, end = new_end) %>%
        gintervals.force_range()

    intervals <- intervals %>%
        mutate(
            start = clean_intervs$start,
            end = clean_intervs$end
        )
    return(intervals)
}

#' Get an intervals set of a gene promoter.
#'
#' @param name The name of the gene ("name2" in UCSC files).
#' @param unify unify the intervals if there are multiple promoters. The start coordinate of the intervals will be the minimum of the starts and the end coordinate of the intervals will be the maximum of the ends.
#' @inheritParams misha.ext::get_promoters
#'
#' @return an intervals set with the promoter of the gene. If \code{unify} is TRUE - the intervals will be unified to a single row.
#'
#' @examples
#' \dontrun{
#' get_gene_promoter("GZMK")
#' }
#'
#' @export
get_gene_promoter <- function(gene, upstream = 500, downstream = 50, unify = FALSE) {
    if (!gintervals.exists("intervs.global.tss")) {
        cli_abort("Intervals set named {.val intervs.global.tss} not found.")
    }
    promoters <- misha.ext::get_promoters(upstream = upstream, downstream = downstream)
    if (!has_name(promoters, "geneSymbol")) {
        cli_abort("Gene symbol not found in the {.val intervs.global.tss} intervals set.")
    }
    res <- promoters %>%
        filter(geneSymbol == gene) %>%
        select(chrom, start, end)

    if (nrow(res) == 0) {
        cli_abort("Gene {.val {gene}} not found in the {.val intervs.global.tss} intervals set.")
    }

    if (unify) {
        res <- res %>%
            summarise(
                chrom = chrom[1],
                start = min(start),
                end = max(end)
            )
    }

    return(res)
}
