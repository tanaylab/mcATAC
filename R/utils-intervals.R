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
