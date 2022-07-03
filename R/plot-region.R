
#' Plot a genomic region
#'
#' @param mct an McTracks object.
#' @param intervals an intervals set with the genomic region to plot (a data frame with a single line). Note that if the start or end coordinates are not divisible by the resolution, the region will be extended to the next resolution interval.
#' @param metacells a vector of metacells to plot. If NULL, all metacells will be plotted.
#'
#' @inheritParams mct_get_mat
#' @inheritParams plot_region_mat
#'
#' @return None. The plot is displayed in a new device.
#'
#' @examples
#' \dontrun{
#' mct_plot_region(mct, intervals = data.frame(chr = "chr5", start = 54974252, end = 55074253))
#' }
#'
#' @export
mct_plot_region <- function(mct, intervals, downsample = TRUE, downsample_n = NULL, metacells = NULL, colors = c("white", "gray", "black", "gold", "gold")) {
    mat <- mct_get_mat(mct, intervals, downsample, downsample_n)
    mat <- mat[, intersect(mct@metacells[mct@order], colnames(mat)), drop = FALSE]

    if (!is.null(metacells)) {
        if (any(metacells %!in% mct@metacells)) {
            cli_abort("The following metacells are not in the McTracks object: {.val {metacells}}")
        }
        mat <- mat[, intersect(metacells, colnames(mat)), drop = FALSE]
    }

    if (has_cell_type(mct) && has_cell_type_colors(mct)) {
        mc_colors <- mct@metadata %>%
            select(metacell, color) %>%
            tibble::deframe()
        mc_colors <- mc_colors[colnames(mat)]
    } else {
        mc_colors <- NULL
    }
    plot_region_mat(mat, mc_colors, colors = colors, intervals = intervals)
}

#' Plot a genomic region given a matrix
#'
#' @param mat a matrix where rows are coordinates and columns are metacells
#' @param mc_colors a vector of colors for the metacells (optional)
#' @param colors color pallette for the ATAC signal
#' @param intervals the plotted intervals (optional)
#'
#' @export
plot_region_mat <- function(mat, mc_colors = NULL, colors = c("white", "gray", "black", "gold", "gold"), intervals = NULL) {
    mat_smooth <- apply(mat, 2, zoo::rollsum, 20, fill = "extend")

    layout(matrix(1:2, nrow = 1), w = c(1, 20))
    par(mar = c(4, 1, 2, 0))
    image(t(as.matrix(1:length(mc_colors), nrow = 1)), col = mc_colors, yaxt = "n", xaxt = "n")

    par(mar = c(4, 0, 2, 2))
    shades <- colorRampPalette(colors)(1000)
    image(mat_smooth, col = shades, yaxt = "n", xaxt = "n")
    if (!is.null(intervals)) {
        axis(1, at = seq(0, 1, l = 11), round(seq(intervals$start[1], intervals$end[1], l = 11) / 1e+6, 3), las = 2)
    }
}
