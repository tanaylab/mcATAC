
#' Plot a genomic region
#'
#' @param mct an McTracks object.
#' @param intervals an intervals set with the genomic region to plot (a data frame with a single line). Note that if the start or end coordinates are not divisible by the resolution, the region will be extended to the next resolution interval.
#' @param metacells a vector of metacells to plot. If NULL, all metacells will be plotted.
#' @param detect_dca mark regions with differential cluster accessibility (DCA)
#' @param hc an hclust object to use for the hierarchical clustering of the metacells. If NULL, metacells will be clustered using \code{mc_hclust_rna}.
#'
#' @inheritDotParams mct_diff_access_on_hc
#' @inheritParams mct_get_mat
#' @inheritParams plot_region_mat
#' @inheritParams mc_hclust_rna
#'
#'
#' @return None. The plot is displayed in a new device.
#'
#' @examples
#' \dontrun{
#' mct_plot_region(mct, intervals = data.frame(chr = "chr5", start = 54974252, end = 55074253))
#' }
#'
#' @export
mct_plot_region <- function(mct, intervals, detect_dca = FALSE, downsample = TRUE, downsample_n = NULL, metacells = NULL, colors = c("white", "gray", "black", "gold", "gold"), hc = NULL, force_cell_type = TRUE, n_pixels = 1000, ...) {
    raw_mat <- mct_get_mat(mct, intervals, downsample, downsample_n)
    if (!is.null(metacells)) {
        if (any(metacells %!in% mct@metacells)) {
            cli_abort("The following metacells are not in the McTracks object: {.val {metacells}}")
        }
        raw_mat <- raw_mat[, intersect(metacells, colnames(raw_mat)), drop = FALSE]
    }

    mat <- raw_mat[, intersect(mct@metacells[mct@order], colnames(raw_mat)), drop = FALSE]

    dca_mat <- NULL
    if (detect_dca) {
        if (is.null(hc)) {
            mct <- mct_subset_metacells(mct, colnames(mat))
            hc <- mc_hclust_rna(mct, force_cell_type = force_cell_type)
        }
        mat <- mat[, hc$label]
        dca_mat <- mct_diff_access_on_hc(mat, hc = hc, ...)
        dca_mat <- dca_mat[, hc$order, drop = FALSE]
        mat <- mat[, hc$order, drop = FALSE]
    }

    if (has_cell_type(mct) && has_cell_type_colors(mct)) {
        mc_colors <- mct@metadata %>%
            select(metacell, color) %>%
            tibble::deframe()
        mc_colors <- mc_colors[colnames(mat)]
    } else {
        mc_colors <- NULL
    }

    plot_region_mat(mat, mc_colors, colors = colors, intervals = intervals, dca_mat = dca_mat, n_pixels = n_pixels)
}

#' Plot a genomic region given a matrix
#'
#' @param mat a matrix where rows are coordinates and columns are metacells
#' @param mc_colors a vector of colors for the metacells (optional)
#' @param colors color pallette for the ATAC signal
#' @param intervals the plotted intervals (optional)
#' @param dca_mat a matrix with the differential cluster accessibility (DCA) for the plotted regions (optional). Output of \code{mct_diff_access_on_hc}.
#' @param n_pixels number of pixels in the plot. The DCA regions would be extended by \code{ceiling(2 * nrow(mat) / n_pixels)}.
#'
#' @export
plot_region_mat <- function(mat, mc_colors = NULL, colors = c("white", "gray", "black", "gold", "gold"), intervals = NULL, dca_mat = NULL, n_pixels = 1000) {
    mat_smooth <- RcppRoll::roll_sum(mat, n = 20, fill = c(0, 0, 0))

    layout(matrix(1:2, nrow = 1), w = c(1, 20))
    par(mar = c(4, 1, 2, 0))
    image(t(as.matrix(1:length(mc_colors), nrow = 1)), col = mc_colors, yaxt = "n", xaxt = "n")

    par(mar = c(4, 0, 2, 2))
    shades <- colorRampPalette(colors)(1000)
    image(mat_smooth, col = shades, yaxt = "n", xaxt = "n")
    if (!is.null(intervals)) {
        axis(1, at = seq(0, 1, l = 11), round(seq(intervals$start[1], intervals$end[1], l = 11) / 1e+6, 3), las = 2)
    }

    if (!is.null(dca_mat)) {
        n_peak_smooth <- ceiling(2 * nrow(mat) / n_pixels)
        cli_alert_info("Extending DCAs with {.val {n_peak_smooth}} bins. This can be tweaked using the {.field n_pixels} paramater.")
        dca_mat[dca_mat == 0] <- -100
        dca_mat <- dca_mat + 3
        dca_mat <- RcppRoll::roll_max(dca_mat, n = n_peak_smooth, fill = c(-97, -97, -97))
        dca_mat <- dca_mat - 3
        dca_mat[dca_mat == -100] <- 0
        dca_cols <- c(rgb(0, 0, 1, 0.4), rgb(0, 0, 1, 0.15), rgb(0, 0, 0, 0), rgb(1, 0, 0, 0.15), rgb(1, 0, 0, 0.4))
        image(dca_mat, col = dca_cols, add = TRUE, breaks = c(-3, -1.5, -0.5, 0.5, 1.5, 3))
    }
}


#' Return a mask matrix with differential cluster accessibility (DCA)
#'
#' @param mat a matrix where rows are coordinates and columns are metacells
#' @param hc an hclust object with the order of the metacells
#' @param sz_frac_for_peak maximal fraction of metacells in a peak
#' @param peak_lf_thresh1,peak_lf_thresh2,trough_lf_thresh1,trough_lf_thresh2 thresholds for the log fold change of the peaks and troughs
#' @param u_reg regularization factor
#'
#' @return a matrix the same dimensions of \code{mat}, with a mask of differential cluster accessibility (DCA). A value of 1 means the log fold change was above
#' \code{peak_lf_thresh1}, and a value of 2 means the log fold change was above \code{peak_lf_thresh2}. The same with -1 and -2 for troughs.
#'
#' @export
mct_diff_access_on_hc <- function(mat, hc, sz_frac_for_peak = 0.25, u_reg = 4, peak_lf_thresh1 = 1, peak_lf_thresh2 = 2, trough_lf_thresh1 = -1,
                                  trough_lf_thresh2 = -2) {
    dca_mat1 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    dca_mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))

    stat_total <- rowSums(mat, na.rm = T)
    n_mc <- ncol(mat)

    hmat <- matrix(nrow = nrow(mat), ncol = 0)
    node_size <- numeric(0)
    subtree <- list()
    for (i in 1:(nrow(hc$merge) - 1)) {
        left <- hc$merge[i, 1]
        right <- hc$merge[i, 2]
        if (left < 0) {
            left_stat <- mat[, -left]
            left_size <- 1
            left_nodes <- c(-left)
        } else {
            left_stat <- hmat[, left]
            left_size <- node_size[left]
            left_nodes <- subtree[[left]]
        }
        if (right < 0) {
            right_stat <- mat[, -right]
            right_size <- 1
            right_nodes <- c(-right)
        } else {
            right_stat <- hmat[, right]
            right_size <- node_size[right]
            right_nodes <- subtree[[right]]
        }
        sz <- left_size + right_size
        stat <- right_stat + left_stat
        back_sz <- n_mc - sz
        back_stat <- stat_total - stat
        hmat <- cbind(hmat, stat)
        nodes <- unlist(c(left_nodes, right_nodes))
        node_size[i] <- left_size + right_size
        subtree[[i]] <- nodes

        reg <- u_reg / min(sz, back_sz)
        lf <- log2((reg + stat / sz) / (reg + back_stat / back_sz))
        if (sz < sz_frac_for_peak * n_mc) {
            if (sum(lf > peak_lf_thresh1) > 0) {
                dca_mat1[lf > peak_lf_thresh1, nodes] <- 1
            }
            if (sum(lf > peak_lf_thresh2) > 0) {
                dca_mat2[lf > peak_lf_thresh2, nodes] <- 2
            }
            if (sum(lf < trough_lf_thresh1) > 0) {
                dca_mat1[lf < trough_lf_thresh1, nodes] <- -1
            }
            if (sum(lf < trough_lf_thresh2) > 0) {
                dca_mat2[lf < trough_lf_thresh2, nodes] <- -2
            }
        }
    }
    dca_mat2[dca_mat2 == 0] <- dca_mat1[dca_mat2 == 0]
    colnames(dca_mat2) <- colnames(mat)
    rownames(dca_mat2) <- rownames(mat)
    return(dca_mat2)
}
