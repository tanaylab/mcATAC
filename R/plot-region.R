#' Plot a genomic region
#'
#' @param mct an McTracks object.
#' @param intervals an intervals set with the genomic region to plot (a data frame with a single line). Note that if the start or end coordinates are not divisible by the resolution, the region will be extended to the next resolution interval.
#' @param metacells a vector of metacells to plot. If NULL, all metacells will be plotted.
#' @param detect_dca mark regions with differential cluster accessibility (DCA)
#' @param hc an hclust object to use for the hierarchical clustering of the metacells when \code{detect_dcs=TRUE}. If NULL, metacells will be clustered using \code{mc_hclust_rna}.
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
#' intervs <- gintervals(5, 54974252, 55074253)
#' mct_plot_region(mct, intervs, gene_annot = TRUE)
#' mct_plot_region(mct, intervs, detect_dca = TRUE, gene_annot = TRUE)
#'
#' get_gene_promoter("GZMK", downstream = 1e5, upstream = 1e5, unify = TRUE) %>%
#'     mct_plot_region(mct, intervs, detect_dca = TRUE, gene_annot = TRUE)
#'
#' intervs %>%
#'     gintervals.zoom_in(4) %>%
#'     mct_plot_region(mct, ., detect_dca = TRUE, gene_annot = TRUE)
#' intervs %>%
#'     gintervals.zoom_out(2) %>%
#'     mct_plot_region(mct, ., detect_dca = TRUE, gene_annot = TRUE)
#' intervs %>%
#'     gintervals.shift_left(1e5) %>%
#'     mct_plot_region(mct, ., detect_dca = TRUE, gene_annot = TRUE)
#' intervs %>%
#'     gintervals.shift_right(1e5) %>%
#'     mct_plot_region(mct, ., detect_dca = TRUE, gene_annot = TRUE)
#' }
#'
#' @export
mct_plot_region <- function(mct, intervals, detect_dca = FALSE, downsample = TRUE, downsample_n = NULL, metacells = NULL, colors = c("white", "gray", "black", "gold"), color_breaks = c(0, 6, 12, 18, 24), hc = NULL, force_cell_type = TRUE, gene_annot = FALSE, n_smooth = 10, n_pixels = 1000, plot_x_axis_ticks = TRUE, gene_annot_pos = "top", flip = FALSE, genes_correlations = NULL, cor_colors=c("blue", "white", "white", "white", "red"), cor_color_breaks=c(-1,-0.05, 0, 0.05, 1), roll_mean = FALSE, ...) {
    gset_genome(mct@genome)
    raw_mat <- mct_get_mat(mct, intervals, downsample, downsample_n)
    if (!is.null(metacells)) {
        if (any(metacells %!in% mct@metacells)) {
            cli_abort("The following metacells are not in the McTracks object: {.val {metacells}}")
        }
        raw_mat <- raw_mat[, intersect(metacells, colnames(raw_mat)), drop = FALSE]
    }

    mat <- raw_mat[, intersect(mct@metacells[mct@order], colnames(raw_mat)), drop = FALSE]

    if (detect_dca && is.null(hc)) {
        if (!has_rna(mct)) {
            cli_abort("Cannot detect DCA without either an hclust object or RNA data.")
        }
        mct <- mct_subset_metacells(mct, colnames(mat))
        hc <- mc_hclust_rna(mct, force_cell_type = force_cell_type)
    }

    if (!is.null(hc)) {
        if (any(hc$label %!in% colnames(mat))) {
            missing_mcs <- setdiff(hc$label, colnames(mat))
            cli_warn("The following metacells are present in the hclust object, but are missing in the matrix (this is probably due to downsampling): {.val {missing_mcs}}")
            hc <- dendextend::prune(hc, missing_mcs)
        }
        mat <- mat[, hc$label]
    }

    dca_mat <- NULL
    if (detect_dca) {
        dca_mat <- mct_diff_access_on_hc(mat, hc = hc, ...)
        dca_mat <- dca_mat[, hc$order, drop = FALSE]
    }
    y_seps <- NULL
    if (!is.null(hc)) {
        mat <- mat[, hc$order, drop = FALSE]
        y_seps = cutree(hc, h=0.1)
    }

    if (has_cell_type(mct) && has_cell_type_colors(mct)) {
        mc_colors <- get_metacell_colors(mct@metadata)
        mc_colors <- mc_colors[colnames(mat)]
    } else {
        mc_colors <- NULL
    }
    if (!is.null(genes_correlations) && has_rna(mct)){
        overlapping_types = intersect(colnames(mct@rna_egc), colnames(mat))
        atac_egc = t(t(mat)/colSums(mat))
        atac_legc = log2(1e-5+atac_egc)
        rna_legc = log2(1e-5+mct@rna_egc)
        cors = tgstat::tgs_cor(
            t(atac_legc[,overlapping_types]), 
            t(rna_legc[toupper(genes_correlations),overlapping_types,drop=FALSE]), 
            pairwise.complete.obs = T
        )
        mat = cors
        mc_colors = rep("white", length(genes_correlations))
        colors = cor_colors
        color_breaks = cor_color_breaks
        y_seps <- NULL
        dca_mat <- NULL
    }
    plot_region_mat(mat, mc_colors, colors = colors, color_breaks = color_breaks, intervals = intervals, resolution = mct@resolution, dca_mat = dca_mat, y_seps=y_seps, n_smooth = n_smooth, gene_annot = gene_annot, n_pixels = n_pixels, genome = mct@genome, plot_x_axis_ticks = plot_x_axis_ticks, gene_annot_pos = gene_annot_pos, roll_mean = roll_mean, flip = flip)
}

#' Plot a genomic region given a matrix
#'
#' @param mat a matrix where rows are coordinates and columns are metacells
#' @param mc_colors a vector of colors for the metacells (optional)
#' @param colors color pallette for the ATAC signal
#' @param color_breaks a vector of breaks for the color palette
#' @param intervals the plotted intervals (optional)
#' @param resolution the resolution of the plotted intervals (optional)
#' @param dca_mat a matrix with the differential cluster accessibility (DCA) for the plotted regions (optional). Output of \code{mct_diff_access_on_hc}.
#' @param n_smooth number of genomic bins to use for smoothing the signal. Signal is smoothed by a rolling sum for each metacell (optional). Default is 20.
#' @param n_pixels number of pixels in the plot. The DCA regions would be extended by \code{ceiling(2 * nrow(mat) / n_pixels)} (optional).
#' @param gene_annot whether to add gene annotations; these annotations rely on the existence of the existence of an intervals set called "intervs.global.tss" and "intervs.global.exon" in the genome's misha directory. (optional)
#' @param genome the genome to use for the gene annotations (optional)
#' @param gene_annot_pos the position of the gene annotations ("top" or "bottom")
#' @param flip whether to flip the coordinates (optional)
#'
#' @export
plot_region_mat <- function(mat, mc_colors = NULL, colors = c("white", "gray", "black", "gold"), color_breaks = c(0, 6, 12, 18, 24), intervals = NULL, resolution = NULL, dca_mat = NULL, y_seps=NULL, y_seps_lty=2, y_seps_lwd=1, n_smooth = 10, n_pixels = 1000, gene_annot = FALSE, genome = NULL, plot_x_axis_ticks = TRUE, gene_annot_pos = "top", flip = FALSE, roll_mean=FALSE) {
    if(roll_mean){
        mat_smooth <- RcppRoll::roll_mean(mat, n = n_smooth, fill = c(0, 0, 0))
    } else {
        mat_smooth <- RcppRoll::roll_sum(mat, n = n_smooth, fill = c(0, 0, 0))
    }

    if (gene_annot) {
        if (is.null(intervals) || is.null(resolution)) {
            cli_abort("If gene annotations are requested, intervals and resolution must be specified")
        }

        # Different layouts depending on annotation position
        if (gene_annot_pos == "top") {
            layout(cbind(c(0, 0, 3), c(1, 2, 4)), widths = c(1, 20), heights = c(2, 0.5, 15))
        } else if (gene_annot_pos == "bottom") {
            layout(cbind(c(3, 0, 0, 0), c(4, 1, 2, 0)), widths = c(1, 20), heights = c(15, 2, 0.5, 0.5))
        }

        par(mar = c(0, 0, 1, 2))
        plot_tss_strip(intervals, flip = flip)
        
        if (gene_annot_pos == "top") {
            par(mar = c(0, 0, 0, 2))
        } else {
            par(mar = c(0, 0, 0, 2))
        }
        gene_annots <- make_gene_annot(intervals, resolution, genome)
        
        if (is.null(gene_annots[["exon_coords"]])) {
            image(as.matrix(rep(0, ncol(mat)), nrow = 1), col = c("white", "black"), breaks = c(-0.5, 0, 1), yaxt = "n", xaxt = "n", frame.plot = FALSE)
        } else {
            exon_mat <- as.matrix(gene_annots[["exon_coords"]])
            if (flip) {
                exon_mat <- exon_mat[nrow(exon_mat):1,, drop = FALSE]
            }

            image(exon_mat, col = c("white", "black"), breaks = c(-0.5, 0, 1), yaxt = "n", xaxt = "n", frame.plot = FALSE)
        }

        top_mar <- 0
        left_mar <- 2
    } else {
        #par(mar = c(0,0,0,0))
        layout(matrix(1:2, nrow = 1), w = c(1, 20))
        top_mar <- 0
        left_mar <- 2
    }

    if (plot_x_axis_ticks && gene_annot_pos == "bottom") {
        plot_x_axis_ticks <- FALSE
        cli_alert_warning("Gene annotations are at the bottom, so x-axis ticks are disabled.")
    }

    if (plot_x_axis_ticks) {
        bottom_mar <- 4
    } else {
        bottom_mar <- 0
    }

    par(mar = c(bottom_mar, left_mar, top_mar, 0))
    image(t(as.matrix(seq_along(mc_colors), nrow = 1)), col = mc_colors, yaxt = "n", xaxt = "n")

    par(mar = c(bottom_mar, 0, top_mar, 2))
    shades <- colorRampPalette(colors)(1000)
    mat_smooth[mat_smooth > max(color_breaks)] <- max(color_breaks)
    shades_breaks <- approx(color_breaks, n = 1001)$y
    if (flip) {
        mat_smooth <- mat_smooth[nrow(mat_smooth):1, , drop = FALSE]
    }
    image(mat_smooth, col = shades, breaks = shades_breaks, yaxt = "n", xaxt = "n")
    if (!is.null(intervals) && plot_x_axis_ticks) {
        axis(1, at = seq(0, 1, l = 11), round(seq(intervals$start[1], intervals$end[1], l = 11) / 1e+6, 3), las = 2)
    }

    if (!is.null(dca_mat)) {
        n_peak_smooth <- ceiling(2 * nrow(mat) / n_pixels)
        cli_alert_info("Extending DCAs with {.val {n_peak_smooth}} bins. This can be tweaked using the {.field n_pixels} paramater.")
        dca_mat <- extend_dca_mat(dca_mat, n_peak_smooth)
        dca_cols <- c(rgb(0, 0, 1, 0.4), rgb(0, 0, 1, 0.15), rgb(0, 0, 0, 0), rgb(1, 0, 0, 0.15), rgb(1, 0, 0, 0.4))
        if (flip) {
            dca_mat <- dca_mat[nrow(dca_mat):1, , drop = FALSE]
        }
        image(dca_mat, col = dca_cols, add = TRUE, breaks = c(-3, -1.5, -0.5, 0.5, 1.5, 3))
    }
    if (!is.null(y_seps)){
        n = length(colnames(mat_smooth))
        a = y_seps[colnames(mat_smooth)]
        y_seps = unname(sapply(split(seq_along(a), a),max))+1
        y_bords = seq(-1/(2*(n-1)), 1+1/(2*(n-1)), length=n+1)[y_seps]
        abline(h=y_bords, lty=y_seps_lty, lwd = y_seps_lwd, add=TRUE, col="grey")
    }
}

plot_tss_strip <- function(intervals, flip = FALSE) {
    plot( # empty plot
        x = intervals$start:intervals$end,
        y = c(0, rep(0.7, intervals$end - intervals$start)),
        ann = FALSE,
        bty = "n",
        type = "n",
        xaxt = "n",
        yaxt = "n"
    )


    tss_df <- gintervals.neighbors1("intervs.global.tss", intervals) %>%
        filter(dist == 0) %>%
        arrange(chrom, start, end, strand, geneSymbol) %>%
        #        distinct(geneSymbol, strand, .keep_all = TRUE) %>%
        select(chrom, tss = start, strand, gene = geneSymbol)

    if (flip) {
        # flip the coordinates relative to the start end end of the intervals
        tss_df$tss <- intervals$end - tss_df$tss + intervals$start
        tss_df$strand <- -1 * tss_df$strand
    }
    if (nrow(tss_df) > 0) {
        text(x = tss_df$tss, y = rep(0.35, length(tss_df$tss)), labels = tss_df$gene, las = 1, cex = 1)
        line_len <- (intervals$end - intervals$start) * 0.01
        plus_tss <- tss_df %>%
            filter(strand == 1) %>%
            pull(tss)
        minus_tss <- tss_df %>%
            filter(strand == -1) %>%
            pull(tss)


        segments(x0 = tss_df$tss, x1 = tss_df$tss, y0 = 0, y1 = 0.2)
        if (length(plus_tss) > 0) {
            arrows(x0 = plus_tss, x1 = plus_tss + line_len, y0 = 0.2, y1 = 0.2, length = 0.05)
        }

        if (length(minus_tss) > 0) {
            arrows(x0 = minus_tss, x1 = minus_tss - line_len, y0 = 0.2, y1 = 0.2, length = 0.05)
        }
    }
}

#' Extend each DCA using a rolling maximum
#'
#' @description This function extends each DCA (differential cluster accesability) using a rolling maximum of size \code{n_peak_smooth} on the rows of \code{dca_mat}.
#'
#' @noRd
extend_dca_mat <- function(dca_mat, n_peak_smooth) {
    dca_mat[dca_mat == 0] <- -100 # set all 0 values to -100 to avoid problems with rolling max
    dca_mat <- dca_mat + 3 # make all the values positive
    dca_mat <- RcppRoll::roll_max(dca_mat, n = n_peak_smooth, fill = c(-97, -97, -97))
    # set the values back to a scale of -2 to 2
    dca_mat <- dca_mat - 3
    dca_mat[dca_mat == -100] <- 0
    return(dca_mat)
}


#' Return a mask matrix with differential cluster accessibility (DCA)
#'
#' @description This function detects differential cluster accessibility (DCA) for each sub-tree of a given `hclust` object with hierarchical clustering of the metacells. \cr
#' For each coordinate (column) in the matrix, the function computes the log fold change between the subtree and the rest of the tree, and returns a mask matrix indicating whether the fold change is above or below the threshold(s). \cr
#' In order to avoid very long DCAs, the functions limits the size of a peak/trough to a fraction of all the metacells (\code{sz_frac_for_peak}).
#'
#'
#' @param mat a matrix where rows are coordinates and columns are metacells
#' @param hc an hclust object with the order of the metacells
#' @param sz_frac_for_peak maximal fraction of metacells in a peak or trough. Default is 0.25.
#' @param peak_lf_thresh1,peak_lf_thresh2,trough_lf_thresh1,trough_lf_thresh2 thresholds for the log fold change of the peaks and troughs
#' @param u_reg regularization factor
#'
#' @return a matrix the same dimensions of \code{mat}, with a mask of differential cluster accessibility (DCA). A value of 1 means the log fold change was above
#' \code{peak_lf_thresh1}, and a value of 2 means the log fold change was above \code{peak_lf_thresh2}. The same with -1 and -2 for troughs.
#'
#' @export
mct_diff_access_on_hc <- function(mat, hc, sz_frac_for_peak = 0.25, u_reg = 4, peak_lf_thresh1 = 1.2, peak_lf_thresh2 = 2, trough_lf_thresh1 = -1, trough_lf_thresh2 = -2, ...) {
    if (length(hc$order) != ncol(mat)) {
        cli_abort("The number of metacells in the matrix and the hclust object do not match.")
    }

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
