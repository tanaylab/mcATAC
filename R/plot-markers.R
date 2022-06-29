#' Plot a heatmap of marker genes fold change over metacells given an McATACPeaks object
#'
#' @description This is a thin wrapper around plot_rna_markers_mat which computes the marker matrix and uses the correct
#' metadata fields in the McATACPeaks object.
#'
#' @param atac_mc a McATACPeaks object with RNA expression (using \code{add_mc_rna})
#' @inheritParams get_rna_markers
#' @inheritParams get_rna_marker_matrix
#' @inheritDotParams plot_rna_markers_mat
#'
#' @examples
#' \dontrun{
#' plot_rna_markers(atac_mc)
#' }
#' @export
plot_rna_markers <- function(atac_mc, n_genes = 100, force_cell_type = TRUE, ...) {
    if (has_cell_type(atac_mc)) {
        metacell_types <- atac_mc@metadata
    } else {
        metacell_types <- NULL
    }

    if (has_cell_type_colors(atac_mc)) {
        cell_type_colors <- atac_mc@metadata %>% distinct(cell_type, color)
    } else {
        cell_type_colors <- NULL
    }

    plot_rna_markers_mat(
        get_rna_marker_matrix(atac_mc, force_cell_type = force_cell_type, n_genes = n_genes),
        metacell_types = metacell_types,
        cell_type_colors = cell_type_colors,
        ...
    )
}

#' Plot RNA and ATAC marker matrix side by side
#'
#' @description plot an rna marker heatmap (right) together with an heatmap of fold change of ATAC peaks that
#' are the most correlated to each marker gene (left, using \code{get_genes_atac_fp}).
#'
#' @return a \code{ggplot} object with the heatmap of RNA markers (right) and the heatmap of ATAC peaks (left).
#' The ATAC matrix is ordered by the RNA matrix.
#'
#' @param row_names show the row names of the heatmaps. Note that the space each plot takes might become uneven due to
#' genes / peaks with longer names.
#'
#' @return a ggplot object with the heatmaps
#'
#' @inheritParams get_rna_markers
#' @inheritParams get_rna_marker_matrix
#' @inheritParams plot_rna_markers_mat
#' @inheritDotParams plot_rna_markers_mat
#'
#' @examples
#' \dontrun{
#' plot_atac_rna_markers(atac_mc)
#' plot_atac_rna_markers(atac_mc, force_cell_type = FALSE)
#' plot_atac_rna_markers(atac_mc, n_genes = 200)
#' }
#'
#' @export
plot_atac_rna_markers <- function(atac_mc, n_genes = 100, force_cell_type = TRUE, plot_legend = TRUE, row_names = FALSE, ...) {
    m_rna <- get_rna_marker_matrix(atac_mc, n_genes = n_genes, force_cell_type = force_cell_type)
    cli_alert("Creating ATAC matrix by finding for each marker gene the ATAC peak that is most correlated to it.")
    m_atac <- get_genes_atac_fp(atac_mc, genes = rownames(m_rna), metacells = colnames(m_rna))

    if (has_cell_type(atac_mc)) {
        metacell_types <- atac_mc@metadata
    } else {
        metacell_types <- NULL
    }

    if (has_cell_type_colors(atac_mc)) {
        cell_type_colors <- atac_mc@metadata %>% distinct(cell_type, color)
    } else {
        cell_type_colors <- NULL
    }

    rna_plot <- plot_rna_markers_mat(
        m_rna,
        metacell_types = metacell_types,
        cell_type_colors = cell_type_colors,
        plot_legend = FALSE,
        row_names = row_names,
        title = "RNA",
        top_annotation_title = "",
        ...
    )

    atac_plot <- plot_rna_markers_mat(
        m_atac,
        metacell_types = metacell_types,
        cell_type_colors = cell_type_colors,
        plot_legend = FALSE,
        row_names = row_names,
        title = "ATAC",
        annotation_title = "",
        ...
    )

    p <- cowplot::plot_grid(atac_plot, rna_plot, nrow = 1)

    if (plot_legend) {
        legend_point_size <- max(1, min(10, 250 / nrow(cell_type_colors)))
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors)) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)))

        p <- cowplot::plot_grid(p, legend, nrow = 1, rel_widths = c(0.8, 0.15))
    }

    return(p)
}

#' Plot a heatmap of marker genes fold change over metacells
#'
#' @param mat a heatmap of marker genes fold change over metacells (e.g. output of \code{get_rna_marker_matrix})
#' @param metacell_types a data frame with a field called "metacell" and a field called "cell_type" (optional)
#' @param cell_type_colors a data frame with a field called "cell_type" and a field called "color" (optional)
#' @param low_color a color for low fold change
#' @param high_color a color for high fold change
#' @param mid_color a color for mid fold change
#' @param midpoint midpoint for the color scale (default: 0)
#' @param min_lfp minimum log2 fold change (default: -3). Value below this will be set to \code{min_lfp}.
#' @param max_lfp maximum log2 fold change (default: 3). Value above this will be set to \code{max_lfp}.
#' @param plot_legend TRUE to plot a legend of the cell types (default: TRUE)
#' @param top_cell_type_bar add a cell type annotation also at the top of the heatmap (default: TRUE)
#' @param gene_colors a named list with a color for each marker gene (optional)
#' @param col_names show column names (metacells)
#' @param interleave show the gene names (rows) interleaved
#' @param vertical_gridlines show vertical gridlines
#' @param annotation_title title for cell type annotation. Default: "Cell type"
#' @param top_annotation_title title for top cell type annotation. Default: "Cell type"
#' @param title title for the plot
#'
#' @return a ggplot object with the heatmap (and legend if \code{plot_legend} is TRUE)
#'
#' @inheritParams tgutil::tgplot_heatmap
#' @examples
#' \dontrun{
#' m_rna <- get_rna_marker_matrix(atac_mc)
#' plot_rna_markers_mat(m_rna, atac_mc@metadata, atac_mc@metadata, col_names = F)
#' }
#' @export
plot_rna_markers_mat <- function(mat,
                                 metacell_types = NULL,
                                 cell_type_colors = NULL,
                                 low_color = "blue",
                                 high_color = "red",
                                 mid_color = "white",
                                 midpoint = 0,
                                 min_lfp = -3,
                                 max_lfp = 3,
                                 plot_legend = TRUE,
                                 top_cell_type_bar = TRUE,
                                 gene_colors = NULL,
                                 col_names = FALSE,
                                 row_names = TRUE,
                                 col_names_orient = "slanted",
                                 plot_right = TRUE,
                                 interleave = TRUE,
                                 vertial_gridlines = FALSE,
                                 annotation_title = "Cell type",
                                 top_annotation_title = "Cell type",
                                 title = NULL) {
    gene_colors <- gene_colors %||% rep("black", nrow(mat))

    if (col_names) {
        col_names <- colnames(mat)
    }

    p_mat <- tgutil::tgplot_heatmap(
        tgutil::clip_vals(mat, min_lfp, max_lfp),
        col_names = col_names,
        col_names_orient = col_names_orient,
        plot_right = plot_right,
        interleave = interleave
    ) +
        scale_fill_gradient2(name = "", low = low_color, high = high_color, mid = mid_color, midpoint = midpoint, limits = c(min_lfp, max_lfp))

    if (row_names) {
        if (interleave) {
            p_mat <- p_mat +
                theme(
                    axis.text.y.right = ggtext::element_markdown(color = gene_colors[seq(2, length(gene_colors), 2)]),
                    axis.text.y.left =  ggtext::element_markdown(color = gene_colors[seq(1, length(gene_colors), 2)])
                )
        } else {
            p_mat <- p_mat +
                theme(axis.text.y = ggtext::element_markdown(color = gene_colors))
        }
    } else {
        p_mat <- p_mat + theme(axis.text.y = element_blank())
    }


    if (vertial_gridlines) {
        p_mat <- p_mat + geom_hline(yintercept = 1:nrow(mat) - 0.5, color = "gray", size = 0.1)
    }

    p_mat <- p_mat + theme(legend.position = "top")

    if (!is.null(title)) {
        p_mat <- p_mat + ggtitle(title)
    }

    if (is.null(metacell_types)) {
        return(p_mat)
    }

    if (is.null(cell_type_colors)) {
        cell_type_colors <- metacell_types %>%
            distinct(cell_type) %>%
            mutate(color = chameleon::distinct_colors(n())$name)
    }

    cell_type_colors <- cell_type_colors %>% distinct(cell_type, color)
    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% mutate(metacell = as.character(metacell)) %>% select(metacell, cell_type), by = "metacell") %>%
        left_join(cell_type_colors, by = "cell_type")

    p_mat <- p_mat %>% tgutil::tgplot_add_axis_annotation(mc_types$color, label = annotation_title, plot_left = FALSE)
    if (top_cell_type_bar) {
        p_mat <- p_mat %>% tgutil::tgplot_add_axis_annotation(mc_types$color, position = "top", label = top_annotation_title, plot_right = FALSE)
    }

    if (plot_legend) {
        legend_point_size <- max(1, min(10, 250 / nrow(cell_type_colors)))
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors)) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)))

        p <- cowplot::plot_grid(p_mat, legend, nrow = 1, rel_widths = c(0.8, 0.15))
        return(p)
    }

    return(cowplot::plot_grid(p_mat))
}
