#' Negation of the %in% operator
#'
#' @export
#' @noRd
`%!in%` <- Negate(`%in%`)


check_files_exist <- function(files) {
    for (file in files) {
        if (!file.exists(file) && !grepl("^http", file)) {
            cli_abort("{.file {file}} file doesn't exist.", call = parent.frame(1))
        }
    }
}

is_sparse_matrix <- function(mat) {
    return(methods::is(mat, "sparseMatrix"))
}

assert_atac_object <- function(obj, param = deparse(substitute(obj)), class = NULL) {
    if (is.null(class)) {
        if (!methods::is(obj, "ATAC")) {
            cli_abort("{.field {param}} must be an ScATAC or McATAC object", call = parent.frame(1))
        }
    } else {
        if (!methods::is(obj, class)) {
            cli_abort("{.field {param}} must be an {.field {class}} object", call = parent.frame(1))
        }
    }
}

#' Function to save pheatmaps to disk while showing them on screen
#'
#' @description \code{pheatmap::pheatmap} accepts a parameter called \code{filename} which
#' saves the pheatmap to disk. However, then the heatmap is not shown on screen.
#' This function is a workaround to show the heatmap on screen and save it to disk.
#'
#' @param dev name of the device to save the pheatmap. e.g. "png" or "pdf"
#' @inheritParams grDevices::png
#'
#' @export
save_pheatmap <- function(x, filename, dev = png, width = 2000, height = 2000, res = 150) {
    dev(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

#' Function to generate generic pheatmap annotation
#'
#' @description Generate a generic pheatmap-compatible annotation when clustering unannotated data and plotting with pheatmap
#'
#' @param clust_vec a clustering vector to generically annotate
#' @param feature_type (optional) type of the feature that was clustered
#' @param feature_name (optional) name of the feature that was clustered
#' @return 2-element list containing a pheatmap-compatible dataframe and a list containing a named vector

#' @noRd
generate_pheatmap_annotation <- function(clust_vec, feature_type = NULL, feature_annotation = NULL) {
    if (is.null(feature_type)) {
        feature_type <- "type"
    }
    if (is.null(feature_annotation)) {
        feature_annotation <- "annotation"
    }
    cts <- unique(clust_vec)

    color_key <- tibble(name = cts, color = chameleon::distinct_colors(length(cts))$name) %>%
        rename(!!feature_annotation := name)
    col_annot <- clust_vec %>%
        enframe(feature_type, feature_annotation) %>%
        column_to_rownames(feature_type)
    ann_colors <- list(feature_annotation = deframe(color_key))
    return(list(col_annot, ann_colors))
}

#' Function to generate generic metacell annotation
#'
#' @description Generate a generic pheatmap-compatible annotation for metacells
#'
#' @param mc_atac (optional) a McATAC object to annotate
#' @param mc_clust (optional) a clustering of metacells, e.g. from gen_atac_mc_clust
#' @param k (optional) parameter for k-means clustering
#'
#' @return 2-element list containing a pheatmap-compatible dataframe and a list containing a named vector
#'
#' @noRd
generate_mc_annotation <- function(mc_atac, mc_clust = NULL, k = 10) {
    if (!is.null(mc_atac)) {
        if (all(has_name(mc_atac@metadata, c("metacell", "cell_type")))) {
            col_annot <- tibble::column_to_rownames(mc_atac@metadata[, c("metacell", "cell_type")], "metacell")
            ann_colors <- list("cell_type" = setNames(unlist(mc_atac@metadata[, "color"]), unlist(mc_atac@metadata[, "cell_type"])))
            mc_annot <- list(col_annot, ann_colors)
        } else {
            cli_alert_warning("McATAC object specified but no metacell annotation exists. Clustering metacells.")
            mc_clust <- gen_atac_mc_clust(atac_mc = mc_atac, use_prior_annot = F, k = k)
            mc_annot <- generate_pheatmap_annotation(mc_clust, feature_type = "metacell", feature_annotation = "cluster")
        }
    } else {
        if (!is.null(mc_clust)) {
            mc_annot <- generate_pheatmap_annotation(mc_clust, feature_type = "metacell", feature_annotation = "cluster")
        } else {
            cli_abort("Error: no McATAC object ({.var mc_atac}) and no metacell clustering ({.var mc_clust}) specified")
        }
    }
    return(mc_annot)
}
