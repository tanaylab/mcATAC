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
