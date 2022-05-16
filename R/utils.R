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

#' @export
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
    ann_colors <- list(deframe(color_key))
    names(ann_colors)[[1]] <- feature_annotation
    return(list(col_annot, ann_colors))
}

#' Function to generate generic metacell annotation
#'
#' @description Generate a generic pheatmap-compatible annotation for metacells
#'
#' @param atac_mc (optional) a McATAC object to annotate
#' @param mc_clust (optional) a clustering of metacells, e.g. from gen_atac_mc_clust
#' @param k (optional) parameter for k-means clustering
#'
#' @return 2-element list containing a pheatmap-compatible dataframe and a list containing a named vector
#'
#' @noRd
generate_mc_annotation <- function(atac_mc, mc_clust = NULL, k = 10) {
    if (!is.null(atac_mc)) {
        if (all(has_name(atac_mc@metadata, c("metacell", "cell_type")))) {
            col_annot <- tibble::column_to_rownames(atac_mc@metadata[, c("metacell", "cell_type")], "metacell")
            ann_colors <- list("cell_type" = setNames(unlist(atac_mc@metadata[, "color"]), unlist(atac_mc@metadata[, "cell_type"])))
            mc_annot <- list(col_annot, ann_colors)
        } else {
            cli_alert_warning("McATAC object specified but no metacell annotation exists. Clustering metacells.")
            mc_clust <- gen_atac_mc_clust(atac_mc = atac_mc, use_prior_annot = F, k = k)
            mc_annot <- generate_pheatmap_annotation(mc_clust, feature_type = "metacell", feature_annotation = "cluster")
        }
    } else {
        if (!is.null(mc_clust)) {
            mc_annot <- generate_pheatmap_annotation(mc_clust, feature_type = "metacell", feature_annotation = "cluster")
        } else {
            cli_abort("Error: no McATAC object ({.var atac_mc}) and no metacell clustering ({.var mc_clust}) specified")
        }
    }
    return(mc_annot)
}

#' For each sparse matrix row compute the sum over a ragged array
#'
#' @param x a dgCMatrix sparse matrix
#' @param index a factor of the same length as the columns of x (would be coerced to a factor by \code{as.factor})
#'
#' @return A sparse dgCMatrix matrix of length(index) X nrow(x) size. Each ‘[i,j]’ element
#'    represents the ‘sum(x[i,which(index==levels(index)[j])])’.
#'
#'
#'
#' @noRd
sparse_matrix_tapply_sum <- function(x, index) {
    if (!is.factor(index)) {
        index <- factor(index)
    }
    groups <- levels(index)

    if (!is.null(colnames(x)) && !is.null(names(index))) {
        if (!all(colnames(x) == names(index))) {
            cli_warn("The column names of the input matrix do not match the index names. The index names would be ignored.")
        }
    }

    if (length(groups) == 1) {
        index_mat <- t(Matrix::Matrix(matrix(rep(1, length(index))), sparse = TRUE))
        rownames(index_mat) <- groups
    } else {
        index_mat <- t(Matrix::sparse.model.matrix(~ 0 + index))
        rownames(index_mat) <- groups
    }


    res <- t(index_mat %*% t(x))

    colnames(res) <- levels(index)
    if (!is.null(rownames(x))) {
        rownames(res) <- rownames(x)
    }
    return(res)
}

#' A helper function to deal with file overwrite
#'
#' @description if the file exists and overwrite is set to FALSE, the function will return an error, and if overwrite is set to TRUE the file will be deleted. If the file does not exist, nothing would happen.
#'
#' @param file name of the file to check
#' @param overwrite whether to overwrite the file
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' fn <- tempfile()
#' overwrite_file(fn, overwrite = FALSE) # this returns an error
#' overwrite_file(fn, overwrite = TRUE)
#' file.exists(fn) # returns FALSE
#' overwrite_file(fn)
#' }
#'
#' @noRd
overwrite_file <- function(file, overwrite) {
    if (file.exists(file)) {
        if (!overwrite) {
            cli_abort("File {.file {file}} already exists. Use 'overwrite = TRUE' to overwrite.")
        } else {
            unlink(file, recursive = TRUE)
        }
    }
}

#' Test if an executable file exists
#'
#' @description tests if an executable file exists by trying to run a command. If the return value is 0 the
#' function returns TRUE, otherwise FALSE.
#'
#' @param command a string with the command
#'
#' @return TRUE if the command returns 0, FALSE otherwise
#'
#' @inheritParams system2
#' @examples
#' bin_exists("echo", "test")
#' bin_exists("tabix", "--version")
#'
#' @noRd
bin_exists <- function(command, args = character(0)) {
    code <- system2(command, args, stdout = FALSE, stderr = FALSE)
    return(code == 0)
}

#' Check if intervals are within genome boundries
#'
#' @param intervals an intervals set
#' @param genome name of the genome (optional)
#'
#' @noRd
intervals_in_genome <- function(intervals, genome = NULL) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }

    if (any(intervals$chrom %!in% gintervals.all()$chrom)) {
        return(FALSE)
    }

    intervals <- intervals %>%
        select(chrom, start, end) %>%
        as.data.frame()

    if (gintervals.force_range(intervals) %>%
        anti_join(intervals, by = c("chrom", "start", "end")) %>%
        nrow() != 0) {
        return(FALSE)
    }

    return(TRUE)
}

#' Set parallel threads
#'
#' @description Set the number of parallel threads to use. mcATAC uses the R function \code{doMC::registerDoMC} to register the parallelization.
#' By default, mcATAC uses 80% of the number of available cores. The options are saved under 'mcatac.parallel' (should we use parallelization, logical) and 'mcatac.parallel.nc' (number of cores to use, integer).
#'
#' @param thread_num number of threads. use '1' for non parallel behavior
#'
#' @return None
#'
#' @examples
#' \donttest{
#' set_parallel(8)
#' }
#' @export
set_parallel <- function(thread_num = max(1, round(parallel::detectCores() * 0.8))) {
    if (thread_num <= 1) {
        options(mcatac.parallel = FALSE)
        cli_alert_info("Parallelization disabled.")
    } else {
        doMC::registerDoMC(thread_num)
        options(mcatac.parallel = TRUE)
        options(mcatac.parallel.nc = thread_num)
        cli_alert_info("Parallelization enabled. Using {.val {thread_num}} threads.")
    }
}
