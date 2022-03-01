#' Construct a new McATAC object
#'
#' @param mat a numeric matrix where rows are peaks, and columns are metacells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations.
#' @description McATAC is a shallow object holding ATAC data over metacells.
#' Minimally it should include a count matrix of peaks over metacells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @export
McATAC <- function(mat,
                   peaks,
                   metadata = NULL) {
    return(make_atac_object(mat, peaks, metadata, "metacell", "McATAC"))
}

#' @export
#' @noRd
print.McATAC <- function(x, ...) {
    cli::cli_text("An McATAC object with {.val {ncol(x$mat)}} metacells {.val {nrow(x$mat)}} ATAC peaks.")
    cli::cli_text("Slots include:")
    cli::cli_dl()
    cli::cli_li(c("{.code $mat}" = "a numeric matrix where rows are peaks, and columns are metacells. Can be a sparse matrix."))
    cli::cli_li(c("{.code $peaks}" = "a misha intervals set with the peak definitions."))
    if (!is.null(x$metadata)) {
        cli::cli_li(c("{.code $metadata}" = "a tibble with a column called 'metacell' and additional metacell annotations."))
    }
}

#' Construct a new ScATAC object
#'
#' @param mat a numeric matrix where rows are peaks, and columns are cells. Can be a sparse matrix.
#' @param metadata data frame with a column called 'cell_id' and additional per-cell annotations.
#' @description ScATAC is a shallow object holding ATAC data over cells.
#' Minimally it should include a count matrix of peaks over cells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @inheritParams McATAC
#' @export
ScATAC <- function(mat,
                   peaks,
                   metadata = NULL) {
    return(make_atac_object(mat, peaks, metadata, "cell_id", "ScATAC"))
}

#' @export
#' @noRd
print.ScATAC <- function(x, ...) {
    cli::cli_text("An ScATAC object with {.val {ncol(x$mat)}} cells {.val {nrow(x$mat)}} ATAC peaks.")
    cli::cli_text("Slots include:")
    cli::cli_dl()
    cli::cli_li(c("{.code $mat}" = "a numeric matrix where rows are peaks, and columns are cells. Can be a sparse matrix."))
    cli::cli_li(c("{.code $peaks}" = "a misha intervals set with the peak definitions."))
    if (!is.null(x$metadata)) {
        cli::cli_li(c("{.code $metadata}" = "a tibble with a column called 'cell_id' and additional cell annotations."))
    }
}

make_atac_object <- function(mat, peaks, metadata, metadata_id_field, class_name) {
    peaks <- PeakIntervals(peaks)
    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
    }
    rownames(mat) <- peak_names(peaks)

    validate_atac_object(mat, peaks, metadata, metadata_id_field)

    if (!is.null(metadata)) {
        metadata <- as_tibble(metadata)
    }

    obj <- list(
        mat = mat,
        peaks = peaks,
        metadata = metadata
    )

    class(obj) <- class_name
    return(obj)
}


#' Validate a McATAC or ScATAC object
#'
#' @param metadata_id_field id field of the metadata. Should be "metacell" for McATAC object and "cell_id" for ScATAC object
#'
#' @inheritParams McATAC
#' @noRd
validate_atac_object <- function(mat, peaks, metadata, metadata_id_field) {
    if (!is.matrix(mat) && !methods::is(mat, "sparseMatrix")) {
        cli_abort("{.field mat} shuold be a matrix or a sparse matrix")
    }

    # make sure the matrix rownames are the peak names
    if (any(rownames(mat) != peak_names(peaks))) {
        cli_abort("rownames of the matrix are not the same as the peak names.")
    }

    if (!is.null(metadata)) {
        if (!is.data.frame(metadata)) {
            cli_abort("{.field metadata} is not a data frame")
        }
        if (!has_name(metadata, metadata_id_field)) {
            cli_abort("{.field metadata} doesn't have the required field {.field {metadata_id_field}}")
        }
        metadata <- as_tibble(metadata)

        # make sure that all cells/metacells exist within the matrix
        missing_cells <- metadata[[metadata_id_field]] %!in% colnames(mat)
        if (any(missing_cells)) {
            missing_cells <- paste(unique(metadata[[metadata_id_field]][missing_cells]), collapse = ", ")
            cli_abort("The following {metadata_id_field}s are missing from {.field mat} colnames: {.val {missing_cells}}")
        }
    }
}
