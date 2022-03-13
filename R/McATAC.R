#' Construct a new McATAC object
#'
#' @param mat a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#' @description McATAC is a shallow object holding ATAC data over metacells.
#' Minimally it should include a count matrix of peaks over metacells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @export
McATAC <- function(mat,
                   peaks,
                   genome,
                   metadata = NULL) {
    return(make_atac_object(mat, peaks, genome, metadata, "metacell", "McATAC"))
}

#' @export
#' @noRd
print.McATAC <- function(x, ...) {
    cli::cli_text("An McATAC object with {.val {ncol(x$mat)}} metacells {.val {nrow(x$mat)}} ATAC peaks.")
    cli::cli_text("Slots include:")
    cli_ul(c("{.code $mat}: a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix."))
    cli_ul(c("{.code $peaks}: a misha intervals set with the peak definitions."))
    if (!is.null(x$metadata)) {
        cli_ul(c("{.code $metadata}: a tibble with a column called 'metacell' and additional metacell annotations."))
    }
}

#' Construct a new ScATAC object
#'
#' @param mat a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix.
#' @param metadata data frame with a column called 'cell_id' and additional per-cell annotations, or the name of a delimited file which contains such annotations.
#' @description ScATAC is a shallow object holding ATAC data over cells.
#' Minimally it should include a count matrix of peaks over cells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @inheritParams McATAC
#' @export
ScATAC <- function(mat,
                   peaks,
                   genome,
                   metadata = NULL) {
    return(make_atac_object(mat, peaks, genome, metadata, "cell_id", "ScATAC"))
}

#' @export
#' @noRd
print.ScATAC <- function(x, ...) {
    cli::cli_text("An ScATAC object with {.val {ncol(x$mat)}} cells {.val {nrow(x$mat)}} ATAC peaks from {.field {x$genome}}.")
    cli::cli_text("Slots include:")
    cli_ul(c("{.code $mat}: a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix."))
    cli_ul(c("{.code $peaks}: a misha intervals set with the peak definitions."))
    cli_ul(c("{.code $genome}: genome assembly of the peaks"))
    if (!is.null(x$metadata)) {
        cli_ul(c("{.code $metadata}: a tibble with a column called 'cell_id' and additional cell annotations."))
    }
}

make_atac_object <- function(mat, peaks, genome, metadata, metadata_id_field, class_name) {
    peaks <- PeakIntervals(peaks)
    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
    }
    rownames(mat) <- peak_names(peaks)

    if (!is.null(metadata)) {
        if (is.character(metadata)) {
            metadata <- tgutil::fread(metadata)
        }
        metadata <- as_tibble(metadata)
    }

    obj <- list(
        mat = mat,
        peaks = peaks,
        genome = genome,
        metadata = metadata
    )

    class(obj) <- class_name
    validate_atac_object(obj)
    return(obj)
}

#' Validate a McATAC or ScATAC object
#'
#' @param obj an McATAC or ScATAC object
#'
#' @noRd
validate_atac_object <- function(obj) {
    assert_that(class(obj) %in% c("McATAC", "ScATAC"))
    if (class(obj) == "McATAC") {
        validate_atac_object_params(obj$mat, obj$peaks, obj$genome, obj$metadata, "metacell")
    } else {
        validate_atac_object_params(obj$mat, obj$peaks, obj$genome, obj$metadata, "cell_id")
    }
}



validate_atac_object_params <- function(mat, peaks, genome, metadata, metadata_id_field) {
    if (is.null(genome)) {
        cli_abort("{.field genome} is missing")
    }

    if (!is.character(genome)) {
        cli_abort("{.field genome} should be a string")
    }

    if (!misha.ext::genome_exists(genome)) {
        cli_alert_warning("genome {.field {genome}} doesn't exist in your {.field misha.ext} configuration file. Note that a few functions may not work unless you add it. Your file is currently at {.file {misha.ext::find_params_yaml()}}")
    }

    if (!is.matrix(mat) && !is_sparse_matrix(mat)) {
        cli_abort("{.field mat} shuold be a matrix or a sparse matrix")
    }

    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
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
