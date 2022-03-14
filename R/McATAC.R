

setClassUnion("any_matrix", c("sparseMatrix", "matrix"))
setClassUnion("data.frame_or_null", c("data.frame", "NULL"))

setOldClass("PeakIntervals")

#' ATAC objects
#'
#'
#' @description ATAC is a shallow object holding ATAC data over cells/metacells. Minimally it should include a count matrix of peaks over cells/metacells, \code{PeakIntervals} which hold the coordinates of the peaks and the id of the genome assembly of the peaks. ScATAC and
#' McATAC extend the ATAC object by adding metadata and additional slots.
#'
#' @slot mat a numeric matrix where rows are peaks and columns are cells/metacells. Can be a sparse matrix.
#' @slot peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations for McATAC, or 'cell_id' and per-cell annotatoins for ScATAC. The constructor can also include or the name of a delimited file which contains such annotations.
#'
#' @exportClass ATAC
ATAC <- setClass(
    "ATAC",
    slots = c(
        mat = "any_matrix",
        peaks = "PeakIntervals",
        genome = "character",
        metadata = "data.frame_or_null"
    ),
    contains = "VIRTUAL"
)

setMethod(
    "initialize",
    signature = "ATAC",
    definition = function(.Object, mat, peaks, genome) {
        .Object <- make_atac_object(.Object, mat, peaks, genome)
        validate_atac_object(.Object)
        return(.Object)
    }
)

make_atac_object <- function(obj, mat, peaks, genome, metadata, metadata_id_field) {
    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
    }
    rownames(mat) <- peak_names(peaks)

    peaks <- PeakIntervals(peaks, genome)
    mat <- mat[peak_names(peaks), ] # remove from matrix peaks that were filtered

    obj@mat <- mat
    obj@peaks <- peaks
    obj@genome <- genome

    validate_atac_object(obj)
    return(obj)
}

#' Validate a McATAC or ScATAC object
#'
#' @param obj an McATAC or ScATAC object
#'
#' @noRd
validate_atac_object <- function(obj) {
    validate_atac_object_params(obj@mat, obj@peaks, obj@genome)
}



validate_atac_object_params <- function(mat, peaks, genome) {
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
}



#' McATAC
#'
#' An ATAC object with data over metacells
#'
#' @rdname ATAC
#' @exportClass McATAC
McATAC <- setClass(
    "McATAC",
    contains = "ATAC"
)


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
setMethod(
    "initialize",
    signature = "McATAC",
    definition = function(.Object, mat, peaks, genome, metadata = NULL) {
        .Object <- make_atac_object(.Object, mat, peaks, genome)
        validate_atac_object(.Object)
        .Object <- add_metadata(.Object, metadata, "metacell")
        return(.Object)
    }
)

add_metadata <- function(obj, metadata, metadata_id_field) {
    if (!is.null(metadata)) {
        if (is.character(metadata)) {
            metadata <- tgutil::fread(metadata)
        }
        if (!is.data.frame(metadata)) {
            cli_abort("{.field metadata} is not a data frame")
        }

        metadata <- as_tibble(metadata)

        if (!has_name(metadata, metadata_id_field)) {
            cli_abort("{.field metadata} doesn't have the required field {.field {metadata_id_field}}")
        }

        metadata <- as_tibble(metadata)

        # make sure that all cells/metacells exist within the matrix
        missing_cells <- metadata[[metadata_id_field]] %!in% colnames(obj@mat)
        if (any(missing_cells)) {
            missing_cells <- paste(unique(metadata[[metadata_id_field]][missing_cells]), collapse = ", ")
            cli_abort("The following {metadata_id_field}s are missing from {.field mat} colnames: {.val {missing_cells}}")
        }
    }

    obj@metadata <- metadata

    return(obj)
}

#' @export
#' @noRd
setMethod(
    "show",
    signature = "McATAC",
    definition = function(object) {
        cli::cli_text("An McATAC object with {.val {ncol(object@mat)}} metacells {.val {nrow(object@mat)}} ATAC peaks.")
        cli::cli_text("Slots include:")
        cli_ul(c("{.code @mat}: a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix."))
        cli_ul(c("{.code @peaks}: a misha intervals set with the peak definitions."))
        cli_ul(c("{.code @genome}: genome assembly of the peaks"))
        if (!is.null(object@metadata)) {
            cli_ul(c("{.code @metadata}: a tibble with a column called 'metacell' and additional metacell annotations."))
        }
    }
)


#' ScATAC
#'
#' An ATAC object with data over cells
#'
#' @rdname ATAC
#' @exportClass ScATAC
ScATAC <- setClass(
    "ScATAC",
    contains = "ATAC"
)


#' Construct a new McATAC object
#'
#' @param mat a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param metadata data frame with a column called 'cell_id' and additional per-cell annotations, or the name of a delimited file which contains such annotations.
#' @description ScATAC is a shallow object holding ATAC data over cells.
#' Minimally it should include a count matrix of peaks over cells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#'
#' @export
setMethod(
    "initialize",
    signature = "ScATAC",
    definition = function(.Object, mat, peaks, genome, metadata = NULL) {
        .Object <- make_atac_object(.Object, mat, peaks, genome)
        validate_atac_object(.Object)
        .Object <- add_metadata(.Object, metadata, "cell_id")
        return(.Object)
    }
)

#' @export
#' @noRd
setMethod(
    "show",
    signature = "ScATAC",
    definition = function(object) {
        cli::cli_text("An ScATAC object with {.val {ncol(object@mat)}} cells {.val {nrow(object@mat)}} ATAC peaks from {.field {object@genome}}.")
        cli::cli_text("Slots include:")
        cli_ul(c("{.code @mat}: a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix."))
        cli_ul(c("{.code @peaks}: a misha intervals set with the peak definitions."))
        cli_ul(c("{.code @genome}: genome assembly of the peaks"))
        if (!is.null(object@metadata)) {
            cli_ul(c("{.code @metadata}: a tibble with a column called 'cell_id' and additional cell annotations."))
        }
    }
)
