setClassUnion("any_matrix", c("sparseMatrix", "matrix"))
setClassUnion("data.frame_or_null", c("data.frame", "NULL"))

setOldClass("PeakIntervals")

#' ATAC objects
#'
#'
#' @description ATAC is a shallow object holding ATAC data over cells/metacells. Minimally it should include a count matrix of peaks over cells/metacells, \code{PeakIntervals} which hold the coordinates of the peaks and the id of the genome assembly of the peaks. ScATAC and
#' McATAC extend the ATAC object by adding metadata and additional slots.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot mat a numeric matrix where rows are peaks and columns are cells/metacells. Can be a sparse matrix.
#' @slot peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations for McATAC, or 'cell_id' and per-cell annotatoins for ScATAC. The constructor can also include or the name of a delimited file which contains such annotations.
#' @slot path original path from which the object was loaded (optional)
#'
#' @exportClass ATAC
ATAC <- setClass(
    "ATAC",
    slots = c(
        id = "character",
        description = "character",
        mat = "any_matrix",
        peaks = "PeakIntervals",
        genome = "character",
        metadata = "data.frame_or_null",
        ignore_peaks = "data.frame",
        ignore_pmat = "dgCMatrix",
        path = "character"
    ),
    contains = "VIRTUAL"
)

setMethod(
    "initialize",
    signature = "ATAC",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, path = "") {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id, description, path = path)
        validate_atac_object(.Object)
        return(.Object)
    }
)

make_atac_object <- function(obj, mat, peaks, genome, id, description, path, metadata, metadata_id_field) {
    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
    }
    rownames(mat) <- peak_names(peaks)

    peaks <- PeakIntervals(peaks, genome)
    mat <- mat[peak_names(peaks), ] # remove from matrix peaks that were filtered

    if (is.null(id)) {
        id <- gsub(" ", "_", randomNames::randomNames(n = 1, name.order = "first.last", name.sep = "_"))
        cli_alert("No id was given, setting id to {.field {id}}")
    }

    description <- description %||% ""

    obj@id <- id
    obj@description <- description
    obj@path <- path
    obj@mat <- mat
    obj@peaks <- peaks
    obj@genome <- genome
    obj@ignore_peaks <- subset(peaks, subset = rep(FALSE, nrow(peaks)))
    obj@ignore_pmat <- methods::as(matrix(0, nrow = 0, ncol = ncol(obj@mat)), "dgCMatrix")
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
#' @slot egc normalized metacell accessibility (fraction of accessibility per metacell scaled to the \code{mc_size_eps_q} quantile of
#' metacell size)
#' @slot fp a matrix showing for each peak (row) the relative enrichment of umis in log2 scale, i.e. \eqn{log2((1 + egc) / median(1 + egc))}
#'
#' @rdname ATAC
#' @exportClass McATAC
McATAC <- setClass(
    "McATAC",
    slots = c(
        egc = "any_matrix",
        fp = "any_matrix",
        mc_size_eps_q = "numeric"
    ),
    contains = "ATAC"
)


#' Construct a new McATAC object
#'
#' @param mat a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k),
#' projection was done using RNA metacells"
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#' @param path path from which the object was loaded (optional)
#'
#' @description McATAC is a shallow object holding ATAC data over metacells.
#' Minimally it should include a count matrix of peaks over metacells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @inheritParams project_atac_on_mc
#' @export
setMethod(
    "initialize",
    signature = "McATAC",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, metadata = NULL, mc_size_eps_q = 0.1, path = "") {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id = id, description = description, path = path)
        validate_atac_object(.Object)
        .Object <- add_metadata(.Object, metadata, "metacell")
        .Object@egc <- calc_mc_egc(.Object, mc_size_eps_q)
        .Object@fp <- calc_mc_fp(.Object)
        .Object@mc_size_eps_q <- mc_size_eps_q
        return(.Object)
    }
)

calc_mc_egc <- function(mcatac, mc_size_eps_q = 0.1) {
    mc_mat <- mcatac@mat
    mc_sum <- colSums(mc_mat, na.rm = TRUE)
    fractions <- t(t(mc_mat) / mc_sum)
    quant_size <- quantile(mc_sum, mc_size_eps_q, na.rm = TRUE)
    cli_li("Setting {.field egc} cell size to {.val {quant_size}} (the {.val {mc_size_eps_q}} quantile of metacell sizes)")
    egc <- fractions * quant_size
    return(egc)
}

calc_mc_fp <- function(mcatac) {
    log_egc <- log2(1 + mcatac@egc)
    mc_fp <- log_egc - sparseMatrixStats::rowMedians(log_egc, na.rm = TRUE)
    return(mc_fp)
}

#' Add per-metacell metadata to an McATAC object
#'
#' @param mcatac an McATAC object
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#'
#' @examples
#' \dontrun{
#' data(mcmd)
#' mc_atac <- add_metadata(mc_atac, mcmd)
#' }
#'
#' @export
add_mc_metadata <- function(mcatac, metadata) {
    add_metadata(mcatac, metadata, "metacell")
}

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
        print_atac_object(object, "McATAC", "metacell", "metacell")
    }
)

print_atac_object <- function(object, object_type, column_type, md_column) {
    cli::cli_text("{.cls {object_type}} object with {.val {ncol(object@mat)}} {column_type}s and {.val {nrow(object@mat)}} ATAC peaks from {.field {object@genome}}.")
    if (object@id != "") {
        cli::cli_text(c("id: {.val {object@id}}"))
    }
    if (object@description != "") {
        cli::cli_text(c("description: {.val {object@description}}"))
    }
    if (object@path != "") {
        cli::cli_text(c("loaded from: {.file {object@path}}"))
    }
    cli::cli_text("Slots include:")
    cli_ul(c("{.code @mat}: a numeric matrix where rows are peaks and columns are {column_type}s. Can be a sparse matrix."))
    cli_ul(c("{.code @peaks}: a misha intervals set with the peak definitions."))
    cli_ul(c("{.code @genome}: genome assembly of the peaks"))
    if (object_type == "McATAC") {
        cli_ul(c("{.code @egc}: a numeric matrix which contains normalized metacell accessibility."))
        cli_ul(c("{.code @fp}: a matrix showing for each peak (row) the relative enrichment of umis in log2 scale."))
    }
    if (!is.null(object@metadata)) {
        cli_ul(c("{.code @metadata}: a tibble with a column called '{md_column}' and additional {column_type} annotations."))
    }
    if (nrow(object@ignore_peaks) > 0) {
        cli_alert("{.val {nrow(object@ignore_peaks)}} peaks are ignored. You can access them at {.code @ignore_peaks} and {.code @ignore_pmat}")
    }
}


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


#' Construct a new ScATAC object
#'
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param mat a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param metadata data frame with a column called 'cell_id' and additional per-cell annotations, or the name of a delimited file which contains such annotations.
#' @param path path from which the object was loaded (optional)
#'
#' @description ScATAC is a shallow object holding ATAC data over cells.
#' Minimally it should include a count matrix of peaks over cells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#'
#' @export
setMethod(
    "initialize",
    signature = "ScATAC",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, metadata = NULL, path = "") {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id = id, description = description, path = path)
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
        print_atac_object(object, "ScATAC", "cell", "cell_id")
    }
)


#' Set ignored (i.e. blacklisted) peaks
#'
#' Given a list of peaks to ignore, this will remove the given peaks from the {.code ignore_peaks} and {.code ignore_pmat} slots. Note that the matrix and peaks would be reordered.
#'
#' @param atac an ScATAC or McATAC object
#' @param ig_peaks a PeakIntervals object, or vector of peak names to ignore
#' @param reset rest current ignore policy if exists. When set to \code{TRUE}, the current ignore policy will be removed, otherwise the new ignore policy will be appended to the current one.
#'
#' @return
#' @examples
#' \dontrun{
#'
#' }
#' @export
atac_ignore_peaks <- function(atac, ig_peaks, reset = FALSE) {
    assert_atac_object(atac)

    if (length(ig_peaks) == 0) {
        cli_alert_warning("Peaks to ignore should be specified (they are either NULL or length 0), returning original object.")
        return(atac)
    }

    if (is.null(dim(ig_peaks)) && length(ig_peaks) > 0 && is.character(ig_peaks)) {
        ig_peaks <- misha.ext::convert_10x_peak_names_to_misha_intervals(ig_peaks)
    }

    ig_peaks <- ig_peaks %>% distinct(chrom, start, end, peak_name)
    n_cur_ig <- nrow(ig_peaks)
    n_prev_ig <- nrow(atac@ignore_peaks)
    if (n_prev_ig > 0) {
        if (!reset) {
            ig_peaks <- bind_rows(atac@ignore_peaks, ig_peaks) %>%
                distinct(chrom, start, end, peak_name)
            cli_alert_warning("Adding to previous ignore policy ({.val {n_prev_ig}} peaks).")
        } else {
            cli_alert_warning("Previous ignore policy ({.val {n_prev_ig}} peaks) is being reset.")
        }
    }

    atac@mat <- rbind(atac@mat, atac@ignore_pmat)
    peaks_merge <- bind_rows(atac@peaks, atac@ignore_peaks)
    cn <- c("chrom", "start", "end", "peak_name")
    new_ord <- with(peaks_merge, order(chrom, start))
    atac@mat <- atac@mat[new_ord, ]
    atac@peaks <- peaks_merge[new_ord, ]
    atac@peaks$temp_intID <- 1:nrow(atac@peaks)
    good_peaks <- anti_join(atac@peaks, ig_peaks, by = cn)
    atac@ignore_peaks <- semi_join(atac@peaks, ig_peaks, by = cn)
    atac@ignore_pmat <- atac@mat[atac@peaks$temp_intID %in% atac@ignore_peaks$temp_intID, ]
    atac@mat <- atac@mat[atac@peaks$temp_intID %in% good_peaks$temp_intID, ]
    atac@peaks <- good_peaks %>%
        select(any_of(cn), everything()) %>%
        select(-temp_intID)
    atac@ignore_peaks <- atac@ignore_peaks %>% select(-temp_intID)

    n_removed_peaks <- nrow(ig_peaks)
    n_good_peaks <- nrow(atac@peaks)
    n_tot_peaks <- n_removed_peaks + n_good_peaks

    if (n_prev_ig == 0) {
        cli_alert_success("Removed {.val {n_removed_peaks}} peaks out of {.val {n_tot_peaks}} {.field ({scales::percent(n_removed_peaks/n_tot_peaks)})}. The object is left with {.val {n_good_peaks}} peaks.")
    } else {
        cli_alert_success("Removed {.val {n_cur_ig}} peaks out of {.val {n_tot_peaks}} {.field ({scales::percent(n_cur_ig/n_tot_peaks)})}. The object is left with {.val {n_good_peaks}} peaks {.field ({scales::percent(n_removed_peaks/n_tot_peaks)})}.")
    }

    return(atac)
}
