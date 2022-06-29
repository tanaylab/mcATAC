setClassUnion("any_matrix", c("sparseMatrix", "matrix"))
setClassUnion("data.frame_or_null", c("data.frame", "NULL"))

setOldClass("PeakIntervals")

#' ATACPeaks objects
#'
#'
#' @description ATACPeaks is a shallow object holding ATAC data over cells/metacells. Minimally it should include a count matrix of peaks over cells/metacells, \code{PeakIntervals} which hold the coordinates of the peaks and the id of the genome assembly of the peaks. ScATACPeaks and
#' McATACPeaks extend the ATACPeaks object by adding metadata and additional slots.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot mat a numeric matrix where rows are peaks and columns are cells/metacells. Can be a sparse matrix.
#' @slot peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot metadata data frame with a column called 'metacell' and additional metacell annotations for McATACPeaks, or 'cell_id' and per-cell annotatoins for ScATACPeaks. The constructor can also include or the name of a delimited file which contains such annotations.
#' @slot path original path from which the object was loaded (optional)
#' @slot promoters are the peaks promoters (optional). When peaks are promoters the peak name would be the gene name instead of coordinates.
#' @slot tad_based are the peaks named based on TADs.
#'
#' @exportClass ATACPeaks
ATACPeaks <- setClass(
    "ATACPeaks",
    slots = c(
        id = "character",
        description = "character",
        mat = "any_matrix",
        peaks = "PeakIntervals",
        genome = "character",
        metadata = "data.frame_or_null",
        ignore_peaks = "PeakIntervals",
        ignore_pmat = "dgCMatrix",
        path = "character",
        promoters = "logical",
        tad_based = "logical"
    ),
    contains = "VIRTUAL"
)

setMethod(
    "initialize",
    signature = "ATACPeaks",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, path = "", tad_based = TRUE, rename = TRUE) {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id, description, path = path, tad_based = tad_based, rename = rename)
        validate_atac_object(.Object)
        return(.Object)
    }
)

make_atac_object <- function(obj, mat, peaks, genome, id, description, path, metadata, metadata_id_field, tad_based, rename = TRUE) {
    if (nrow(mat) != nrow(peaks)) {
        cli_abort("Number of peaks is not equal to the matrix rows.")
    }

    peaks <- PeakIntervals(peaks, genome)
    if (!is.null(rownames(mat)) && has_name(peaks, "peak_name")) {
        mat <- mat[peaks$peak_name, ] # filter out peaks that do not exists in peak intervals
    }
    gset_genome(genome)
    if (rename || !has_name(peaks, "peak_name")) {
        peaks <- peaks %>% select(-any_of("peak_name"))
        peaks$peak_name <- peak_names(peaks, tad_based = tad_based)
    }
    rownames(mat) <- peaks$peak_name

    if (is.null(id)) {
        id <- gsub(" ", "_", randomNames::randomNames(n = 1, name.order = "first.last", name.sep = "_"))
        cli_alert("No id was given, setting id to {.field {id}}")
    }

    description <- description %||% ""

    if (path != "") {
        path <- as.character(normalizePath(path))
    }

    obj@id <- id
    obj@description <- description
    obj@path <- path
    obj@mat <- mat
    obj@peaks <- peaks
    obj@genome <- genome
    obj@ignore_peaks <- subset(peaks, subset = rep(FALSE, nrow(peaks)))
    obj@ignore_pmat <- methods::as(matrix(0, nrow = 0, ncol = ncol(obj@mat)), "dgCMatrix")
    obj@promoters <- FALSE
    obj@tad_based <- tad_based
    validate_atac_object(obj)
    return(obj)
}

#' Validate a McATACPeaks or ScATACPeaks object
#'
#' @param obj an McATACPeaks or ScATACPeaks object
#'
#' @noRd
validate_atac_object <- function(obj) {
    validate_atac_object_params(obj@mat, obj@peaks, obj@genome, obj@promoters)
}



validate_atac_object_params <- function(mat, peaks, genome, promoters = FALSE) {
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
    if (any(rownames(mat) != peak_names(peaks, tad_based = TRUE, promoters = promoters))) {
        cli_abort("rownames of the matrix are not the same as the peak names.")
    }
}



#' McATACPeaks
#'
#' An ATACPeaks object with data over metacells
#'
#' @slot egc normalized metacell accessibility: fraction of accessibility per metacell scaled to the \code{mc_size_eps_q} quantile of
#' metacell size. Accessibility is normalized by peak length.
#' @slot fp a matrix showing for each peak (row) the relative enrichment of umis in log2 scale, i.e. \eqn{log2((1 + egc) / median(1 + egc))}
#' @slot mc_size_eps_q quantile of MC size (in UMIs) to scale the number of UMIs per metacell. See \code{project_atac_on_mc}
#' @slot rna_egc normalized gene expression per gene per metacell (optional). Can be created using \code{add_mc_rna}
#'
#' @rdname ATACPeaks
#' @exportClass McATACPeaks
McATACPeaks <- setClass(
    "McATACPeaks",
    slots = c(
        egc = "any_matrix",
        fp = "any_matrix",
        mc_size_eps_q = "numeric",
        cell_to_metacell = "data.frame_or_null",
        rna_egc = "any_matrix"
    ),
    contains = "ATACPeaks"
)


#' Construct a new McATACPeaks object
#'
#' @param mat a numeric matrix where rows are peaks and columns are metacells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k),
#' projection was done using RNA metacells"
#' @param cell_to_metacell a data frame mapping 'cell_id' to 'metacell' (optional). See \code{project_atac_on_mc}
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#' @param path path from which the object was loaded (optional)
#' @param tad_based whether the peak names are TAD-based
#'
#' @description McATACPeaks is a shallow object holding ATAC data over metacells.
#' Minimally it should include a count matrix of peaks over metacells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#' @inheritParams project_atac_on_mc
#' @export
setMethod(
    "initialize",
    signature = "McATACPeaks",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, metadata = NULL, cell_to_metacell = NULL, mc_size_eps_q = 0.1, path = "", tad_based = TRUE) {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id = id, description = description, path = path, tad_based = tad_based, rename = FALSE)
        validate_atac_object(.Object)
        .Object <- add_metadata(.Object, metadata, "metacell")
        .Object@egc <- calc_mc_egc(.Object, mc_size_eps_q)
        .Object@fp <- calc_mc_fp(.Object)
        .Object@mc_size_eps_q <- mc_size_eps_q
        .Object@rna_egc <- matrix(0, nrow = 0, ncol = ncol(.Object@mat), dimnames = list(NULL, colnames(.Object@mat)))
        .Object@cell_to_metacell <- cell_to_metacell
        return(.Object)
    }
)

calc_mc_egc <- function(mcatac, mc_size_eps_q = 0.1) {
    mc_mat <- mcatac@mat
    peak_len <- mcatac@peaks$end - mcatac@peaks$start
    mc_mat <- mc_mat / peak_len
    mc_sum <- colSums(mc_mat, na.rm = TRUE)
    fractions <- t(t(mc_mat) / mc_sum)
    quant_size <- quantile(colSums(mcatac@mat, na.rm = TRUE), mc_size_eps_q, na.rm = TRUE)
    cli_li("Setting {.field egc} cell size to {.val {quant_size}} (the {.val {mc_size_eps_q}} quantile of metacell sizes)")
    egc <- fractions * quant_size
    return(as.matrix(egc))
}

calc_mc_fp <- function(mcatac) {
    log_egc <- as.matrix(log2(1 + mcatac@egc))
    mc_fp <- log_egc - sparseMatrixStats::rowMedians(log_egc, na.rm = TRUE)
    return(mc_fp)
}

#' @export
#' @noRd
setMethod(
    "show",
    signature = "McATACPeaks",
    definition = function(object) {
        print_atac_object(object, "McATACPeaks", "metacell", "metacell")
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
        cli::cli_text(c("Loaded from: {.file {object@path}}"))
    }
    if (object@promoters) {
        cli::cli_text(c("A promoter object (peaks are on promoters)"))
    }
    cli::cli_text("Slots include:")
    cli_ul(c("{.code @mat}: a numeric matrix where rows are peaks and columns are {column_type}s. Can be a sparse matrix."))
    cli_ul(c("{.code @peaks}: a misha intervals set with the peak definitions."))
    cli_ul(c("{.code @genome}: genome assembly of the peaks"))
    if (object_type == "McATACPeaks") {
        cli_ul(c("{.code @egc}: a numeric matrix which contains normalized metacell accessibility."))
        cli_ul(c("{.code @fp}: a matrix showing for each peak (row) the relative enrichment of umis in log2 scale."))
        if (has_rna(object)) {
            cli_ul(c("{.code @rna_egc}: a numeric matrix which contains normalized RNA expression per gene (rows) per metacell (columns)."))
        }
    }
    if (!is.null(object@metadata)) {
        cli_ul(c("{.code @metadata}: a tibble with a column called '{md_column}' and additional {column_type} annotations."))
    }
    if (nrow(object@ignore_peaks) > 0) {
        cli_alert("{.val {nrow(object@ignore_peaks)}} peaks are ignored. You can access them at {.code @ignore_peaks} and {.code @ignore_pmat}")
    }
}


#' ScATACPeaks
#'
#' An ATACPeaks object with data over cells
#'
#' @rdname ATACPeaks
#' @exportClass ScATACPeaks
ScATACPeaks <- setClass(
    "ScATACPeaks",
    contains = "ATACPeaks"
)


#' Construct a new ScATACPeaks object
#'
#' @param id an identifier for the object, e.g. "pbmc".
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @param mat a numeric matrix where rows are peaks and columns are cells. Can be a sparse matrix.
#' @param peaks misha intervals set. Can contain a field named 'peak_name' with a unique name per peak. Both the names and intervals should be unique (a peak cannot appear more than once).
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param metadata data frame with a column called 'cell_id' and additional per-cell annotations, or the name of a delimited file which contains such annotations.
#' @param path path from which the object was loaded (optional)
#' @param tad_based wether the peak names are TAD-based
#'
#' @description ScATACPeaks is a shallow object holding ATAC data over cells.
#' Minimally it should include a count matrix of peaks over cells, and \code{PeakIntervals} which hold the coordinates
#' of the peaks.
#'
#'
#' @export
setMethod(
    "initialize",
    signature = "ScATACPeaks",
    definition = function(.Object, mat, peaks, genome, id = NULL, description = NULL, metadata = NULL, path = "", tad_based = TRUE) {
        .Object <- make_atac_object(.Object, mat, peaks, genome, id = id, description = description, path = path, tad_based = tad_based)
        validate_atac_object(.Object)
        .Object <- add_metadata(.Object, metadata, "cell_id")
        return(.Object)
    }
)

#' @export
#' @noRd
setMethod(
    "show",
    signature = "ScATACPeaks",
    definition = function(object) {
        print_atac_object(object, "ScATACPeaks", "cell", "cell_id")
    }
)


#' Set ignored (i.e. blacklisted) peaks
#'
#' Given a list of peaks to ignore, this will remove the given peaks from the {.code ignore_peaks} and {.code ignore_pmat} slots. Note that the matrix and peaks would be reordered.
#'
#' @param atac an ScATACPeaks or McATACPeaks object
#' @param ig_peaks a PeakIntervals object, or vector of peak names to ignore
#' @param reset rest the current ignore policy if exists. When set to \code{TRUE}, the current ignore policy will be removed, otherwise the new ignore policy will be appended to the current one.
#'
#' @return an ScATACPeaks or McATACPeaks object with the new ignore policy.
#' @examples
#' \dontrun{
#' max_peak_length <= 1000
#' peak_stats <- get_peak_coverage_stats(atac_sc, scale = 100)
#' too_long_peaks <- peak_stats$peak_name[peak_stats$len > max_peak_length]
#' atac_sc_filtered <- atac_ignore_peaks(atac_sc, too_long_peaks)
#' }
#' @export
atac_ignore_peaks <- function(atac, ig_peaks, reset = FALSE) {
    assert_atac_object(atac)

    if (length(ig_peaks) == 0) {
        cli_alert_warning("Peaks to ignore should be specified (they are either NULL or length 0), returning original object.")
        return(atac)
    }

    peaks_merge <- bind_rows(atac@peaks, atac@ignore_peaks)
    # if ig_peaks is a character vector, we need to convert it to a PeakIntervals object by selecting the
    # peaks with the given names
    if (is.null(dim(ig_peaks)) && length(ig_peaks) > 0 && is.character(ig_peaks)) {
        ig_peaks <- peaks_merge %>% filter(peak_name %in% ig_peaks)
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
    cn <- c("chrom", "start", "end", "peak_name")
    new_ord <- with(peaks_merge, order(chrom, start))
    atac@mat <- atac@mat[new_ord, ]
    atac@peaks <- peaks_merge[new_ord, ]
    atac@peaks$temp_intID <- 1:nrow(atac@peaks)
    good_peaks <- anti_join(atac@peaks, ig_peaks, by = cn)
    atac@ignore_peaks <- semi_join(atac@peaks, ig_peaks, by = cn)
    atac@ignore_pmat <- atac@mat[atac@peaks$temp_intID %in% atac@ignore_peaks$temp_intID, ]
    atac@mat <- atac@mat[atac@peaks$temp_intID %in% good_peaks$temp_intID, ]
    if (class(atac) == "McATACPeaks") {
        atac@egc <- atac@egc[atac@peaks$temp_intID %in% good_peaks$temp_intID, ]
        atac@fp <- atac@fp[atac@peaks$temp_intID %in% good_peaks$temp_intID, ]
    }
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

#' Remove cells
#'
#' This function removes a given list of cells from the {.code mat} slot.
#'
#' @param atac an ScATACPeaks object
#' @param ig_cells a vector of cell names (10X barcodes) to removes
#'
#' @return an scATAC object with the given cells removed.
#'
#' @examples
#' \dontrun{
#' cs <- Matrix::colSums(atac_sc@mat)
#' big_cells <- names(cs)[which(cs >= quantile(cs, 0.98))]
#' atac_sc_filtered <- atac_ignore_cells(atac_sc, big_cells)
#' }
#' @export
atac_ignore_cells <- function(atac, ig_cells) {
    assert_atac_object(atac)
    cells_in <- colnames(atac@mat)[colnames(atac@mat) %!in% ig_cells]
    if (length(ig_cells) == 0) {
        cli_alert_warning("Cells to ignore should be specified (they are either NULL or length 0), returning original object.")
        return(atac)
    }
    atac@mat <- atac@mat[, cells_in]
    if (nrow(atac@ignore_pmat) > 0) {
        atac@ignore_pmat <- atac@ignore_pmat[, cells_in]
    } else {
        atac@ignore_pmat <- methods::as(matrix(0, nrow = 0, ncol = ncol(atac@mat)), "dgCMatrix")
    }
    return(atac)
}
