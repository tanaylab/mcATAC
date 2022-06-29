#' Create an ScATACPeaks,McATACPeaks object from matrix and peaks
#'
#' @param mat a matrix/sparse matrix of counts
#' @param peaks an intervals set where each row corresponds to a row in \code{mat}. Can contain
#' an additional column named 'peak_name' with peak names, in which case the rownames of \code{mat}
#' are expected to equal to this column.
#' @param class should the output object be a ScATACPeaks or McATACPeaks
#' @inheritParams import_from_10x
#' @inheritParams project_atac_on_mc
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_matrix(mat, peaks, genome = "hg38")
#' atac_mc <- import_from_matrix(mat, peaks, genome = "hg38", class = "McATACPeaks")
#' }
#'
#' @export
import_from_matrix <- function(mat, peaks, genome, id = NULL, description = NULL, metadata = NULL, class = "ScATACPeaks", rm_zero_peaks = TRUE, tad_based = TRUE, zero_based = FALSE, mc_size_eps_q = 0.1) {
    if (nrow(peaks) != nrow(mat)) {
        cli_abort("Number of rows in peaks and matrix do not match")
    }
    atac_mat <- mat

    if (rm_zero_peaks) {
        non_zero_peaks <- Matrix::rowSums(atac_mat) > 0
        if (sum(!non_zero_peaks) > 0) {
            cli_alert_info("Removed {.val {sum(!non_zero_peaks)}} all-zero peaks")
        }
    } else {
        non_zero_peaks <- rep(TRUE, nrow(atac_mat))
    }

    peaks <- peaks[non_zero_peaks, ]
    atac_mat <- atac_mat[non_zero_peaks, ]

    if (!zero_based) {
        peaks <- peaks %>%
            mutate(start = start - 1, end = end - 1)
    }

    if (is.null(rownames(atac_mat)) || !has_name(peaks, "peak_name")) {
        peaks <- peaks %>%
            mutate(peak_name = misha.ext::convert_misha_intervals_to_10x_peak_names(.))
        rownames(atac_mat) <- peaks$peak_name
    }

    if (class == "ScATACPeaks") {
        res <- new("ScATACPeaks", atac_mat, peaks, genome = genome, id = id, description = description, metadata = metadata, tad_based = tad_based)
        cli_alert_success("successfully imported to an ScATACPeaks object with {.val {ncol(res@mat)}} cells and {.val {nrow(res@mat)}} ATAC peaks")
    } else if (class == "McATACPeaks") {
        res <- new("McATACPeaks", mat = atac_mat, peaks = peaks, genome = genome, metadata = metadata, cell_to_metacell = NULL, mc_size_eps_q = mc_size_eps_q, id = id, description = description, tad_based = tad_based)
        cli_alert_success("Created a new McATACPeaks object with {.val {ncol(res@mat)}} metacells and {.val {nrow(res@mat)}} ATAC peaks.")
    } else {
        cli_abort("{.field class} should be either ScATACPeaks or McATACPeaks")
    }

    return(res)
}

#' Create an ScATACPeaks,McATACPeaks object from an h5ad file
#'
#' @param file name of an h5ad file with ATAC data
#' @param class is the file storing ScATACPeaks or McATACPeaks. If NULL - the class would be determined by the 'class' field in the 'uns' part
#' of the h5ad file, if exists and otherwise the class would be McATACPeaks.
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10". If NULL - the assembly would be determined by the 'genome' field in the 'uns' part of the h5ad file.
#' @param id an identifier for the object, e.g. "pbmc". If NULL - the id would be determined by the 'id' field in the 'uns' part of the
#' h5ad file, and if this doesn't exist - a random id would be assigned.
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)". If NULL - the id would be determined by the 'description' field in the 'uns' part of the h5ad file
#'
#' @param tad_based whether to name peaks based on TADs if possible (when an intervals set named "intervs.global.tad_names" exists). When FALSE - peaks are named based on their coordinates. if NULL - the value would be determined by the 'tad_based' field in the 'uns' part of the h5ad file, and if this doesn't exist - the value would be TRUE.
#'
#' @description Reads an ATAC object from an h5ad file. Peak data is taken from the 'X' section and metadata is taken from 'obs'.
#' The 'var' section can contain a special field called 'ignore' which marks peaks that should be ignored.
#'
#' @return an ScATACPeaks/McATACPeaks object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' export_to_h5ad(atac_sc, "pbmc_data/atac_sc.h5ad")
#' atac_sc_loaded <- import_from_h5ad("pbmc_data/atac_sc.h5ad")
#' }
#'
#' @export
import_from_h5ad <- function(file, class = NULL, genome = NULL, id = NULL, description = NULL, tad_based = NULL) {
    check_files_exist(file)
    cli_li("Reading {.file {file}}")
    adata <- anndata::read_h5ad(file)

    mat <- t(adata$X)

    peaks <- adata$var
    peaks <- peaks %>%
        mutate(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))
    rownames(peaks) <- NULL
    peaks <- as_tibble(peaks)

    metadata <- adata$obs
    metadata <- metadata %>% mutate(across(where(is.factor), as.character))
    if (ncol(metadata) == 0) {
        metadata <- NULL
    }

    if (!is.null(metadata) && !is.null(rownames(metadata))) {
        colnames(mat) <- rownames(metadata)
    }

    if (is.null(class)) {
        if (!is.null(adata$uns[["class"]])) {
            class <- adata$uns[["class"]]
            if (class %!in% c("McATACPeaks", "ScATACPeaks")) {
                cli_alert_warning("Unknown class name - creating an McATACPeaks object")
                class <- "McATACPeaks"
            }
        } else {
            class <- "McATACPeaks"
        }
    }

    if (is.null(genome)) {
        if (!is.null(adata$uns[["genome"]])) {
            genome <- adata$uns[["genome"]]
        } else {
            cli_abort("h5ad file doesn't have the {.field genome} field at the {.file uns} section. Please provide the genome assembly explicitly using the {.field genome} paramter (e.g. {.code genome = \"mm10\"})")
        }
    }

    if (is.null(id) && !is.null(adata$uns[["id"]])) {
        id <- adata$uns[["id"]]
    }

    if (is.null(tad_based) && !is.null(adata$uns[["tad_based"]])) {
        tad_based <- adata$uns[["tad_based"]]
    } else {
        tad_based <- TRUE
    }

    if (is.null(description) && !is.null(adata$uns[["description"]])) {
        description <- adata$uns[["description"]]
    }

    if (class == "McATACPeaks") {
        if (!is.null(adata$uns[["mc_size_eps_q"]])) {
            mc_size_eps_q <- adata$uns[["mc_size_eps_q"]]
        } else {
            mc_size_eps_q <- 0.1
            cli_alert_warning("h5ad file doesn't have the {.field mc_size_eps_q} at the {.file uns} section. Using the default: {.val {mc_size_eps_q}}")
        }

        if (!is.null(adata$uns[["cell_to_metacell"]])) {
            cell_to_metacell <- adata$uns[["cell_to_metacell"]]
        } else {
            cell_to_metacell <- NULL
        }

        res <- new("McATACPeaks", mat = mat, peaks = peaks, genome = genome, id = id, description = description, metadata = metadata, cell_to_metacell = cell_to_metacell, mc_size_eps_q = mc_size_eps_q, path = file, tad_based = tad_based)

        if (!is.null(adata$uns[["rna_egc"]]) && !is.null(adata$uns[["rna_mcs"]]) && !is.null(adata$uns[["rna_gene_names"]])) {
            rna_egc <- adata$uns[["rna_egc"]]
            colnames(rna_egc) <- adata$uns[["rna_mcs"]]
            rownames(rna_egc) <- adata$uns[["rna_gene_names"]]
            res <- add_mc_rna(res, rna_egc)
        }
    } else {
        res <- new("ScATACPeaks", mat, peaks, genome, id, description, metadata, path = file, tad_based = tad_based)
    }

    if (has_name(peaks, "ignore")) {
        peaks_to_remove <- peaks$peak_name[peaks$ignore]
        res <- atac_ignore_peaks(res, peaks_to_remove)
        res@peaks <- res@peaks %>% select(-ignore)
        res@ignore_peaks <- res@ignore_peaks %>% select(-ignore)
    }

    if (class == "McATACPeaks") {
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res@mat)}} metacells and {.val {nrow(res@mat)}} ATAC peaks")
    } else {
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res@mat)}} cells and {.val {nrow(res@mat)}} ATAC peaks")
    }

    return(res)
}

#' Create an ScATACPeaks object from 10x directory
#'
#'
#' @param dir 10x directory. Should have the following files: 'matrix.mtx', 'barcodes.tsv' and 'features.tsv'.
#' If you used 'Cell Ranger ARC', this should be the 'outs/filtered_feature_bc_matrix' directory. \cr
#' For more details see: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/overview/welcome
#' @param matrix_fn if \code{dir} is missing, the filename of the matrix to import ("matrix.mtx" or "matrix.mtx.gz")
#' @param cells_fn if \code{dir} is missing, the filename of the cells to import ("barcodes.tsv" or "barcodes.tsv.gz")
#' @param features_fn if \code{dir} is missing, the filename of the features to import ("features.tsv" or "features.tsv.gz")
#' @param tad_based whether to name peaks based on TADs if possible (when an intervals set named "intervs.global.tad_names" exists). When FALSE - peaks are named based on their coordinates.
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param id an identifier for the object, e.g. "pbmc". If NULL - a random id would be assigned.
#' @param description (Optional) description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)".
#' @param metadata (Optional) data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#' @param zero_based are the coordinates of the features 0-based? Note that when this is FALSE (default) the coordinates of the created object would be different from the ones in \code{features_fn} by 1 bp.
#'
#'
#' @return an ScATACPeaks object
#'
#' @examples
#' \dontrun{
#' atac <- import_from_10x("./pbmc_data", genome = "hg38")
#' atac
#' }
#'
#' @export
import_from_10x <- function(dir, genome, id = NULL, description = NULL, metadata = NULL, matrix_fn = file.path(dir, "matrix.mtx"), cells_fn = file.path(dir, "barcodes.tsv"), features_fn = file.path(dir, "features.tsv"), tad_based = TRUE, zero_based = FALSE) {
    if (!file.exists(matrix_fn)) {
        matrix_fn <- glue("{matrix_fn}.gz")
    }
    if (!file.exists(cells_fn)) {
        cells_fn <- glue("{cells_fn}.gz")
    }
    if (!file.exists(features_fn)) {
        features_fn <- glue("{features_fn}.gz")
    }
    check_files_exist(c(matrix_fn, cells_fn, features_fn))

    cli_li("Importing matrix")
    mat <- tgutil::fread_mm(matrix_fn, row.names = features_fn, col.names = cells_fn)
    cli_alert_info("Imported a matrix of {.val {ncol(mat)}} cells and {.val {nrow(mat)}} features")

    cli_li("Importing features")
    features <- tgutil::fread(features_fn, col.names = c("name", "name2", "type", "chrom", "start", "end")) %>%
        as_tibble()
    stopifnot(all(features$name == rownames(mat)))

    atac_peaks <- features %>%
        filter(type == "Peaks") %>%
        select(chrom, start, end, peak_name = name)

    atac_mat <- mat[atac_peaks$peak_name, ]
    non_zero_inds <- which(Matrix::rowSums(atac_mat) > 0)
    atac_peaks <- atac_peaks[non_zero_inds, ]
    n_all_zero <- nrow(atac_mat) - length(non_zero_inds)
    atac_mat <- atac_mat[non_zero_inds, ]

    if (n_all_zero > 0) {
        cli_alert_info("Removed {.val {nrow(atac_mat)}} ATAC peaks which were all zero")
    }

    if (!zero_based) {
        atac_peaks <- atac_peaks %>%
            mutate(start = start - 1, end = end - 1) %>%
            mutate(peak_name = misha.ext::convert_misha_intervals_to_10x_peak_names(.))
        rownames(atac_mat) <- atac_peaks$peak_name
    }

    cli_alert_info("{.val {nrow(atac_mat)}} ATAC peaks")
    res <- new("ScATACPeaks", atac_mat, atac_peaks, genome = genome, id = id, description = description, metadata = metadata, path = matrix_fn, tad_based = tad_based)
    cli_alert_success("successfully imported to an ScATACPeaks object with {.val {ncol(res@mat)}} cells and {.val {nrow(res@mat)}} ATAC peaks")

    return(res)
}
