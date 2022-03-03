

#' Create an ScATAC,McATAC object from an h5ad file
#'
#' @param file name of an h5ad file with ATAC data
#' @param class is the file storing ScATAC or McATAC. If NULL - the class would be determined by the 'class' field in the 'uns' part
#' of the h5ad file, if exists, and otherwise the class would be
#' McATAC.
#'
#' @return an ScATAC/McATAC object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data")
#' export_to_h5ad(atac_sc, "pbmc_data/atac_sc.h5ad")
#' atac_sc_loaded <- import_from_h5ad("pbmc_data/atac_sc.h5ad")
#' }
#'
#' @export
import_from_h5ad <- function(file, class = NULL) {
    check_files_exist(file)
    cli_ul("Reading {.file {file}}")
    adata <- anndata::read_h5ad(file)

    mat <- t(adata$X)

    peaks <- adata$var
    metadata <- adata$obs
    if (ncol(metadata) == 0) {
        metadata <- NULL
    }

    if (is.null(class)) {
        if (!is.null(adata$uns[["class"]])) {
            class <- adata$uns[["class"]]
            if (class %!in% c("McATAC", "ScATAC")) {
                cli_alert_warning("Unknown class name - creating an McATAC object")
                class <- "McATAC"
            }
        } else {
            class <- "McATAC"
        }
    }

    if (class == "McATAC") {
        res <- McATAC(mat, peaks, metadata)
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res$mat)}} metacells and {.val {nrow(res$mat)}} ATAC peaks")
    } else {
        res <- ScATAC(mat, peaks, metadata)
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res$mat)}} cells and {.val {nrow(res$mat)}} ATAC peaks")
    }

    return(res)
}

#' Create an ScATAC object from 10x directory
#'
#' @param dir 10x directory. Should have the following files: 'matrix.mtx', 'barcodes.tsv' and 'features.tsv'.
#' See 'Cell Ranger Arc' for more details: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/overview/welcome
#' @param matrix_fn if \code{dir} is missing, the filename of the matrix to import ("matrix.mtx")
#' @param cells_fn if \code{dir} is missing, the filename of the cells to import ("barcodes.tsv")
#' @param features_fn if \code{dir} is missing, the filename of the features to import ("features.tsv")
#'
#' @return an ScATAC object
#'
#' @examples
#' \dontrun{
#' atac <- import_from_10x("./pbmc_data")
#' atac
#' }
#'
#' @export
import_from_10x <- function(dir = NULL, matrix_fn = file.path(dir, "matrix.mtx"), cells_fn = file.path(dir, "barcodes.tsv"), features_fn = file.path(dir, "features.tsv")) {
    check_files_exist(c(matrix_fn, cells_fn, features_fn))

    cli_ul("Importing matrix")
    mat <- tgutil::fread_mm(matrix_fn, row.names = features_fn, col.names = cells_fn)
    cli_alert_info("Imported a matrix of {.val {ncol(mat)}} cells and {.val {nrow(mat)}} features")

    cli_ul("Importing features")
    features <- tgutil::fread(features_fn, col.names = c("name", "name2", "type", "chrom", "start", "end")) %>%
        as_tibble()
    stopifnot(all(features$name == rownames(mat)))

    atac_peaks <- features %>%
        filter(type == "Peaks") %>%
        select(chrom, start, end, peak_name = name)

    atac_mat <- mat[atac_peaks$peak_name, ]
    cli_alert_info("{.val {nrow(atac_mat)}} ATAC peaks")

    res <- ScATAC(atac_mat, atac_peaks)
    cli_alert_success("Succesfully imported to an ScATAC object with {.val {ncol(atac_mat)}} cells and {.val {nrow(atac_mat)}} ATAC peaks")

    return(res)
}
