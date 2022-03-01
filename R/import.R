

#' Create an ScATAC object from an h5ad file
#'
#' @param file name of an h5ad file with ATAC data
#'
#' @return an ScATAC object
#'
#' @export
import_from_h5ad <- function(file) {

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

    cli_alert("Importing matrix")
    mat <- tgutil::fread_mm(matrix_fn, row.names = features_fn, col.names = cells_fn)
    cli_alert_info("Imported a matrix of {.val {ncol(mat)}} cells and {.val {nrow(mat)}} features")

    cli_alert("Importing features")
    features <- tgutil::fread(features_fn, col.names = c("name", "name2", "type", "chrom", "start", "end")) %>%
        as_tibble()
    stopifnot(all(features$name == rownames(mat)))

    atac_peaks <- features %>%
        filter(type == "Peaks") %>%
        select(chrom, start, end, peak_name = name)

    atac_mat <- mat[atac_peaks$peak_name, ]
    cli_alert_info("{.val {nrow(atac_mat)}} ATAC peaks")

    res <- ScATAC(atac_mat, atac_peaks)
    cli_alert_success("Succesfully imported to an ScATAC object with {.val {ncol(atac_mat)}} cells {.val {nrow(atac_mat)}} ATAC peaks")

    return(res)
}
