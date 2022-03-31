

#' Create an ScATAC,McATAC object from an h5ad file
#'
#' @param file name of an h5ad file with ATAC data
#' @param class is the file storing ScATAC or McATAC. If NULL - the class would be determined by the 'class' field in the 'uns' part
#' of the h5ad file, if exists and otherwise the class would be McATAC.
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10". If NULL - the assembly would be determined by the 'genome' field in the 'uns' part of the h5ad file.
#' @param id an identifier for the object, e.g. "pbmc". If NULL - the id would be determined by the 'id' field in the 'uns' part of the
#' h5ad file, and if this doesn't exist - a random id would be assigned.
#' @param description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)". If NULL - the id would be determined by the 'description' field in the 'uns' part of the h5ad file
#'
#' @description Reads an ATAC object from an h5ad file. Peak data is taken from the 'X' section and metadata is taken from 'obs'.
#' The 'var' section can contain a special field called 'ignore' which marks peaks that should be ignored.
#'
#' @return an ScATAC/McATAC object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' export_to_h5ad(atac_sc, "pbmc_data/atac_sc.h5ad")
#' atac_sc_loaded <- import_from_h5ad("pbmc_data/atac_sc.h5ad")
#' }
#'
#' @export
import_from_h5ad <- function(file, class = NULL, genome = NULL, id = NULL, description = NULL) {
    check_files_exist(file)
    cli_ul("Reading {.file {file}}")
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

    if (is.null(description) && !is.null(adata$uns[["description"]])) {
        description <- adata$uns[["description"]]
    }

    if (class == "McATAC") {
        if (!is.null(adata$uns[["mc_size_eps_q"]])) {
            mc_size_eps_q <- adata$uns[["mc_size_eps_q"]]
        } else {
            mc_size_eps_q <- 0.1
            cli_alert_warning("h5ad file doesn't have the {.field mc_size_eps_q} at the {.file uns} section. Using the default: {.val {mc_size_eps_q}")
        }
        res <- new("McATAC", mat, peaks, genome, id, description, metadata, mc_size_eps_q = mc_size_eps_q, path = file)
        if (!is.null(adata$uns[["rna_egc"]]) && !is.null(adata$uns[["rna_mcs"]]) && !is.null(adata$uns[["rna_gene_names"]])) {
            rna_egc <- adata$uns[["rna_egc"]]
            colnames(rna_egc) <- adata$uns[["rna_mcs"]]
            rownames(rna_egc) <- adata$uns[["rna_gene_names"]]
            res <- add_mc_rna(res, rna_egc)
        }
    } else {
        res <- new("ScATAC", mat, peaks, genome, id, description, metadata, path = file)
    }

    if (has_name(peaks, "ignore")) {
        peaks_to_remove <- peaks$peak_name[peaks$ignore]
        res <- atac_ignore_peaks(res, peaks_to_remove)
        res@peaks <- res@peaks %>% select(-ignore)
        res@ignore_peaks <- res@ignore_peaks %>% select(-ignore)
    }

    if (class == "McATAC") {
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res@mat)}} metacells and {.val {nrow(res@mat)}} ATAC peaks")
    } else {
        cli_alert_success("Successfully loaded an {.var {class}} object with {.val {ncol(res@mat)}} cells and {.val {nrow(res@mat)}} ATAC peaks")
    }

    return(res)
}

#' Create an ScATAC object from 10x directory
#'
#' @param dir 10x directory. Should have the following files: 'matrix.mtx', 'barcodes.tsv' and 'features.tsv'.
#' See 'Cell Ranger Arc' for more details: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/overview/welcome
#' @param matrix_fn if \code{dir} is missing, the filename of the matrix to import ("matrix.mtx" or "matrix.mtx.gz")
#' @param cells_fn if \code{dir} is missing, the filename of the cells to import ("barcodes.tsv" or "barcodes.tsv.gz")
#' @param features_fn if \code{dir} is missing, the filename of the features to import ("features.tsv" or "features.tsv.gz")
#' @param genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @param id an identifier for the object, e.g. "pbmc". If NULL - a random id would be assigned.
#' @param description (Optional) description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)".
#' @param metadata (Optional) data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#'
#'
#' @return an ScATAC object
#'
#' @examples
#' \dontrun{
#' atac <- import_from_10x("./pbmc_data", genome = "hg38")
#' atac
#' }
#'
#' @export
import_from_10x <- function(dir, genome, id = NULL, description = NULL, metadata = NULL, matrix_fn = file.path(dir, "matrix.mtx"), cells_fn = file.path(dir, "barcodes.tsv"), features_fn = file.path(dir, "features.tsv")) {
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
    non_zero_inds <- which(Matrix::rowSums(atac_mat) > 0)
    atac_peaks <- atac_peaks[non_zero_inds, ]
    n_all_zero <- nrow(atac_mat) - length(non_zero_inds)
    atac_mat <- atac_mat[non_zero_inds, ]

    if (n_all_zero > 0) {
        cli_alert_info("Removed {.val {nrow(atac_mat)}} ATAC peaks which were all zero")
    }

    cli_alert_info("{.val {nrow(atac_mat)}} ATAC peaks")
    res <- new("ScATAC", atac_mat, atac_peaks, genome = genome, id = id, description = description, metadata = metadata, path = matrix_fn)
    cli_alert_success("successfully imported to an ScATAC object with {.val {ncol(res@mat)}} cells and {.val {nrow(res@mat)}} ATAC peaks")

    return(res)
}
