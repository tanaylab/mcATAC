#' Given metacells (usually from RNA data), project ATAC counts to get a McATAC object
#'
#' @description Given cell to metacell association, summarise atac counts to generate a McATAC object.
#' \code{project_atac_on_mc_from_metacell1} projects given a metacell1 'mc' object, while \code{project_atac_on_mc_from_h5ad}
#' uses the output of the 'metacells' python package (metacell2).
#'
#'
#' @param atac an ScATAC object
#' @param cell_to_metacell a data frame with a column named "cell_id" with cell id and another column named "metacell" with the metacell the cell is
#' part of.
#' @param metadata per-metacell metadata. A data frame with a column called 'metacell' and additional metacell annotations.
#' @param min_int_frac (optional) minimal expected fraction of intersection of barcodes (cell names) in ScATAC
#' @param mc_size_eps_q (optional) quantile of MC size (in UMIs) to add as epsilon for regularizing accessibility across MCs
#'
#' @return an McATAC object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' atac_mc <- project_atac_on_mc_from_metacell1(atac_sc, "pbmc_data/scdb", "rna")
#' }
#'
#' @export
project_atac_on_mc <- function(atac, cell_to_metacell = NULL, metadata = NULL, min_int_frac = 0.5, mc_size_eps_q = 0.1) {
    cell_to_metacell <- deframe(cell_to_metacell)
    assert_that(all(names(cell_to_metacell) %in% colnames(atac@mat)))
    sc_mat <- atac@mat[, colnames(atac@mat) %in% names(cell_to_metacell), drop = FALSE]
    n_removed_cells <- ncol(atac@mat) - ncol(sc_mat)
    if (n_removed_cells > 0) {
        cli_alert_info("{.val {n_removed_cells}} cells (out of {.val {ncol(atac@mat)}}) do not have a metacell and have been removed.")
        if ((ncol(atac@mat) - n_removed_cells) <= round(min_int_frac * ncol(atac@mat))) {
            cli_abort("Intersect of ATAC mat colnames and mc names is less than {.field {scales::percent(min_int_frac)}}. Make sure you are projecting the right objects. To override - set {.code min_int_frac=0}")
        }
    }
    non_zero_peaks <- which(Matrix::rowSums(sc_mat) > 0)
    sc_sizes <- Matrix::colSums(sc_mat)
    mc_sizes <- tapply(sc_sizes, cell_to_metacell, sum)
    eps <- quantile(mc_sizes, mc_size_eps_q)
    mc_mat <- t(tgs_matrix_tapply(sc_mat[non_zero_peaks,], cell_to_metacell, mean))
    # mc_mat <- t(apply(mc_mat, 1, function(x) log2((x + eps)/median(x + eps))))
    assert_that(are_equal(atac@peaks$peak_name[non_zero_peaks], rownames(mc_mat)))
    assert_that(all(colnames(mc_mat) %in% cell_to_metacell))

    # TODO: deal with cell metadata
    # Naive solution - tabulate metadata per metacell and concatenate...

    mc_atac <- new("McATAC", mc_mat, atac@peaks[non_zero_peaks,], atac@genome, metadata)
    cli_alert_success("Created a new McATAC object with {.val {ncol(mc_atac@mat)}} metacells and {.val {nrow(mc_atac@mat)}} ATAC peaks.")

    return(mc_atac)
}

#' @param atac an ScATAC object
#' @param scdb a metacell1 \code{scdb} path
#' @param mc_id id of the metacell object within \code{scdb}
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_metacell1 <- function(atac, scdb, mc_id, metadata = NULL) {
    metacell::scdb_init(scdb, force_reinit = TRUE)
    rna_mc <- metacell::scdb_mc(mc_id)
    cell_to_metacell <- rna_mc@mc %>%
        enframe("cell_id", "metacell") %>%
        as_tibble()
    if (!is.null(rna_mc@colors)) {
        md <- enframe(rna_mc@colors, name = 'metacell', value = 'color')
    }
    if (!is.null(rna_mc@color_key)) {
        md <- enframe(rna_mc@colors, name = 'metacell', value = 'color')
        md$cell_type <- rna_mc@color_key$cell_type[match(md$color, rna_mc@color_key$color)]
    }
    return(project_atac_on_mc(atac, cell_to_metacell, metadata = md))
}

#' @param atac an ScATAC object
#' @param h5ad_file name of an h5ad file which is the output of 'metacells' python package.
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_h5ad <- function(atac, h5ad_file, metadata = NULL) {
    cli_alert_info("Reading {.file {h5ad_file}}")
    adata <- anndata::read_h5ad(h5ad_file)
    if (!("metacell" %in% colnames(adata$obs))) {
        cli_abort("h5ad object doesn't have a {.field metacell} field.")
    }

    cell_to_metacell <- adata$obs %>%
        rownames_to_column("cell_id") %>%
        select(cell_id, metacell) %>%
        as_tibble()

    # TODO: deal with cell metadata
    # TODO: test

    return(project_atac_on_mc(atac, cell_to_metacell))
}
