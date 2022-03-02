#' Given metacells (usually from RNA data), project ATAC counts to get a McATAC object
#'
#' @description Given cell to metacell association, summarise atac counts to generate a McATAC object.
#' \code{project_atac_on_mc_from_metacell1} projects given a metacell1 'mc' object, while \code{project_atac_on_mc_from_h5ad}
#' uses the output of the 'metacells' python package (metacell2).
#'
#'
#' @param atac an ScATAC object
#' @param cell_to_metacell a data frame with a column named "cell_id" with cell id and
#' another column named "metacell" with the metacell the cell is part of.
#'
#' @return an McATAC object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("raw")
#' atac_mc <- project_atac_on_mc_from_metacell1(atac_sc, "raw/scdb", "rna")
#' }
#'
#' @export
project_atac_on_mc <- function(atac, cell_to_metacell) {
    cell_to_metacell <- deframe(cell_to_metacell)
    sc_mat <- atac$mat[, colnames(atac$mat) %in% cell_to_metacell, drop = FALSE]
    mc_mat <- tgs_matrix_tapply(sc_mat, cell_to_metacell, sum)
    # TODO: deal with metadata (?)
    return(McATAC(mc_mat, atac$peaks))
}

#'
#' @param scdb a metacell1 \code{scdb} path
#' @param mc_id id of the metacell object within \code{scdb}
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_metacell1 <- function(atac, scdb, mc_id) {
    metacell::scdb_init(scdb, force_reinit = TRUE)
    rna_mc <- metacell::scdb_mc(mc_id)
    cell_to_metacell <- rna_mc@mc %>%
        enframe("cell_id", "metacell") %>%
        as_tibble()

    return(project_atac_on_mc(atac, cell_to_metacell))
}

#'
#' @param h5ad_file name of an h5ad file which is the output of 'metacells' python package.
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_h5ad <- function(atac, h5ad_file) {

}
