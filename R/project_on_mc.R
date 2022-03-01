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
#' @export
project_atac_on_mc <- function(atac, cell_to_metacell) {

}

#'
#' @param mc a metacell1 \code{mc} object
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_metacell1 <- function(atac, mc) {

}

#'
#' @param h5ad_file name of an h5ad file which is the output of 'metacells' python package.
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_h5ad <- function(atac, h5ad_file) {

}
