#' Given metacells (usually from RNA data), project ATAC counts to get a McATAC object
#'
#' @description Given cell to metacell association, summarise atac counts to generate a McATAC object.
#' \code{project_atac_on_mc_from_metacell1} projects given a metacell1 'mc' object, while \code{project_atac_on_mc_from_h5ad}
#' uses the output of the 'metacells' python package (metacell2).
#'
#'
#' @param atac an ScATAC object
#' @param mc a metacell1 \code{mc} object
#' @param cell_to_metacell a data frame with a column named "cell_id" with cell id and
#' @param MIN_INT_FRAC (optional) minimal expected fraction of intersection of barcodes (cell names) in ScATAC and mc/cell_to_metacell
#' another column named "metacell" with the metacell the cell is part of.
#'
#' @export
project_atac_on_mc <- function(atac, mc = NULL, cell_to_metacell = NULL, MIN_INT_FRAC = 0.1) {
    if (!is.null(mc) && is.null(cell_to_metacell)) {
        cells_both = intersect(colnames(atac$mat), names(mc@mc))
        mc_filt = mc@mc[names(mc@mc) %in% cells_both]
    }
    if (is.null(mc) && !is.null(cell_to_metacell)) {
        cells_both = intersect(colnames(atac$mat), cell_to_metacell$cell_id)
        c2m_ord = cell_to_metacell[match(cells_both, cell_to_metacell$cell_id),]
        mc_filt = setNames(c2m_ord$metacell, c2m_ord$cell_id)
    }
    if (length(cells_both) <= round(0.1*length(colnames(atac$mat)))) {
                warning(glue::glue('Intersect of ATAC mat colnames and mc names is less than {round(100*MIN_INT_FRAC, 2)}%. Make sure you are projecting the right objects.'))
            }
    atac_mat_filt = as.matrix(atac$mat[,cells_both])
    mc_filt = mc_filt[cells_both]
    mc_atac = tgs_matrix_tapply(atac_mat_filt, mc_filt, mean)
    return(t(mc_atac))
}

#' @param atac an ScATAC object
#' @param mc a metacell1 \code{mc} object
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_metacell1 <- function(atac, mc = NULL, cell_to_metacell = NULL) {
    
    
}

#' @param atac an ScATAC object
#' @param h5ad_file name of an h5ad file which is the output of 'metacells' python package.
#'
#' @export
#' @rdname project_atac_on_mc
project_atac_on_mc_from_h5ad <- function(atac, h5ad_file) {

}
