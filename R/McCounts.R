#' McCounts object
#'
#' @description a McCounts object is a collection of sparse matrices where rows are genomic coordinates and columns are metacells. It is derived from ScCounts using \code{project_counts_on_mc}. The metacell names are store in the \code{cell_names} slot, and the only additional slot is \code{cell_to_metacell} which contains the mapping from single cell names to metacell names.
#'
#' @slot cell_to_metacell
#'
#' @include ScCounts.R
#' @rdname ScCounts
#' @exportClass McCounts
McCounts <- setClass(
    "McCounts",
    slots = c(
        cell_to_metacell = "data.frame_or_null"
    ),
    contains = "ScCounts"
)

setMethod(
    "initialize",
    signature = "McCounts",
    definition = function(.Object, cell_to_metacell, ...) {
        .Object@cell_to_metacell <- cell_to_metacell
        callNextMethod(.Object, ...)
        return(.Object)
    }
)
#

#' @export
#' @noRd
setMethod(
    "show",
    signature = "McCounts",
    definition = function(object) {
        print_counts_object(object, "McCounts", "metacell")
    }
)

#' Given metacells (usually from RNA data), project ATAC counts to get a McCounts object
#' 
#' @description Given cell to metacell association, summarise atac read counts to generate a McCounts object. This can 
#' take a while - around 5 minutes using 24 cores on the PBMC dataset.
#'
#' @param sc_counts A ScCounts object
#' @param cell_to_metacell a data frame with a column named "cell_id" with cell id and another column named "metacell" with the metacell the cell is part of.
#' @param num_cores The number of cores to use.
#'
#' @return A McCounts object
#'
#' @examples
#' \dontrun{
#' data(cell_to_metacell_pbmc_example)
#' project_counts_on_mc(sc_counts, cell_to_metacell_pbmc_example)
#' }
#'
#' @export
project_counts_on_mc <- function(sc_counts, cell_to_metacell, num_cores = parallel::detectCores()) {
    cell_to_metacell <- deframe(cell_to_metacell)
    cell_to_metacell <- cell_to_metacell[cell_to_metacell %!in% c(-2, -1, 0)] # make sure you don't have any outlier metacells
    assert_that(all(names(cell_to_metacell) %in% sc_counts@cell_names))


    removed_cells <- setdiff(sc_counts@cell_names, names(cell_to_metacell))
    if (length(removed_cells) > 0) {
        cli_alert_info("{.val {length(removed_cells)}} cells (out of {.val {length(sc_counts@cell_names)}}) do not have a metacell and have been removed.")
    }

    cell_to_metacell <- cell_to_metacell[intersect(names(cell_to_metacell), sc_counts@cell_names)]

    doMC::registerDoMC(num_cores)
    new_data <- plyr::llply(sc_counts@data, function(sc_mat) {
        cells <- intersect(colnames(sc_mat), names(cell_to_metacell))
        mc_mat <- t(tgs_matrix_tapply(sc_mat[, cells], cell_to_metacell[cells], sum))
        return(mc_mat)
    }, .parallel = TRUE)

    res <- new("McCounts", data = new_data, cell_names = sort(unique(cell_to_metacell)), genome = sc_counts@genome, chromosomes = sc_counts@chromosomes, id = sc_counts@id, description = sc_counts@description, path = sc_counts@path, cell_to_metacell = enframe(cell_to_metacell, "cell_id", "metacell"))

    return(res)
}
