#' Add per-metacell gene expression data to an McATAC object
#'
#' @param mcatac an McATAC object
#' @param mc_rna a \code{metacell1} 'mc' object or a \code{metacells} metacell UMI matrix (a matrix where each row is a gene and each column is a metacell)
#'
#' @examples
#' \dontrun{
#' data(mcmd)
#' metacell::scdb_init("pbmc_data/scdb", force_reinit = TRUE)
#' mc_rna <- metacell::scdb_mc("rna")
#' mc_atac <- add_mc_rna(mc_atac, mc_rna)
#' }
#'
#' @export
add_mc_rna <- function(mcatac, mc_rna) {
    assert_atac_object(mcatac, class = "McATAC")

    if (class(mc) == "tgMCCov") {
        egc <- mc_rna@e_gc
    } else if (is.matrix(mc_rna) || is_sparse_matrix(mc_rna)) {
        egc <- t(t(mc_rna) / colSums(mc_sum))
    } else {
        cli_abort("mc_rna must be a tgMCCov object (from the metacell package) or a matrix.")
    }

    rna_mc_not_in_atac <- colnames(egc)[colnames(egc) %!in% colnames(mcatac@mat)]
    if (length(rna_mc_not_in_atac) > 0) {
        cli_warn("{.field mc_rna} contains {.field {length(rna_mc_not_in_atac)}} metacells not present in the McATAC object: {.val {rna_mc_not_in_atac}}")
    }

    atac_mc_not_in_rna <- colnames(mcatac@mat)[colnames(mcatac@mat) %!in% colnames(egc)]
    if (length(atac_mc_not_in_rna) > 0) {
        cli_warn("McATAC object contains {.field {length(atac_mc_not_in_rna)}} metacells not present in {.field mc_rna}: {.val {atac_mc_not_in_rna}}")
    }

    both_mcs <- intersect(colnames(mcatac@mat), colnames(egc))
    if (length(both_mcs) == 0) {
        cli_abort("No metacells in common between the McATAC object and {.field mc_rna}.")
    }

    mcatac@rna_egc <- egc[, both_mcs]
    return(mcatac)
}

has_rna <- function(mc_atac) {
    return(is.null(nrow(mc_atac@rna_egc)) || nrow(mc_atac@rna_egc) == 0)
}
