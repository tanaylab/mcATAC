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

    if ("tgMCCov" %in% class(mc_rna)) {
        egc <- mc_rna@e_gc
    } else if (is.matrix(mc_rna) || is_sparse_matrix(mc_rna)) {
        egc <- t(t(mc_rna) / colSums(mc_rna))
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

#' Does the McATAC object contain per-metacell gene expression data?
#'
#' @param mc_atac an McATAC object
#'
#' @return TRUE if the McATAC object contains per-metacell gene expression data, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_rna(mc_atac)
#' }
#'
#' @export
has_rna <- function(mc_atac) {
    return(!is.null(nrow(mc_atac@rna_egc)) && nrow(mc_atac@rna_egc) > 0)
}

#' Get the RNA expression matrix from a McATAC object
#'
#' @param atac_mc a McATAC object with RNA expression (using \code{add_mc_rna})
#' @param genes list of genes to match. Default (NULL): all genes
#' @param rm_zeros remove genes with no RNA expression in any metacell. Default: TRUE
#'
#' @return a matrix with RNA expression values for each gene (rows) and metacell (columns)
#'
#' @examples
#' \dontrun{
#' rna_mat <- get_rna_matrix(atac_mc)
#' rna_mat <- get_rna_matrix(atac_mc, genes = c("MESP1", "MESP2", "PF4"))
#' rna_mat <- get_rna_matrix(atac_mc, rm_zeros = FALSE)
#' }
#'
#' @export
get_rna_matrix <- function(atac_mc, genes = NULL, rm_zeros = TRUE) {
    if (!has_rna(atac_mc)) {
        cli_abort("{.val {atac_mc}} does not contain RNA.")
    }
    rna_mat <- atac_mc@rna_egc
    if (!is.null(genes)) {
        if (any(genes %!in% rownames(rna_mat))) {
            missing_genes <- genes[genes %!in% rownames(rna_mat)]
            cli_abort("Genes {.val {missing_genes}} are not a subset of {.field atac_mc@rna_egc}.")
        }
        rna_mat <- rna_mat[genes, ]
    }

    if (rm_zeros) {
        rna_mat <- rm_zero_expr_genes(rna_mat)
    }
    return(rna_mat)
}

rm_zero_expr_genes <- function(rna_mat) {
    f <- rowSums(rna_mat, na.rm = TRUE) == 0
    if (sum(f) > 0) {
        cli_alert("removing {.field {sum(f)}} genes with no RNA expression in any metacell.")
        rna_mat <- rna_mat[!f, ]
    }
    return(rna_mat)
}


#' Match every gene with the k ATAC peaks most correlated to it
#'
#' @param atac_mc a McATAC object with RNA expression (using \code{add_mc_rna})
#' @param k number of peaks to match for each gene. Default: 1
#'
#' @return a tibble with the following columns:
#' \itemize{
#'  \item{gene: }{Gene name}
#'  \item{peak: }{ATAC peak name}
#'  \item{cor: }{Correlation between ATAC and RNA}
#'  \item{rank: }{Rank of the correlation (for the gene)}
#' }
#'
#' @examples
#' \dontrun{
#' rna_atac_cor_knn(atac_mc, k = 1)
#' rna_atac_cor_knn(atac_mc, k = 1, genes = c("MESP1", "MESP2", "PF4"))
#' }
#'
#' @inheritParams get_rna_matrix
#' @inheritParams tgstat::tgs_cor_knn
#' @export
rna_atac_cor_knn <- function(atac_mc, k = 1, genes = NULL, rm_zeros = TRUE, spearman = TRUE, pairwise.complete.obs = TRUE) {
    assert_atac_object(atac_mc, "McATAC")

    rna_mat <- get_rna_matrix(atac_mc, genes = genes, rm_zeros = rm_zeros)

    knn_df <- tgs_cor_knn(t(rna_mat), t(atac_mc@mat), knn = k, spearman = spearman, pairwise.complete.obs = pairwise.complete.obs)

    knn_df <- knn_df %>%
        rename(gene = col1, peak = col2) %>%
        as_tibble()

    return(knn_df)
}
