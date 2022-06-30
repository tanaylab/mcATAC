#' Add per-metacell gene expression data to an McATAC object
#'
#' @param atac_mc an McATAC object
#' @param mc_rna a \code{metacell1} 'mc' object or a \code{metacells} metacell UMI matrix (a matrix where each row is a gene and each column is a metacell)
#'
#' @examples
#' \dontrun{
#' data(mcmd)
#' metacell::scdb_init("pbmc_data/scdb", force_reinit = TRUE)
#' mc_rna <- metacell::scdb_mc("rna")
#' atac_mc <- add_mc_rna(atac_mc, mc_rna)
#' }
#'
#' @export
add_mc_rna <- function(atac_mc, mc_rna) {
    assert_atac_object(atac_mc, class = "McATAC")

    if ("tgMCCov" %in% class(mc_rna)) {
        egc <- mc_rna@e_gc
    } else if (is.matrix(mc_rna) || is_sparse_matrix(mc_rna)) {
        egc <- t(t(mc_rna) / colSums(mc_rna))
    } else {
        cli_abort("mc_rna must be a tgMCCov object (from the metacell package) or a matrix.")
    }

    rna_mc_not_in_atac <- colnames(egc)[colnames(egc) %!in% colnames(atac_mc@mat)]
    if (length(rna_mc_not_in_atac) > 0) {
        cli_warn("{.field mc_rna} contains {.field {length(rna_mc_not_in_atac)}} metacells not present in the McATAC object: {.val {rna_mc_not_in_atac}}")
    }

    atac_mc_not_in_rna <- colnames(atac_mc@mat)[colnames(atac_mc@mat) %!in% colnames(egc)]
    if (length(atac_mc_not_in_rna) > 0) {
        cli_warn("McATAC object contains {.field {length(atac_mc_not_in_rna)}} metacells not present in {.field mc_rna}: {.val {atac_mc_not_in_rna}}")
    }

    both_mcs <- intersect(colnames(atac_mc@mat), colnames(egc))
    if (length(both_mcs) == 0) {
        cli_abort("No metacells in common between the McATAC object and {.field mc_rna}.")
    }

    atac_mc@rna_egc <- egc[, both_mcs]
    return(atac_mc)
}

#' Does the McATAC object contain per-metacell gene expression data?
#'
#' @param atac_mc an McATAC object
#'
#' @return TRUE if the McATAC object contains per-metacell gene expression data, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_rna(atac_mc)
#' }
#'
#' @export
has_rna <- function(atac_mc) {
    return(!is.null(nrow(atac_mc@rna_egc)) && nrow(atac_mc@rna_egc) > 0)
}


#' Get the normalized RNA expression matrix (egc) from a McATAC object
#'
#' @description Get RNA expression data, normalized by the total RNA expression in each metacell.
#'
#' @param atac_mc a McATAC object with RNA expression (using \code{add_mc_rna})
#' @param genes list of genes to match. Default (NULL): all genes
#' @param rm_zeros remove genes with no RNA expression in any metacell. Default: TRUE
#' @param epsilon regularization factor added to the log normalized expression
#'
#' @return a matrix with normalized RNA expression values for each gene (rows) and metacell (columns)
#'
#' @examples
#' \dontrun{
#' rna_mat <- get_rna_egc(atac_mc)
#' rna_mat <- get_rna_egc(atac_mc, genes = c("MESP1", "MESP2", "PF4"))
#' rna_mat <- get_rna_egc(atac_mc, rm_zeros = FALSE)
#' }
#'
#' @export
get_rna_egc <- function(atac_mc, genes = NULL, rm_zeros = TRUE, epsilon = 1e-5) {
    validate_atac_object(atac_mc)
    if (!has_rna(atac_mc)) {
        cli_abort("{.val {atac_mc}} does not contain RNA.")
    }
    rna_mat <- atac_mc@rna_egc
    if (!is.null(genes)) {
        if (any(genes %!in% rownames(rna_mat))) {
            missing_genes <- genes[genes %!in% rownames(rna_mat)]
            cli_abort("Genes {.val {missing_genes}} are not a subset of {.field atac_mc@rna_egc}.")
        }
        rna_mat <- rna_mat[genes, , drop = FALSE]
    }

    if (rm_zeros) {
        rna_mat <- rm_zero_expr_genes(rna_mat)
    }

    if (!is.null(epsilon)) {
        rna_mat <- rna_mat + epsilon
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


#' Get enrichment of normalized counts of gene expression over the median
#'
#' @description
#' The function first transforms the UMI matrix to fractions per metacell, and then calculates the enrichment of each gene
#' over the median (dividing the gene expression by the median).
#'
#' @return a matrix with normalized counts of gene expression for each gene (rows) and metacell (columns)
#'
#' @inheritParams get_rna_egc
#'
#' @examples
#' \dontrun{
#' get_rna_fp(atac_mc)
#' get_rna_fp(atac_mc, "GNLY")
#' get_rna_fp(atac_mc, "GNLY", epsilon = NULL)
#' }
#'
#' @export
get_rna_fp <- function(atac_mc, genes = NULL, rm_zeros = TRUE, epsilon = 1e-5) {
    rna_egc <- get_rna_egc(atac_mc, genes = genes, rm_zeros = rm_zeros, epsilon = epsilon)
    rna_fp <- rna_egc / apply(rna_egc, 1, median, na.rm = TRUE)
    return(rna_fp)
}

#' Calculate gene markers of metacells gene expression
#'
#' @description
#' The function first removes genes without sufficient expression in any metacell \code{minimal_max_log_fraction}, or without sufficient fold
#' change over the median (\code{minimal_relative_expression}), and then ranks the fold change of genes within each metacell. The markers
#' are then the genes with the highest rank, up to \code{n_genes} genes.
#'
#' @param n_genes maximal number of genes to return. Default: 100
#' @param minimal_max_log_fraction take only genes with at least one value
#' (in log fraction units - normalized egc) above this threshold
#' @param minimal_relative_log_fraction take only genes with at least one value with relative
#' log fraction (mc_fp) above this this value
#' @param fold_change_reg regularization factor for the fold change calculation (fold_change would be changed to
#' \code{fold_change = fold_change + fold_change_reg})
#'
#' @inheritParams get_rna_egc
#'
#' @examples
#' \dontrun{
#' get_rna_markers(atac_mc)
#' }
#'
#' @export
get_rna_markers <- function(atac_mc, n_genes = 100, minimal_max_log_fraction = -13, minimal_relative_log_fraction = 2, fold_change_reg = 0.1, genes = NULL, rm_zeros = TRUE, epsilon = 1e-5) {
    mc_egc <- log2(get_rna_egc(atac_mc, genes = genes, rm_zeros = rm_zeros, epsilon = epsilon))

    max_log_fractions_of_genes <- sparseMatrixStats::rowMaxs(mc_egc, na.rm = TRUE)
    interesting_genes_mask <- max_log_fractions_of_genes >= minimal_max_log_fraction
    cli_alert("removing {.val {sum(!interesting_genes_mask)}} genes with no RNA expression (log2) of above {.field {minimal_max_log_fraction}} in any metacell.")

    mc_egc <- mc_egc[interesting_genes_mask, , drop = FALSE]
    fold_matrix <- sweep(mc_egc, 1, sparseMatrixStats::rowMedians(mc_egc, na.rm = TRUE))
    fold_matrix <- fold_matrix + fold_change_reg
    max_relative_log_fractions_of_genes <- sparseMatrixStats::rowMaxs(fold_matrix, na.rm = TRUE)
    interesting_genes_mask <- max_relative_log_fractions_of_genes >= minimal_relative_log_fraction
    cli_alert("removing {.val {sum(!interesting_genes_mask)}} genes with no fold change (log2) of above {.field {minimal_relative_log_fraction}} in any metacell.")
    fold_matrix[interesting_genes_mask, , drop = FALSE]

    cli_alert_info("{.val {nrow(fold_matrix)}} genes left for consideration.")

    if (nrow(fold_matrix) < n_genes) {
        cli_abort("not enough genes left to return {.val {n_genes}} markers. Please try to set a lower value for {.field minimal_max_log_fraction} or {.field minimal_relative_log_fraction}.")
    }

    if (nrow(fold_matrix) > n_genes) {
        gene_ranks <- sparseMatrixStats::colRanks(fold_matrix, useNames = TRUE, preserveShape = TRUE, ties.method = "first")
        max_gene_ranks <- sparseMatrixStats::rowMaxs(gene_ranks, na.rm = TRUE, useNames = TRUE)
        markers <- sort(max_gene_ranks, decreasing = TRUE) %>%
            head(n = n_genes) %>%
            names()
    } else {
        markers <- rownames(fold_matrix)
    }

    cli_alert_success("{.val {length(markers)}} marker genes selected.")

    return(markers)
}

#' Get enrichment matrix for marker genes
#'
#' @param atac_mc a McATAC object with RNA expression (using \code{add_mc_rna})
#' @param markers a list of marker genes. If NULL - the function uses \code{get_rna_markers} with default parameters which can be overridden
#' using the ellipsis \code{...}.
#' @param force_cell_type do not split cell types when ordering the metacells. Default: TRUE
#' @inheritParams get_rna_fp
#' @inheritDotParams get_rna_markers
#'
#' @return A matrix with log2 normalized counts of gene expression for each marker gene (rows) and metacell (columns). \cr
#' The columns are columns are clustered using \code{hclust} with "ward" linkage, on the euclidean distance
#' of the pearson correlation matrix between the columns. Then, the gene with the highest number of metacells with fold change (abs) above 1 (gene1) is selected together with the gene most anti-correlated to it (gene2), and the \code{hclust} is ordered using the difference (gene1 - gene2). \cr
#' If \code{force_cell_type} is TRUE and \code{atac_mc@metadata} has a field named "cell_type", the columns are only ordered only within each
#' cell type.
#'
#'
#' @examples
#' \dontrun{
#' marker_mat <- get_rna_marker_matrix(atac_mc)
#' marker_mat <- get_rna_marker_matrix(atac_mc, n_genes = 100)
#' }
#'
#' @export
get_rna_marker_matrix <- function(atac_mc, markers = NULL, force_cell_type = TRUE, rm_zeros = TRUE, epsilon = 1e-5, ...) {
    if (is.null(markers)) {
        markers <- get_rna_markers(atac_mc, ...)
    }
    rna_fp <- get_rna_fp(atac_mc, genes = markers, rm_zeros = rm_zeros, epsilon = epsilon)

    if (force_cell_type && has_cell_type(atac_mc)) {
        metacell_types <- atac_mc@metadata %>% select(metacell, cell_type)
    } else {
        metacell_types <- NULL
    }

    rna_fp <- order_marker_matrix(rna_fp, metacell_types = metacell_types)

    cli_alert_success("marker matrix of {.val {nrow(rna_fp)}} genes x {.val {ncol(rna_fp)}} metacells created.")
    return(log2(rna_fp))
}

order_marker_matrix <- function(mat, metacell_types = NULL) {
    g_ncover <- rowSums(abs(mat) > 1, na.rm = TRUE)
    main_mark <- names(g_ncover)[which.max(g_ncover)]
    f <- mat[main_mark, , drop = FALSE] < 0.25
    if (sum(f) > 0.2 * ncol(mat)) {
        g_score <- apply(mat[, f, drop = FALSE] > 1, 1, sum)
    } else {
        g_score <- -tgs_cor(t(mat), t(mat[main_mark, , drop = FALSE]))[, 1]
    }
    second_mark <- names(g_score)[which.max(g_score)]
    cli_alert_info("Ordering metacells based on {.file {main_mark}} vs {.file {second_mark}}")

    zero_mcs <- colSums(abs(mat) > 0) < 2
    if (any(zero_mcs)) {
        mat_all <- mat
        mat <- mat_all[, !zero_mcs, drop = FALSE]
        mat_zero <- mat_all[, zero_mcs, drop = FALSE]
    }

    if (ncol(mat) == 0) { # all the metacells do not have enough non-zero values
        return(1:ncol(mat_all))
    }

    hc <- stats::hclust(tgs_dist(tgs_cor(mat, pairwise.complete.obs = TRUE)), method = "ward.D2")

    d <- stats::reorder(
        stats::as.dendrogram(hc),
        mat[main_mark, ] - mat[second_mark, ],
        agglo.FUN = mean
    )
    ord <- as.hclust(d)$order

    if (any(zero_mcs)) {
        mc_order <- c(colnames(mat)[ord], colnames(mat_zero))
        mat <- mat_all
    } else {
        mc_order <- colnames(mat)[ord]
    }

    if (!is.null(metacell_types)) {
        cli_alert_info("Maintaining metacell order within cell types")
        metacell_types <- metacell_types %>% mutate(metacell = as.character(metacell))
        ord_df <- tibble(metacell = colnames(mat)) %>%
            mutate(orig_ord = 1:n()) %>%
            left_join(
                tibble(metacell = as.character(mc_order)) %>% mutate(glob_ord = 1:n()),
                by = "metacell"
            ) %>%
            left_join(
                metacell_types %>% select(metacell, cell_type),
                by = "metacell"
            )

        ord <- ord_df %>%
            mutate(orig_ord = 1:n()) %>%
            group_by(cell_type) %>%
            mutate(ct_ord = mean(glob_ord)) %>%
            ungroup() %>%
            arrange(ct_ord, glob_ord) %>%
            pull(orig_ord)
    }

    mat <- mat[, ord, drop = FALSE]

    gene_ord <- order(apply(mat, 1, which.max))
    mat <- mat[gene_ord, , drop = FALSE]
    return(mat)
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
#' @inheritParams get_rna_egc
#' @inheritParams tgstat::tgs_cor_knn
#' @export
rna_atac_cor_knn <- function(atac_mc, k = 1, genes = NULL, rm_zeros = TRUE, spearman = TRUE, pairwise.complete.obs = TRUE) {
    assert_atac_object(atac_mc, "McATAC")

    rna_mat <- get_rna_egc(atac_mc, genes = genes, rm_zeros = rm_zeros)

    if ("matrix" %!in% class(atac_mc@mat)) {
        atac_mc@mat <- as.matrix(atac_mc@mat)
    }
    knn_df <- tgs_cor_knn(t(rna_mat), t(atac_mc@mat), knn = k, spearman = spearman, pairwise.complete.obs = pairwise.complete.obs)

    knn_df <- knn_df %>%
        rename(gene = col1, peak = col2) %>%
        as_tibble()

    return(knn_df)
}

#' Return a relative fold change matrix of ATAC peaks for a list of genes
#'
#' @description
#' This function returns a relative fold change matrix of ATAC peaks for a list of genes matched using \code{rna_atac_cor_knn}.
#'
#' @param atac_mc a McATAC object
#' @param genes a list of genes.
#' @param metacell select only a subset of the metacells.
#'
#' @inheritParams rna_atac_cor_knn
#'
#' @examples
#' \dontrun{
#' marker_genes <- get_rna_markers(atac_mc)
#' atac_fp <- get_genes_atac_fp(atac_mc, genes = marker_genes)
#' rna_fp <- get_rna_marker_matrix(atac_mc, genes = marker_genes)
#' }
#'
#' @export
get_genes_atac_fp <- function(atac_mc, genes = NULL, metacells = NULL, rm_zeros = TRUE, spearman = TRUE, pairwise.complete.obs = TRUE) {
    assert_atac_object(atac_mc, "McATAC")
    knn_df <- rna_atac_cor_knn(atac_mc, genes = genes, rm_zeros = rm_zeros, spearman = spearman, pairwise.complete.obs = pairwise.complete.obs)

    atac_fp <- atac_mc@fp[knn_df$peak, ]

    if (!is.null(metacells)) {
        atac_fp <- atac_fp[, metacells, drop = FALSE]
    }

    return(atac_fp)
}
