#' Plot a scatter of gene expression vs an atac profile
#'
#' @param mc_atac a McATAC object
#' @param gene name of the gene to plot
#' @param mc_rna a \code{metacell1} 'mc' object or a \code{metacells} metacell UMI matrix (a matrix where each row is a gene and each column is a metacell). Can be NULL if \code{mc_atac} already contains the gene expression data (added by \code{add_mc_rna}).
#' @param peak name of the peak to plot. If NULL - the promoter of \code{gene} would be shown. You can get the peak names from \code{atac@peaks$peak_name}
#' @param max_dist_to_promoter_peak how far from \code{gene}'s TSS to search for a promoter-proximal peak
#' @param normalize_atac whether to use normalized atac profiles (default: TRUE)
#' @param eps_rna added regularization when calculating log expression (Default: 1e-5). Promoter ATAC signal has been observed empirically to often be linearly correlated with log of gene expression.
#'
#' @return a ggplot object with the scatter plot
#'
#' @examples
#' \dontrun{
#' p1 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4")
#'
#' ## Plot gene vs some peak of interest; here - expression of CD4 vs. promoter ATAC signal of CD8
#' tss <- gintervals.load("intervs.global.tss")
#' cd8_tss <- tss[which(tss$geneSymbol == "CD8")[[1]], ]
#' peaks_near_cd8_tss <- misha::gintervals.neighbors(pbmc_atac_mc@peaks, cd8_tss, mindist = -2e+3, maxdist = 2e+3)
#' nearest_peak <- peaks_near_cd8_tss[which.min(peaks_near_cd8_tss$dist), ]
#' p2 <- plot_atac_rna(mc_atac = pbmc_atac_mc, mc_rna = pbmc_rna_mc, gene = "CD4", peak = nearest_peak[, "peak_name"])
#' }
#' @export
plot_atac_rna <- function(mc_atac, gene, mc_rna = NULL, peak = NULL, max_dist_to_promoter_peak = 5e+2, normalize_atac = TRUE, eps_rna = 1e-5, tss_intervals = "intervs.global.tss") {
    if (has_rna(mc_atac)) {
        if (is.null(mc_rna)) {
            cli_abort("No gene expression data in the McATAC object. Use {.code add_mc_rna()} to add it or provide the {.field mc_rna} parameter.")
        } else {
            mc_atac <- add_mc_rna(mc_atac, mc_rna)
        }
    }

    rv <- log2(mc_atac@rna_egc[gene, ] + eps_rna)

    if (is.null(peak)) {
        if (!gintervals.exists(tss_intervals)) {
            cli_abort("{.val {tss_intervals}} intervals do not exist. You can either define the peak explicitly or import them from ucsc")
        }
        tss <- gintervals.load(tss_intervals)
        tss_gene <- tss[tss$geneSymbol == gene, c("chrom", "start", "end", "geneSymbol")]
        if (nrow(tss_gene) == 0) {
            cli_abort("Peak not supplied and gene TSS location not found")
        } else if (nrow(tss_gene) > 1) {
            cli_alert("The gene {.val {gene}} has {.val {nrow(tss_gene)}} alternative promoters. Summing the ATAC signal from all of them")
        }

        nei_peaks_tss <- misha.ext::gintervals.neighbors1(
            tss_gene[, 1:3],
            as.data.frame(mc_atac@peaks),
            mindist = -max_dist_to_promoter_peak,
            maxdist = max_dist_to_promoter_peak,
            maxneighbors = 500
        ) %>%
            filter(!is.na(peak_name))

        peak <- unique(nei_peaks_tss$peak_name)
        peak_range <- misha.ext::convert_10x_peak_names_to_misha_intervals(peak) %>%
            summarise(chrom = chrom[1], start = min(start), end = max(end)) %>%
            peak_names()
        peak_str <- glue("{peak_range} ({length(peak)} peaks)")

        if (length(peak) == 0) {
            cli_abort("{.field peak} was not supplied and no relevant promoter peak was found for gene {.val {gene}}. Check if your ATAC matrix should have a peak for this gene's promoter.")
        } else if (length(peak) > 1) {
            cli_alert("The gene {.val {gene}} has multiple ({.val {length(peak)}}) peaks within {.val {max_dist_to_promoter_peak}} bp of its TSS. Summing the ATAC signal from all of them.")
        }
    } else {
        peak_str <- peak
    }

    if (normalize_atac) {
        av <- colSums(mc_atac@egc[peak, ])
        ylab <- "ATAC (normalized)"
    } else {
        av <- colSums(mc_atac@mat[peak, ])
        ylab <- "ATAC (not normalized)"
    }

    # take only metacells which exist in both atac and rna
    both_mcs <- intersect(names(av), names(rv))
    av <- av[both_mcs]
    rv <- rv[both_mcs]

    cor_r2 <- cor.test(rv, av)$estimate

    df <- tibble(
        metacell = names(av),
        atac = av,
        rna = rv
    )

    if (all(has_name(mc_atac@metadata, c("cell_type", "color")))) {
        df <- df %>%
            left_join(mc_atac@metadata %>%
                mutate(metacell = as.character(metacell)) %>%
                select(metacell, cell_type, color), by = "metacell")
        gg <- ggplot(df, aes(x = rna, y = atac, fill = cell_type)) +
            geom_point(shape = 21) +
            scale_fill_manual(name = "Cell type", values = get_cell_type_colors(mc_atac@metadata))
    } else {
        gg <- ggplot(df, aes(x = rna, y = atac)) +
            geom_point()
    }

    gg <- gg +
        labs(
            title = gene,
            subtitle = glue("ATAC of {peak_str} vs. RNA, R^2 = {round(cor_r2, digits=2)}"),
            x = "log2(gene expression)",
            y = ylab,
            caption = glue("object id: {mc_atac@id}")
        ) +
        theme(plot.subtitle = ggtext::element_markdown())

    return(gg)
}

#' Plot a correlation matrix of ATAC metacells
#'
#' @param mc_atac McATAC object
#' @param sp_f whether to use Spearman correlation (default) or Pearson
#'
#' @return p a pheatmap of ATAC metacell correlations
#'
#' @examples
#' \dontrun{
#' p1 <- plot_atac_atac_cor(my_atac_mc)
#' p2 <- plot_atac_atac_cor(my_atac_mc, sp_f = F)
#' }
#' @export
plot_atac_atac_cor <- function(mc_atac, sp_f = TRUE) {
    if (all(!grepl("cell_type", colnames(mc_atac@metadata)))) {
        if (all(!grepl("cluster_k_", colnames(mc_atac@peaks)))) {
            k <- round(ncol(mc_atac@mat) / 10)
            message(glue::glue('There is no "cell_type" or "cluster" field in metadata, clustering with k == {k}'))
            mc_atac <- gen_atac_mc_clust(mc_atac, k = k, use_prior_annot = F)
            clust_vec <- unlist(mc_atac@metadata[, paste0("cluster_k_", k)])
        } else {
            clust_vec <- unlist(mc_atac@metadata[, grep("cluster_k_", colnames(mc_atac@metadata))[[1]]])
        }
    } else {
        clust_vec <- unlist(mc_atac@metadata[, "cell_type"])
    }
    col_annot <- as.data.frame(list("cluster" = clust_vec))
    rownames(col_annot) <- 1:nrow(mc_atac@metadata)
    if (any(grepl("color", colnames(mc_atac@metadata)))) {
        col_key <- unique(mc_atac@metadata[, c("cell_type", "color")])
        ann_colors <- list("cluster" = setNames(col_key$color, col_key$cell_type))
    } else {
        ann_colors <- list("cluster" = setNames(sample(
            grep("white|gray|grey", colors(), v = T, inv = T),
            length(unique(clust_vec))
        ), unique(clust_vec)))
    }
    cor_mat <- tgs_cor(mc_atac@mat, spearman = sp_f, pairwise.complete.obs = TRUE)
    p <- pheatmap::pheatmap(cor_mat[order(clust_vec), order(clust_vec)],
        cluster_rows = FALSE, cluster_cols = FALSE,
        show_rownames = FALSE, show_colnames = FALSE,
        annotation_col = col_annot, annotation_colors = ann_colors
    )
    return(p)
}

#' Plot a cross-correlation matrix between RNA metacells and ATAC metacell scores
#'
#' @description use peak gene annotations to match between RNA metacells and ATAC metacells and plot
#' a cross-correlation matrix.
#'
#' @param mc_atac McATAC object
#' @param rna_mat an RNA metacell count matrix, where metacells are in columns and genes are in rows
#' @param gene_field name of a field in \code{mc_atac@peaks} which contains the gene names. If NULL - the peaks would be
#' transformed to promoter peaks and the gene names would be taken from the promoter gene names.
#'
#' @export
plot_atac_rna_cor <- function(mc_atac, rna_mat) {

}

#' Plot normalized accessibility of peaks over metacells, ordered by clustering
#'
#' @param mc_atac McATAC object
#' @param mc_atac_clust output of \code{gen_atac_mc_clust}
#' @param peak_clust output of \code{gen_atac_peak_clust}
#' @examples
#' \dontrun{
#'
#' }
#' @export
plot_atac_peak_map <- function(mc_atac, mc_atac_clust, peak_clust) {
    # the central heat map showing normalized accessibility of peaks over metacells, ordered by clustering
}
