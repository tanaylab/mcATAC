
#' Select only peaks that are close to a gene promoter
#'
#' @description Select only peaks that are close to a gene promoter. In case of alternative promoters or multiple peaks close to the same promoter, the reads from all the peaks and promoters are summed.
#'
#' @param atac a McATACPeaks or ScATACPeaks object
#' @param id_field the field in \code{tss_intervals} containing the gene names. Default: "geneSymbol"
#'
#' @return the McATACPeaks or ScATACPeaks object with new peaks which are a sum of the peaks within the defined promoters. The rownames of the \code{atac@mat} and the field 'peak_name' in \code{atac@peaks} would contain the name of the promoter (from \code{id_field}).
#' Note that previous peak metadata and previous ignored peaks are dropped.
#'
#' @inheritParams misha.ext::get_promoters
#'
#' @examples
#' \dontrun{
#' atac_mc_promoters <- gen_promoter_features(atac_mc)
#' }
#'
#' @export
gen_promoter_features <- function(atac, upstream = 500, downstream = 50, tss_intervals = "intervs.global.tss", id_field = "geneSymbol") {
    assert_atac_object(atac)
    if (atac@promoters) {
        cli_abort("{.val {atac}} peaks are already promoters.")
    }

    promoters <- misha.ext::get_promoters(upstream = upstream, downstream = downstream, tss_intervals = tss_intervals)

    peak_map <- atac@peaks %>%
        as.data.frame() %>%
        select(-any_of("strand")) %>%
        gintervals.neighbors1(promoters) %>%
        filter(dist == 0) %>%
        select(chrom, start, end, peak_name, chrom1, start1, end1, strand, !!id_field) %>%
        rename(promoter_name := !!id_field)

    new_mat <- tgs_matrix_tapply(t(atac@mat[peak_map$peak_name, ]), peak_map$promoter_name, sum, na.rm = TRUE)

    new_peaks <- peak_map %>%
        group_by(promoter_name) %>%
        summarise(chrom = chrom[1], start = min(start), end = min(end), strand = strand[1], .groups = "drop") %>%
        select(chrom:strand, peak_name = promoter_name)
    new_mat <- new_mat[new_peaks$peak_name, ]

    atac@peaks <- PeakIntervals(new_peaks, atac@genome)
    atac@mat <- new_mat[atac@peaks$peak_name, ]
    atac@promoters <- TRUE

    if (nrow(atac@ignore_peaks) > 0) {
        cli_alert_warning("Removing previous {.field @ignore_peaks} and {.field @ignore_pmat}")
        atac@ignore_peaks <- subset(atac@peaks, subset = rep(FALSE, nrow(atac@peaks)))
        atac@ignore_pmat <- methods::as(matrix(0, nrow = 0, ncol = ncol(atac@mat)), "dgCMatrix")
    }

    cli_alert_success("Created a new {.class {class(atac)}} object with {.field {nrow(atac@peaks)}} promoters.")

    return(atac)
}

#' Calculate the cross-correlation between promoter accessibility and RNA expression
#'
#' @param atac_mc a McATACPeaks object with promoters (using \code{gen_promoter_features}) and RNA expression (using \code{add_mc_rna})
#' @param match_genes calculate correlation only between genes and promoters which have both accesability and RNA expression.
#' Matching is done based on name. Default: FALSE
#'
#' @inheritParams tgstat::tgs_cor
#' @inheritParams get_rna_egc
#'
#' @return a correlation matrix where rows are promoters and columns are genes
#'
#' @export
calc_prom_rna_cor <- function(atac_mc, genes = NULL, rm_zeros = TRUE, match_genes = FALSE, spearman = FALSE, pairwise.complete.obs = TRUE) {
    assert_atac_object(atac_mc, "McATACPeaks")
    if (!atac_mc@promoters) {
        cli_abort("{.val {atac_mc}} does not contain promoters.")
    }
    rna_mat <- get_rna_egc(atac_mc, genes = genes, rm_zeros = rm_zeros)

    atac_mat <- atac_mc@mat
    if (match_genes) {
        both_genes <- intersect(rownames(atac_mat), rownames(rna_mat))
        atac_mat <- atac_mat[both_genes, ]
        rna_mat <- rna_mat[both_genes, ]
        cli_alert("Calculating cross-correlation between {.field {length(both_genes)}} genes which have both ATAC and RNA data.")
    } else {
        cli_alert("Calculating cross-correlation between {.field {nrow(atac_mat)}} promoters and {.field {nrow(rna_mat)}} genes.")
    }

    cm <- tgs_cor(t(atac_mat), t(rna_mat), pairwise.complete.obs = pairwise.complete.obs, spearman = spearman)

    return(cm)
}
