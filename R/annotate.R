
#' Annotate ATAC peaks
#'
#' @param atac a McATAC or ScATAC object
#'
#' @return the original \code{atac} object, with the peaks annotated using \code{annotate_intervals}
#' @examples
#' \dontrun{
#' pbmc_atac_mc <- annotate_peaks(pbmc_atac_mc)
#' }
#' @export
annotate_peaks <- function(atac) {
    if (!methods::is(atac, "ATAC")) {
        cli_abort("{.field atac} is not an ScATAC or McATAC object. You can annotate intervals directly using the {.code annotate_intervals} function.")
    }

    atac@peaks <- annotate_intervals(atac@peaks, atac@genome)
    return(atac)
}

#' Annotate an intervals set
#'
#' @description ATAC peaks are classified into either 'promoter', 'intronic', 'exonic', 'ig_proximal' (intergenic-proximal), 'ig_distal' (intergenic-distal) or 'desert'.
#' Promoter peaks are those less than \code{min_proximal} away from a TSS.
#' Intronic peaks are within gene bodies, but not in an exon.
#' Exonic peaks overlap exons.
#' Intergenic-proximal peaks are less than \code{max_proximal} away from a TSS, but not within a gene body
#' Intergenic-distal peaks are more than \code{min_distal} away from some TSS, and less than \code{max-distal} from any TSS, and not within a gene body.
#' Desert peaks are at least \code{max_distal} from any known TSS (as defined by the \code{tss} interval set).
#'
#' @param intervals the intervals set to annotate
#' @param genome which genome/genomic database to use (e.g. 'mm10', 'hg19')
#' @param min_proximal the minimum distance from a TSS, under which a peak is considered to be a promoter peak, and beyond (or at) which it is a proximal enhancer.
#' @param max_proximal the minimum distance from a TSS, under which a peak is considered to be a proximal enhancer peak, and beyond (or at) which it is not.
#' @param min_distal the minimum distance from some TSS, only above which a peak is considered to be of a distal enhancer. Usually equal to \code{max_proximal}, unless there is good reason
#' @param max_distal the maximum distance from any TSS, under which a peak is considered to be of a distal enhancer, and above which it is a 'desert' peak.
#' @param exonic_peak_dist a parameter that allows proximity (without overlap) of peaks to exons that still classifies them as exonic peaks, as active transcription usually confers accessibility in the vicinity of exons, and this might confound a possible assignment of regulatory activity to an exon.
#'
#' @return intervals - the input intervals set with added fields {.code peak_type} and {.code closest_tss}
#' @examples
#' \dontrun{
#' my_intervals <- annotate_intervals(my_intervals, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
#' table(my_intervals$peak_type)
#' my_intervals[which(toupper(my_intervals$closest_tss) == "PCNA"), ]
#' }
#' @export
annotate_intervals <- function(intervals, genome,
                               min_proximal = 1e+03, max_proximal = 2e+04,
                               min_distal = 2e+04, max_distal = 1e+06,
                               exonic_peak_dist = 0,
                               tss = gintervals.load("intervs.global.tss"),
                               exons = gintervals.load("intervs.global.exon")) {
    if (missing(genome)) {
        cli_abort("Please Specify genome. Look for slot '@genome' in relevant McATAC/ScATAC object.")
    }
    misha.ext::gset_genome(genome)
    cn <- c("chrom", "start", "end", "peak_name")
    orig_class <- class(intervals)
    orig_fields <- colnames(intervals)
    intervals <- as.data.frame(intervals)

    if (!has_name(intervals, "intervalID")) {
        intervals <- intervals %>% mutate(intervalID = 1:n())
    }
    cli_alert_info("Finding peaks in promoters...")
    nei_peak_prom <- gintervals.neighbors(intervals, tss, mindist = -min_proximal, maxdist = min_proximal)
    cli_alert_info("Finding peaks in exons...")
    nei_peak_exon <- gintervals.neighbors(intervals, exons, mindist = -exonic_peak_dist, maxdist = exonic_peak_dist)
    cli_alert_info("Finding peaks in introns...")
    gene_body_df <- get_gene_body_df(tss, exons)
    nei_peak_gb <- gintervals.neighbors(intervals, gene_body_df, maxdist = 0, mindist = 0)

    nei_peak_prom_all <- gintervals.neighbors(intervals, tss, mindist = -max_distal, maxdist = max_distal)
    closest_tss <- deframe(nei_peak_prom_all[, c("peak_name", "geneSymbol")])
    closest_tss[!(intervals$peak_name %in% names(closest_tss))] <- NA
    closest_tss <- closest_tss[order(match(names(closest_tss), intervals$peak_name))]

    prom_peaks <- nei_peak_prom$peak_name
    exon_peaks <- nei_peak_exon$peak_name[!(nei_peak_exon$peak_name %in% prom_peaks)]
    intron_peaks <- nei_peak_gb$peak_name[!(nei_peak_gb$peak_name %in% union(exon_peaks, prom_peaks))]
    intID_left <- intervals$peak_name[!(intervals$peak_name %in% union(prom_peaks, union(exon_peaks, intron_peaks)))]

    cli_alert_info("Finding proximal intergenic peaks...")
    nei_peak_tss_prox <- gintervals.neighbors(intervals[intervals$peak_name %in% intID_left, ], tss, mindist = min_proximal, maxdist = max_proximal)
    nei_peak_tss_prox_neg <- gintervals.neighbors(intervals[intervals$peak_name %in% intID_left, ], tss, maxdist = -min_proximal, mindist = -max_proximal)
    nei_peak_prox_all <- anti_join(unique(rbind(nei_peak_tss_prox[, cn], nei_peak_tss_prox_neg[, cn])),
        intervals[union(prom_peaks, union(exon_peaks, intron_peaks)), cn],
        by = c("chrom", "start", "end")
    )
    ig_prox_peaks <- nei_peak_prox_all$peak_name
    intID_left <- intervals$peak_name[!(intervals$peak_name %in% unique(c(prom_peaks, exon_peaks, intron_peaks, ig_prox_peaks)))]
    cli_alert_info("Finding distal intergenic peaks...")
    nei_peak_dist <- gintervals.neighbors(intervals[intervals$peak_name %in% intID_left, ], tss, mindist = min_distal, maxdist = max_distal)
    nei_peak_dist_neg <- gintervals.neighbors(intervals[intervals$peak_name %in% intID_left, ], tss, maxdist = -min_distal, mindist = -max_distal)
    nei_peak_dist_all <- anti_join(unique(rbind(nei_peak_dist[, cn], nei_peak_dist_neg[, cn])),
        intervals[union(ig_prox_peaks, union(prom_peaks, union(exon_peaks, intron_peaks))), cn],
        by = c("chrom", "start", "end")
    )
    ig_dist_peaks <- nei_peak_dist_all$peak_name
    desert_peaks <- intervals$peak_name[!(intervals$peak_name %in%
        unique(c(prom_peaks, exon_peaks, intron_peaks, ig_prox_peaks, ig_dist_peaks)))]

    res <- c(
        setNames(rep("promoter", length(prom_peaks)), prom_peaks),
        setNames(rep("exonic", length(exon_peaks)), exon_peaks),
        setNames(rep("intronic", length(intron_peaks)), intron_peaks),
        setNames(rep("ig_proximal", length(ig_prox_peaks)), ig_prox_peaks),
        setNames(rep("ig_distal", length(ig_dist_peaks)), ig_dist_peaks),
        setNames(rep("desert", length(desert_peaks)), desert_peaks)
    )

    class(intervals) <- orig_class
    intervals <- intervals %>%
        select(any_of(orig_fields)) %>%
        mutate(
            peak_annot = res[order(match(names(res), intervals$peak_name))],
            closest_tss = closest_tss
        )
    return(intervals)
}


#' Generate dataframe of approximate gene bodies (edges of first and last exons) - backend function
#'
#' @param tss tss interval set built-in to misha
#' @param exons exon interval set built-in to misha
#'
#' @return gene_body_df a misha intervals set approximating starts and ends of genes
#'
#' @noRd
get_gene_body_df <- function(tss, exons) {
    kgid_both <- intersect(unique(exons$kgID), unique(tss$kgID))
    exons_filt <- exons[exons$kgID %in% kgid_both, ]
    genes_start <- tapply(exons_filt$start, exons_filt$kgID, min)
    genes_end <- tapply(exons_filt$end, exons_filt$kgID, max)
    gene_body_df <- as.data.frame(do.call(
        "cbind",
        list(
            "chrom" = as.character(exons_filt$chrom)[match(names(genes_start), exons_filt$kgID)],
            "start" = genes_start,
            "end" = genes_end,
            "strand" = exons_filt$strand[match(names(genes_start), exons_filt$kgID)],
            "mRNA" = exons_filt$mRNA[match(names(genes_start), exons_filt$kgID)],
            "geneSymbol" = exons_filt$geneSymbol[match(names(genes_start), exons_filt$kgID)]
        )
    ))
    gene_body_df[, 2:4] <- apply(gene_body_df[, 2:4], 2, as.numeric)
    gene_body_df <- gene_body_df[with(gene_body_df, order(chrom, start)), ]
    return(gene_body_df)
}
