#' Annotate ATAC peaks
#'
#' @param atac a McPeaks or ScPeaks object
#'
#' @return the original \code{atac} object, with the peaks annotated using \code{annotate_intervals}
#' @examples
#' \dontrun{
#' pbmc_atac_mc <- annotate_peaks(pbmc_atac_mc)
#' }
#' @export
annotate_peaks <- function(atac) {
    if (!methods::is(atac, "ATACPeaks")) {
        cli_abort("{.field atac} is not an ScPeaks or McPeaks object. You can annotate intervals directly using the {.code annotate_intervals} function.")
    }

    atac@peaks <- annotate_intervals(atac@peaks, atac@genome)
    return(atac)
}

#' Annotate an intervals set
#'
#' @description ATAC peaks are classified into either 'promoter', 'intronic', 'exonic', 'ig_proximal' (intergenic-proximal), 'ig_distal' (intergenic-distal) or 'desert'. \cr
#' Promoter peaks are those less than \code{min_proximal} away from a TSS. \cr
#' Intronic peaks are within gene bodies, but not in an exon. \cr
#' Exonic peaks overlap exons. \cr
#' Intergenic-proximal peaks are less than \code{max_proximal} away from a TSS, but not within a gene body. \cr
#' Intergenic-distal peaks are more than \code{min_distal} away from some TSS, and less than \code{max-distal} from any TSS, and not within a gene body. \cr
#' Desert peaks are at least \code{max_distal} from any known TSS (as defined by the \code{tss} interval set).
#'
#' @param intervals the intervals set to annotate
#' @param genome which genome/genomic database to use (e.g. 'mm10', 'hg19')
#' @param min_proximal the minimum distance from a TSS, under which a peak is considered to be a promoter peak, and beyond (or at) which it is a proximal enhancer.
#' @param max_proximal the minimum distance from a TSS, under which a peak is considered to be a proximal enhancer peak, and beyond (or at) which it is not.
#' @param min_distal the minimum distance from some TSS, only above which a peak is considered to be of a distal enhancer. Usually equal to \code{max_proximal}, unless there is good reason
#' @param max_distal the maximum distance from any TSS, under which a peak is considered to be of a distal enhancer, and above which it is a 'desert' peak.
#' @param exonic_peak_dist a parameter that allows proximity (without overlap) of peaks to exons that still classifies them as exonic peaks, as active transcription usually confers accessibility in the vicinity of exons, and this might confound a possible assignment of regulatory activity to an exon.
#' @param tss the TSS interval set or the name of the TSS interval set to use.
#' @param exons the exons interval set or the name of the exons interval set to use.
#'
#' @return intervals - the input intervals set with added fields {.code peak_type}, {.code closest_tss}, {.code closest_exon_gene} and {.code gene_body_gene}.
#' @examples
#' \dontrun{
#' my_intervals <- annotate_intervals(my_intervals, "mm10", min_proximal = 1e+03, max_proximal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 5e+2)
#' table(my_intervals$peak_type)
#' my_intervals[which(toupper(my_intervals$closest_tss) == "PCNA"), ]
#' }
#' @export
annotate_intervals <- function(intervals, genome,
                               min_proximal = 1e+03,
                               max_proximal = 2e+04,
                               min_distal = max_proximal,
                               max_distal = 1e+06,
                               exonic_peak_dist = 0,
                               tss = "intervs.global.tss",
                               exons = "intervs.global.exon") {
    if (missing(genome)) {
        cli_abort("Please Specify genome. Look for slot '@genome' in relevant McPeaks/ScPeaks object.")
    }
    misha.ext::gset_genome(genome)

    cn <- c("chrom", "start", "end", "peak_name")
    orig_class <- class(intervals)
    orig_fields <- colnames(intervals)

    gene_body_df <- get_gene_body_df(tss, exons)

    df <- as.data.frame(intervals) %>%
        gintervals.neighbors1(tss, fields = "geneSymbol") %>%
        rename(closest_tss = geneSymbol, tss_d = dist) %>%
        gintervals.neighbors1(exons, fields = "geneSymbol") %>%
        rename(closest_exon_gene = geneSymbol, exons_d = dist) %>%
        gintervals.neighbors1(gene_body_df, fields = "geneSymbol") %>%
        rename(gene_body_gene = geneSymbol, gene_body_d = dist)

    df$peak_type <- NA

    add_annotation <- function(df, column, mindist, maxdist, name) {
        df <- df %>%
            mutate(peak_type = ifelse(is.na(peak_type) & abs(!!sym(column)) <= maxdist & abs(!!sym(column)) >= mindist, name, peak_type))
    }

    res <- df %>%
        add_annotation("tss_d", "promoter", mindist = 0, maxdist = min_proximal) %>%
        add_annotation("exons_d", "exonic", mindist = exonic_peak_dist, maxdist = exonic_peak_dist) %>%
        add_annotation("gene_body_d", "intronic", mindist = 0, maxdist = 0) %>%
        add_annotation("tss_d", "ig_proximal", mindist = min_proximal, maxdist = max_proximal) %>%
        add_annotation("tss_d", "ig_distal", mindist = min_distal, maxdist = max_distal) %>%
        tidyr::replace_na(replace = list(peak_type = "desert"))

    class(intervals) <- orig_class
    intervals <- intervals %>%
        select(any_of(orig_fields)) %>%
        left_join(res %>% select(chrom, start, end, peak_type, closest_tss, closest_exon_gene, gene_body_gene), by = c("chrom", "start", "end"))

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
    if (is.character(tss)) {
        tss <- gintervals.load(tss)
    }
    if (is.character(exons)) {
        exons <- gintervals.load(exons)
    }
    kgid_both <- intersect(unique(exons$kgID), unique(tss$kgID))
    exons_filt <- exons[exons$kgID %in% kgid_both, ]
    genes_start <- tapply(exons_filt$start, exons_filt$kgID, min)
    genes_end <- tapply(exons_filt$end, exons_filt$kgID, max)
    gene_body_df <- data.frame(
        "chrom" = as.character(exons_filt$chrom)[match(names(genes_start), exons_filt$kgID)],
        "start" = genes_start,
        "end" = genes_end,
        "strand" = exons_filt$strand[match(names(genes_start), exons_filt$kgID)],
        "mRNA" = exons_filt$mRNA[match(names(genes_start), exons_filt$kgID)],
        "geneSymbol" = exons_filt$geneSymbol[match(names(genes_start), exons_filt$kgID)]
    )
    gene_body_df <- gene_body_df %>%
        mutate_at(vars(start, end, strand), as.numeric) %>%
        arrange(chrom, start)
    return(gene_body_df)
}

#' Name ATAC peaks from TAD names
#'
#' @description
#' The peaks are named in the following pattern "{TAD_name}_{distance_from_TAD_start_in_kbp}".
#' For the TAD naming scheme see: https://github.com/tanaylab/tadbris. \cr
#' When two peaks are within the same kbp of each other, they are appended a number to the name, according to
#' their order on the genome (except for the first peak), i.e.: "{TAD_name}_{distance_from_TAD_start_in_kbp}.{number}".
#' See also \code{make.unique}.
#'
#' @param atac a McPeaks or ScPeaks a PeakIntervals object or a misha intervals set
#' @return the original \code{atac} object, with the field \code{peak_name} replaced by the new TAD-referenced name
#' @examples
#' \dontrun{
#' atac_mc <- name_enhancers(atac_mc)
#' my_peaks <- name_enhancers(my_peaks)
#' }
#' @export
name_enhancers <- function(atac, tad_names = gintervals.load("intervs.global.tad_names")) {
    cl <- class(atac)
    if (cl[[1]] %in% c("ScPeaks", "McPeaks")) {
        peaks <- atac@peaks
    } else if ("PeakIntervals" %in% cl) {
        peaks <- atac
    } else if ("data.frame" %in% cl && all(c("chrom", "start", "end") %in% colnames(atac))) {
        peaks <- atac
    } else {
        cli_abort("Class of {.var atac} is not recognized (should be either ScPeaks, McPeaks, PeakIntervals object or misha intervals).")
    }
    if (has_name(peaks, "peak_original_order")) {
        cli_abort("{.field peaks} shouldn't have a field named 'peak_original_order'. Sorry about that.")
    }
    peaks <- as.data.frame(peaks) %>%
        misha.ext::gintervals.neighbors1(tad_names) %>%
        mutate(
            dist_diff = start - start1,
            peak_name = glue("{tad_name}_{kb}", kb = round(dist_diff, -3) / 1e+3),
            peak_original_order = 1:n()
        ) %>%
        arrange(chrom, start, end) %>% # we want to make sure the naming is stable regardless of the order of the peaks
        mutate(peak_name = make.unique(peak_name)) %>%
        arrange(peak_original_order) %>% # restore the original order
        select(-(chrom1:dist_diff), -peak_original_order) %>%
        PeakIntervals()

    if (cl[[1]] %in% c("ScPeaks", "McPeaks")) {
        atac@peaks <- peaks
    } else if (cl[[1]] == "PeakIntervals") {
        atac <- peaks
    } else if ("data.frame" %in% cl && all(c("chrom", "start", "end") %in% colnames(atac))) {
        atac <- peaks
    }

    return(atac)
}

#' Generate dataframe of approximate gene bodies (edges of first and last exons) - backend function
#'
#' @param intervals query interval set
#' @param feature_set reference feature set to query
#' @param maxdist (optional) maximal absolute distance to search for feature
#'
#' @return gene_body_df a misha intervals set approximating starts and ends of genes
#'
#' @noRd
get_nearest_feature_instance <- function(intervals, feature_set, maxdist = 1e+6) {
    nei_peak_feat <- gintervals.neighbors(intervals, feature_set, mindist = -maxdist, maxdist = maxdist)
    closest_feat <- deframe(nei_peak_feat[, c("peak_name", "geneSymbol")])
    closest_feat[intervals$peak_name[!(intervals$peak_name %in% names(closest_feat))]] <- NA
    closest_feat <- closest_feat[order(match(names(closest_feat), intervals$peak_name))]
    return(closest_feat)
}
