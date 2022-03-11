
#' Annotate ATAC peaks
#'
#' @param atac a McATAC or ScATAC object
#'
#' @return the original \code{atac} object, with the peaks annotated using \code{annotate_intervals}
#'
#' @export
annotate_peaks <- function(atac) {
    if (!has_name("peaks", atac)) {
        cli_abort("{.field atac} doesn't have a field called {.field peaks}. You can annotate intervals directly using the {.code annotate_intervals} function.")
    }

    atac$peaks <- annotate_intervals(atac$peaks)

    return(atac)
}

#' Annotate an intervals set
#' ATAC peaks are classified into either 'promoter', 'intronic', 'exonic', 'ig_proximal' (intergenic-proximal), 'ig_distal' (intergenic-distal) or 'desert'.
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
#' @return an intervals set with the following annotations: 'promoter', 'intronic', 'exonic', 'ig_proximal' (intergenic-proximal), 'ig_distal' (intergenic-distal).
#'
#' @export
annotate_intervals <- function(intervals, genome, min_proximal = 1e+03, max_proximal = 2e+04, 
                min_distal = 2e+04, max_distal = 1e+06, exonic_peak_dist = 0) {
    # TODO: update peak annotation data frame with neighbor gene, distance, more (exons, enhancers etc.). The annotation might be an optional paramter.

    # TODO: check whether a genome is loaded
    # TODO: check whether 'tss' and 'exons' interval sets are already loaded
    
    gsetroot('/home/aviezerl/mm10')
    tss = gintervals.load('intervs.global.tss')
    exons = gintervals.load('intervs.global.exon')
    
    if (length(grep('intervalID', colnames(intervals))) == 0) {intervals$intervalID = 1:nrow(intervals)}
    nei_peak_prom = gintervals.neighbors(intervals, tss, mindist = -min_proximal, maxdist = min_proximal)
    nei_peak_exon = gintervals.neighbors(intervals, exons, mindist = -exonic_peak_dist, maxdist = exonic_peak_dist)
    nei_peak_gb = gintervals.neighbors(intervals, gene_body_df, maxdist = 0, mindist = 0)
    prom_peaks = nei_peak_prom$intervalID
    exon_peaks = nei_peak_exon$intervalID[!(nei_peak_exon$intervalID %in% prom_peaks)]
    intron_peaks = nei_peak_gb$intervalID[!(nei_peak_gb$intervalID %in% union(exon_peaks, prom_peaks))]
    intID_left = intervals$intervalID[!(intervals$intervalID %in% union(prom_peaks, union(exon_peaks, intron_peaks)))]
    nei_peak_tss_prox = gintervals.neighbors(intervals[intervals$intervalID %in% intID_left,], tss, mindist = min_proximal, maxdist = max_proximal)
    nei_peak_tss_prox_neg = gintervals.neighbors(intervals[intervals$intervalID %in% intID_left,], tss, maxdist = -min_proximal, mindist = -max_proximal)
    nei_peak_prox_all = dplyr::anti_join(unique(rbind(nei_peak_tss_prox[,1:4], nei_peak_tss_prox_neg[,1:4])), 
        intervals[union(prom_peaks, union(exon_peaks, intron_peaks)),1:4], by = c('chrom', 'start', 'end', 'intervalID'))
    ig_prox_peaks = nei_peak_prox_all$intervalID
    intID_left = intervals$intervalID[!(intervals$intervalID %in% unique(c(prom_peaks, exon_peaks, intron_peaks, ig_prox_peaks)))]
    nei_peak_dist = gintervals.neighbors(intervals[intervals$intervalID %in% intID_left,], tss, mindist = min_distal, maxdist = max_distal)
    nei_peak_dist_neg = gintervals.neighbors(intervals[intervals$intervalID %in% intID_left,], tss, maxdist = -min_distal, mindist = -max_distal)
    nei_peak_dist_all = dplyr::anti_join(unique(rbind(nei_peak_dist[,1:4], nei_peak_dist_neg[,1:4])), 
        intervals[union(ig_prox_peaks, union(prom_peaks, union(exon_peaks, intron_peaks))),1:4], by = c('chrom', 'start', 'end', 'intervalID'))
    ig_dist_peaks = nei_peak_dist_all$intervalID
    desert_peaks = intervals$intervalID[!(intervals$intervalID %in% 
                                    unique(c(prom_peaks, exon_peaks, intron_peaks, ig_prox_peaks, ig_dist_peaks)))]
    res = c(setNames(rep('promoter', length(prom_peaks)), prom_peaks),
            setNames(rep('exonic', length(exon_peaks)), exon_peaks),
            setNames(rep('intronic', length(intron_peaks)), intron_peaks),
            setNames(rep('ig_proximal', length(ig_prox_peaks)), ig_prox_peaks),
            setNames(rep('ig_distal', length(ig_dist_peaks)), ig_dist_peaks),
            setNames(rep('desert', length(desert_peaks)), desert_peaks)
    )
    intervals$peak_annot = res[order(as.numeric(names(res)))]
    return(intervals)
}


#' Generate dataframe of approximate gene bodies (edges of first and last exons)
#'
#' @param tss tss interval set built-in to misha
#' @param exons exon interval set built-in to misha
#'
#' @return gene_body_df a misha intervals set approximating starts and ends of genes
#'
#' @export
annotate_peaks <- function(tss, exons) {
    kgid_both = intersect(unique(exons$kgID), unique(tss$kgID))
    exons_filt = exons[exons$kgID %in% kgid_both,]
    genes_start = tapply(exons_filt$start, exons_filt$kgID, min)
    genes_end = tapply(exons_filt$end, exons_filt$kgID, max)
    gene_body_df = as.data.frame(do.call('cbind', 
                                        list('chrom' = as.character(exons_filt$chrom)[match(names(genes_start), exons_filt$kgID)], 
                                            'start' = genes_start,
                                            'end' = genes_end,
                                            'strand' = exons_filt$strand[match(names(genes_start), exons_filt$kgID)],
                                            'mRNA' = exons_filt$mRNA[match(names(genes_start), exons_filt$kgID)],
                                            'geneSymbol' = exons_filt$geneSymbol[match(names(genes_start), exons_filt$kgID)])
                                        )) 
    gene_body_df[,2:4] = apply(gene_body_df[,2:4], 2, as.numeric) 
    gene_body_df = gene_body_df[with(gene_body_df, order(chrom, start)),]
    return(gene_body_df)
}