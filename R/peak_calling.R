#' Import marginal of ATAC counts to a misha track
#'
#'
#' @param file Name of the 'bigwig' file to import, such as the 'atac_cut_sites.bigwig' from the 10x pipeline.
#' For more details, see:
#' https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bigwig
#' @param genome Genome name, such as 'hg19' or 'mm10'.
#' @param overwrite Overwrite the existing track if it exists.
#' @param wig_temp_dir Temporary directory to store the intermediate wig files.
#'
#' @examples
#' \dontrun{
#' import_atac_marginal("./atac_cut_sites.bigwig", "pbmc_atac.marginal", "Marignal counts from PBMC ATAC", genome = "hg38")
#' }
#'
#' @inheritParams misha::gtrack.import
#' @export
import_atac_marginal <- function(file, track, description, genome, binsize = 20, overwrite = FALSE, wig_temp_dir = tempdir()) {
    gset_genome(genome)
    if (gtrack.exists(track)) {
        if (overwrite) {
            cli_alert_warning("Removing previous track {.val {track}}")
            gtrack.rm(track, force = TRUE)
            gdb.reload()
            assert_that(!gtrack.exists(track))
        } else {
            cli_abort("Track {.val {track}} already exists. Use {.code overwrite = TRUE} to overwrite it.")
        }
    }
    cli_li("Converting {.field bigwig} file to {.field wig}")
    wig_file <- paste0(tempfile(), ".wig")
    bigwig_to_wig(file, wig_file, genome, wig_temp_dir = wig_temp_dir)

    cli_li("Creating {.field misha} track")
    misha.ext::gtrack.create_dirs(track, showWarnings = FALSE)
    gtrack.import(track, description, wig_file, binsize = binsize)
    cli_alert_success("Successfully imported marginal ATAC counts to {.val {track}} from {.file {file}}")
    cli_text("To generate a new set of peaks, please run {.field call_peaks(\"{track}\")}")
}

#' Convert bigwig file to wig file
#'
#' @description Convert a bigwig file to wig file while keeping only chromosomes that exists in the misha database.
#'
#' @param bigwig_file Name of the 'bigwig' file to convert, such as the 'atac_cut_sites.bigwig' from the 10x pipeline.
#' @param wig_file Name of the 'wig' file to create.
#' @param genome Genome name, such as 'hg19' or 'mm10'.
#' @param wig_temp_dir Temporary directory to store the intermediate wig files.
#'
#' @return name of the wig_file (invisibly)
#'
#' @noRd
bigwig_to_wig <- function(bigwig_file, wig_file, genome, wig_temp_dir = tempdir()) {
    gset_genome(genome)
    bin <- system.file("exec/bigWigToWig", package = "misha")
    chroms <- gintervals.all()$chrom
    cli::cli_progress_bar("Converting", total = length(chroms) + 1)
    prefix <- gsub(".bigwig$", "", basename(tempfile()))
    for (chrom in chroms) {
        out_file <- paste0(wig_temp_dir, "/", prefix, "_", chrom, ".wig")
        cmd <- paste0(bin, " -chrom=", chrom, " ", bigwig_file, " ", out_file)
        system(cmd)
        cli::cli_progress_update()
    }
    cmd <- paste0("cat ", wig_temp_dir, glue("/{prefix}*.wig > "), wig_file)
    system(cmd)
    cli::cli_progress_update()
    cli::cli_progress_done()
    system(glue("rm -f {wig_temp_dir}/{prefix}*.wig"))
    invisible(wig_file)
}

#' Call peaks from ATAC marginals track
#'
#' @description Peaks are called by screening for genomic regions with a number of UMIs above a quantile of the marginal ATAC counts,
#' and above \code{min_umis}. Peaks that are longer than \code{max_peak_size} would be splitted equally into smaller peaks.
#'
#' @param marginal_track Name of the 'misha' track to call peaks from. You can create it using {.code import_atac_marginal}.
#' @param quantile_thresh Quantile of the marginal track above which peaks are called.
#' @param min_umis Minimal number of UMIs to call a peak.
#' @param max_peak_size Maximal size of a peak. Peaks above this size would be splitted equally into smaller peaks.
#' @param min_peak_size Minimal size of a peak.
#' @param genome Genome name, such as 'hg19' or 'mm10'. If NULL - the current misha database is used.
#'
#' @return an intervals set with the called peaks
#'
#' @examples
#' \dontrun{
#' peaks <- call_peaks("pbmc_atac.marginal", quantile_thresh = 0.95, min_umis = 8, max_peak_size = 600, genome = "hg38")
#' }
#'
#' @export
call_peaks <- function(marginal_track, quantile_thresh = 0.9, min_umis = 8, max_peak_size = 600, min_peak_size = 50, genome = NULL) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }
    thresh <- max(gquantiles(marginal_track, quantile_thresh), min_umis)

    cli::cli_alert_info("Coverage threshold: {.val {round(thresh, digits=3)}}")


    df <- gscreen(glue("{marginal_track} >= thresh"), intervals = gintervals.all())

    # split peaks that are longer than max_peak_size

    return(df)
}

plot_marginal_coverage <- function(marginal_track, scope, peaks = NULL, show_thresh = TRUE, quantile_thresh = 0.9, min_umis = 8, genome = NULL, log_scale = TRUE) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }

    thresh <- NULL
    if (show_thresh) {
        thresh <- max(gquantiles(marginal_track, quantile_thresh), min_umis)
    }

    ggdata <- gextract(marginal_track, scope, colnames = "counts")

    if (!is.null(peaks)) {
        peaks_f <- peaks %>%
            gintervals.neighbors1(scope) %>%
            filter(dist == 0) %>%
            mutate(start = pmax(scope$start, start), end = pmin(scope$end, end))
    }

    p <- ggdata %>%
        ggplot(aes(x = start, y = counts)) +
        geom_col() +
        xlim(c(scope$start, scope$end))

    if (!is.null(peaks)) {
        p <- p + geom_rect(data = peaks_f, inherit.aes = FALSE, aes(xmin = start, xmax = end), ymin = 0, ymax = max(ggdata$counts), alpha = 0.1)
    }

    if (!is.null(thresh)) {
        p <- p + geom_hline(yintercept = thresh, color = "blue", linetype = "dashed")
    }

    if (log_scale) {
        p <- p +
            scale_y_log10()
    }


    return(p)
}
