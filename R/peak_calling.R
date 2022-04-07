#' Import marginal of ATAC counts to a misha track
#'
#'
#' @param file Name of the 'bigwig' file to import, such as the 'atac_cut_sites.bigwig' from the 10x pipeline.
#' For more details, see:
#' https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/
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
import_atac_marginal <- function(file, track, description, genome, binsize = 200, overwrite = FALSE, wig_temp_dir = tempdir()) {
    gset_genome(genome)
    if (gtrack.exists(track)) {
        if (overwrite) {
            cli_alert_warning("Removing previous track {.val track}")
            gtrack.rm(track, force = TRUE)
        } else {
            cli_abort("Track '{.var {track}}' already exists. Use {.code overwrite = TRUE} to overwrite it.")
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
    prefix <- gsub(".bigwig$", "", basename(bigwig_file))
    for (chrom in chroms) {
        out_file <- paste0(wig_temp_dir, "/", prefix, "_", chrom, ".wig")
        cmd <- paste0(bin, " -chrom=", chrom, " ", bigwig_file, " ", out_file)
        system(cmd)
        withr::defer(unlink(out_file))
        cli::cli_progress_update()
    }
    cmd <- paste0("cat ", wig_temp_dir, "/*.wig > ", wig_file)
    system(cmd)
    cli::cli_progress_update()
    cli::cli_progress_done()
    invisible(wig_file)
}
