#' Import marginal of ATAC counts to a misha track
#'
#'
#' @param file Name of the 'bigwig' file to import, such as the 'atac_cut_sites.bigwig' from the 10x pipeline.
#' For more details, see:
#' https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bigwig
#' @param track Name of the track to create.
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


#' Calculate the threshold for calling peaks
#'
#' @description Threshold is calculated by taking the maximum between the \code{quantile_thresh} quantile of the ATAC marginal and \code{min_umis}.
#'
#' @param marginal_track Name of the 'misha' track to call peaks from. You can create it using \code{import_atac_marginal}.
#' @param quantile_thresh Quantile threshold to use.
#' @param min_umis Minimum number of UMIs to use.
#' @param genome Genome name, such as 'hg19' or 'mm10'. If NULL - the current genome is used.
#' @param seed random seed for reproducibility (\code{misha::gquantiles} sometimes samples the data.
#'
#' @return The number of UMIs to use for calling peaks.
#'
#' @examples
#' \dontrun{
#' get_quantile_cov_thresh("pbmc_atac.marginal", quantile_thresh = 0.9, min_umis = 8, genome = "hg38")
#' }
#'
#' @export
get_quantile_cov_thresh <- function(marginal_track, quantile_thresh, min_umis, genome = NULL, seed = 60427) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }
    withr::with_seed(seed, thresh <- max(gquantiles(marginal_track, quantile_thresh), min_umis))
    return(thresh)
}

#' Call peaks from ATAC marginals track
#'
#' @description Peaks are called by screening for genomic regions with a number of UMIs above a quantile of the marginal ATAC counts,
#' and above \code{min_umis}. Peaks that are longer than \code{max_peak_size} would be splitted equally into smaller peaks.
#'
#' @param marginal_track Name of the 'misha' track to call peaks from. You can create it using \code{import_atac_marginal}.
#' @param quantile_thresh Quantile of the marginal track above which peaks are called.
#' @param min_umis Minimal number of UMIs to call a peak.
#' @param genome Genome name, such as 'hg19' or 'mm10'. If NULL - the current misha database is used.
#' @param split_peaks Split peaks that are longer than \code{max_peak_size} into smaller peaks using \code{split_long_peaks}.
#' @param seed random seed for reproducibility (\code{misha::gquantiles} sometimes samples the data)
#' @inheritParams split_long_peaks
#'
#' @return an intervals set with the called peaks
#'
#' @examples
#' \dontrun{
#' peaks <- call_peaks("pbmc_atac.marginal", quantile_thresh = 0.9, min_umis = 8, max_peak_size = 600, genome = "hg38")
#' }
#'
#' @export
call_peaks <- function(marginal_track, quantile_thresh = 0.9, min_umis = 8, split_peaks = TRUE, target_size = 500, max_peak_size = 1e3, very_long = 5e3, min_peak_size = 200, genome = NULL, seed = 60427) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }
    thresh <- get_quantile_cov_thresh(marginal_track, quantile_thresh, min_umis, genome = genome, seed = seed)
    cli::cli_alert_info("Coverage threshold: {.val {round(thresh, digits=3)}}")


    df <- gscreen(glue("{marginal_track} >= thresh"), intervals = gintervals.all())

    if (split_peaks) {
        df <- split_long_peaks(marginal_track, peaks = df, target_size = target_size, max_peak_size = max_peak_size, very_long = very_long, min_peak_size = min_peak_size)
    }

    return(df)
}


#' Plot marginal coverage of an interval
#'
#' @description Plot the marginal coverage of an interval, with the peaks marked (optionally).
#'
#' @param marginal_track Name of the 'misha' track to plot. You can create it using \code{import_atac_marginal}.
#' @param interval An interval to plot.
#' @param peaks An intervals set with the peaks to mark, e.g. output of \code{call_peaks}.
#' @param expand Expand the plotting area by this number of bp.
#' @param show_threshold Show the coverage threshold as a dashed line.
#' @param quantile_thresh,min_umis paramters needed to calculate the threshold. See \code{call_peaks}.
#' @param genome Genome name, such as 'hg19' or 'mm10'. If NULL - the current misha database is used.
#' @param thresh Threshold to use. If NULL - the threshold is calculated using \code{get_quantile_cov_thresh}.
#' @param log_scale Use log scale for the y axis.
#'
#' @examples
#' \dontrun{
#' peaks_raw <- call_peaks("pbmc_atac.marginal", split_peaks = FALSE, quantile_thresh = 0.9, min_umis = 8, max_peak_size = 600, genome = "hg38")
#' peaks_split <- call_peaks("pbmc_atac.marginal", split_peaks = TRUE, target_size = 500, quantile_thresh = 0.9, min_umis = 8, max_peak_size = 600, genome = "hg38")
#' plot_marginal_coverage("pbmc_atac.marginal", interval = peaks_raw[967, ], peaks = peaks_split, expand = 1000, show_thresh = TRUE, quantile_thresh = 0.9, min_umis = 8, genome = "hg38")
#'
#' # cache the threshold in order to plot multiple intervals
#' thresh <- get_quantile_cov_thresh("pbmc_atac.marginal", 0.9, 8, genome = "hg38", seed = 60427)
#' plot_marginal_coverage("pbmc_atac.marginal", interval = peaks_raw[967, ], peaks = peaks_split, expand = 1000, show_thresh = TRUE, thresh = thresh, genome = "hg38")
#' plot_marginal_coverage("pbmc_atac.marginal", interval = peaks_raw[900, ], peaks = peaks_split, expand = 1000, show_thresh = TRUE, thresh = thresh, genome = "hg38")
#' }
#'
#' @export
plot_marginal_coverage <- function(marginal_track, interval, peaks = NULL, expand = 1e3, show_thresh = TRUE, quantile_thresh = 0.9, min_umis = 8, genome = NULL, seed = 60427, thresh = get_quantile_cov_thresh(marginal_track, quantile_thresh, min_umis, genome = genome, seed = seed), log_scale = TRUE) {
    if (!is.null(genome)) {
        gset_genome(genome)
    }

    scope <- interval %>%
        mutate(start = start - expand, end = end + expand) %>%
        as.data.frame() %>%
        gintervals.force_range()

    if (nrow(scope) > 1) {
        cli_alert_warning("More than one interval in the scope, taking the first one")
        scope <- scope %>% slice(1)
    }

    ggdata <- gextract(marginal_track, scope, colnames = "counts")
    if (!is.null(peaks)) {
        peaks_f <- peaks %>%
            mutate(intervalID = factor(1:n())) %>%
            as.data.frame() %>%
            gintervals.neighbors1(scope) %>%
            filter(dist == 0) %>%
            mutate(start = pmax(start, scope$start), end = pmin(end, scope$end))
    }

    p <- ggdata %>%
        ggplot(aes(x = start, y = counts)) +
        geom_col() +
        xlim(c(scope$start, scope$end))

    if (!is.null(peaks)) {
        p <- p + geom_rect(data = peaks_f, inherit.aes = FALSE, aes(xmin = start, xmax = end, fill = intervalID), ymin = 0, ymax = max(ggdata$counts, na.rm = TRUE), alpha = 0.4) + chameleon::scale_fill_chameleon() + guides(fill = "none")
    }

    if (show_thresh) {
        p <- p + geom_hline(yintercept = thresh, color = "blue", linetype = "dashed")
    }

    if (log_scale) {
        p <- p +
            scale_y_log10()
    }

    p <- p + ggtitle(glue("{scope$chrom}: {scope$start}-{scope$end}"))

    return(p)
}

split_peaks_arbitrarily <- function(peaks, max_peak_size) {
    peaks <- peaks %>% select(chrom, start, end)
    bad_peaks <- peaks %>%
        filter(end - start > max_peak_size)
    splitted_peaks <- giterator.intervals(intervals = bad_peaks, iterator = max_peak_size)
    peaks_f <- peaks %>%
        anti_join(bad_peaks, by = c("chrom", "start", "end")) %>%
        bind_rows(splitted_peaks) %>%
        arrange(chrom, start, end)
    return(peaks_f)
}

#' Split long peaks into smaller peaks
#'
#' @description Peaks are splitted into smaller peaks if they are longer than \code{target_size}. The splitting is done by first checking
#' if the peak has a length above \code{very_long}. If it does, the peak is splitted arbitrarily into chunks of size that is approximately \code{very_long}. \cr
#' Then, the chunks are splitted into smaller chunks of approximately size \code{target_size}. This splitting is done by first finding the
#' best offset to start from and then splitting the peaks into chunks of size \code{target_size} starting from that offset. \cr
#' Detection of the best offset is done by substracting the mean of each interval from each coverage, removing values which became
#' zero, and then correlating thr marginal coverage with simulated 'triangle' peaks starting at different offsets.
#'
#'
#' @param marginal_track Name of the 'misha' track with the marginal coverage. You can create it using \code{import_atac_marginal}.
#' @param peaks An intervals set with the peaks to split.
#' @param target_size The target size of peaks.
#' @param max_peak_size Peaks above this size would be splitted into smaller peaks.
#' @param very_long Peaks above this size would be splitted arbitrarily into smaller peaks before fitting the best offset.
#' @param min_peak_size Peaks below this size would be discarded.
#'
#' @examples
#' \dontrun{
#' split_peaks <- split_long_peaks("pbmc_atac.marginal", peaks = peaks, target_size = 500, max_peak_size = 1e3, very_long = 5e3, min_peak_size = 20)
#' }
#'
#' @export
split_long_peaks <- function(marginal_track, peaks, target_size = 500, max_peak_size = 1e3, very_long = 5e3, min_peak_size = NULL) {
    if (!gtrack.exists(marginal_track)) {
        cli_abort("Track {.val {marginal_track}} does not exist")
    }

    if (target_size > max_peak_size) {
        cli_abort("{.field target_size} should be less or equal to {.field max_peak_size}")
    }

    if (very_long < max_peak_size) {
        cli_abort("{.field very_long} should be higher or equal to {.field max_peak_size}")
    }

    cli_li("Splitting peaks which are longer than {.val {very_long}} into smaller peaks of size {.val {max_peak_size}}")
    peaks <- split_peaks_arbitrarily(peaks, very_long)

    peaks_to_split <- peaks %>%
        filter(end - start > max_peak_size)
    other_peaks <- peaks %>%
        filter(end - start <= max_peak_size)

    binsize <- gtrack.info(marginal_track)$bin.size
    cli_li("Extracting {.val {marginal_track}}")
    withr::with_options(list("gmax.data.size" = 1e9), {
        peaks_to_split$intervalID <- 1:nrow(peaks_to_split)
        mar_df <- gextract(marginal_track, intervals = peaks_to_split, colnames = "cov", iterator = binsize) %>%
            as_tibble()
        mar_df <- mar_df %>%
            group_by(intervalID) %>%
            mutate(i = 1:n(), max_i = max(i)) %>%
            ungroup()
    })
    max_len <- max(peaks$end - peaks$start)
    max_bins <- ceiling(max_len / binsize)

    cli_li("Finding best offset for each peak")
    tri_mat <- get_peak_triangles(target_size, binsize, max_bins)

    cov_mat <- mar_df %>%
        select(intervalID, i, cov) %>%
        tidyr::spread(i, cov) %>%
        column_to_rownames("intervalID") %>%
        as.matrix()

    cov_mat <- cov_mat - rowMeans(cov_mat, na.rm = TRUE)
    cov_mat[cov_mat < 0] <- NA

    tri_cors <- tgs_cor(t(cov_mat), tri_mat, spearman = TRUE, pairwise.complete.obs = TRUE)

    row_best_match <- apply(tri_cors, 1, which.max)
    best_match_tri <- t(tri_mat[, row_best_match])
    rownames(best_match_tri) <- rownames(cov_mat)
    colnames(best_match_tri) <- colnames(cov_mat)

    offset_df <- best_match_tri %>%
        as.data.frame() %>%
        rownames_to_column("intervalID") %>%
        tidyr::gather(i, ind, -intervalID) %>%
        mutate(intervalID = as.numeric(intervalID), i = as.numeric(i)) %>%
        as_tibble()

    cli_li("Splitting peaks into chunks of {.val {target_size}}")
    mar_df <- mar_df %>% left_join(offset_df, by = c("intervalID", "i"))
    mar_df <- mar_df %>%
        left_join(peaks_to_split %>% mutate(intervalID = 1:n()) %>% select(chrom1 = chrom, start1 = start, end1 = end, intervalID), by = "intervalID")

    # we have a cut point at the start of the peak, at the end of the peak, and at each start of a triangle
    mar_df <- mar_df %>%
        mutate(cut_point = ifelse(ind == 1 | i == 1 | i == max_i, 0, 1))

    mar_df <- mar_df %>%
        # when we have two successive cut points, take only the second one
        mutate(cut_point = ifelse(intervalID == lead(intervalID) & ind == lead(ind) & cut_point == 0, 1, cut_point)) %>%
        # The last point became NA, so we have to define it as a cut point
        tidyr::replace_na(replace = list(cut_point = 0)) %>%
        # instead of using an entire bin as a cut point - we use a single bp
        mutate(end = ifelse(cut_point == 0 & i != max_i, start + 1, end)) %>%
        mutate(start = ifelse(cut_point == 0 & i == max_i, end - 1, start))

    # make sure that there are no gaps
    mar_df <- mar_df %>%
        group_by(intervalID) %>%
        mutate(start = pmin(lag(end) + 1, start, na.rm = TRUE)) %>%
        ungroup()

    # use gscreen to extract the new peaks from the cutting points
    temp_track <- paste0("temp.", basename(tempfile()))
    gdir.create("temp", showWarnings = FALSE)
    gtrack.create_sparse(temp_track, intervals = mar_df %>% select(chrom, start, end), values = mar_df$cut_point, description = "temp")
    withr::defer(gtrack.rm(temp_track, force = TRUE))
    splitted_peaks <- gscreen(glue("{temp_track} > 0"), intervals = peaks_to_split)

    # merge the splitted peaks with the original peaks
    splitted_peaks <- splitted_peaks %>%
        gintervals.neighbors1(peaks_to_split) %>%
        # again make sure that there are no gaps
        group_by(intervalID) %>%
        mutate(start = pmin(lag(end + 1), start, na.rm = TRUE)) %>%
        mutate(i = 1:n(), max_i = max(i)) %>%
        ungroup() %>%
        # make sure that there are no gaps at the beggining and end of the original peak
        mutate(start = ifelse(i == 1, start1, start), end = ifelse(i == max_i, end1, end))


    new_peaks <- bind_rows(
        other_peaks,
        splitted_peaks %>% select(chrom, start, end)
    ) %>% arrange(chrom, start, end)

    if (!is.null(min_peak_size)) {
        new_peaks <- new_peaks %>%
            mutate(l = end - start) %>%
            filter(l >= min_peak_size)
    }

    return(new_peaks)
}

get_peak_triangles <- function(target_size, binsize, max_bins) {
    tri_binsize <- floor(target_size / binsize)
    if (tri_binsize %% 2 == 0) {
        peak_side_n <- (tri_binsize - 2) / 2
        triangle <- c(1:peak_side_n, peak_side_n + 1, peak_side_n + 1, peak_side_n:1)
    } else {
        peak_side_n <- (tri_binsize - 1) / 2
        triangle <- c(1:peak_side_n, peak_side_n + 1, peak_side_n:1)
    }

    tri_vec <- rep(triangle, max_bins * 100)
    tri_mat <- sapply(1:length(triangle), function(x) tri_vec[x:(x + max_bins - 1)])
    return(tri_mat)
}
