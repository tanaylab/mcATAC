#' McCounts object
#'
#' @description a McCounts object is a collection of sparse matrices where rows are genomic coordinates and columns are metacells. It is derived from ScCounts using \code{scc_project_on_mc}. The metacell names are store in the \code{cell_names} slot, and the only additional slot is \code{cell_to_metacell} which contains the mapping from single cell names to metacell names.
#'
#' @slot cell_to_metacell
#'
#' @include ScCounts.R
#' @rdname ScCounts
#' @exportClass McCounts
McCounts <- setClass(
    "McCounts",
    slots = c(
        cell_to_metacell = "data.frame_or_null"
    ),
    contains = "ScCounts"
)

setMethod(
    "initialize",
    signature = "McCounts",
    definition = function(.Object, cell_to_metacell, ...) {
        .Object@cell_to_metacell <- cell_to_metacell
        return(callNextMethod(.Object, ...))
    }
)


#' @export
#' @noRd
setMethod(
    "show",
    signature = "McCounts",
    definition = function(object) {
        print_counts_object(object, "McCounts", "metacell")
    }
)

#' @export
#' @noRd
setMethod("Ops", c("ScCounts", "ScCounts"), function(e1, e2) {
    if (e1@genome != e2@genome) {
        cli_abort("Genomes are not the same.")
    }
    if (any(e1@genomic_bins$name != e2@genomic_bins$name)) {
        cli_abort("Genomic bins are not the same.")
    }

    both <- intersect(e1@cell_names, e2@cell_names)
    only_e1 <- setdiff(e1@cell_names, e2@cell_names)
    only_e2 <- setdiff(e2@cell_names, e1@cell_names)
    all_cells <- unique(c(both, only_e1, only_e2))

    matrices <- plyr::llply(e1@genomic_bins$name, function(bin) {
        m1 <- e1@data[[bin]]
        m2 <- e2@data[[bin]]

        intersect_m <- m1[, both, drop = FALSE] + m2[, both, drop = FALSE]
        m <- cbind(intersect_m, m1[, only_e1, drop = FALSE], m2[, only_e2, drop = FALSE])
        m <- m[, all_cells, drop = FALSE]
        m
    }, .parallel = getOption("mcatac.parallel"))

    names(matrices) <- e1@genomic_bins$name

    res <- new("ScCounts",
        data = matrices,
        genome = e1@genome,
        genomic_bins = e1@genomic_bins,
        cell_names = all_cells
    )

    return(res)
})

#' @export
#' @noRd
setMethod("Ops", c("McCounts", "McCounts"), function(e1, e2) {
    scc <- callNextMethod()
    cell_to_metacell <- rbind(e1@cell_to_metacell, e2@cell_to_metacell) %>%
        distinct(cell_id, metacell)
    res <- new("McCounts",
        data = scc@data,
        genome = scc@genome,
        genomic_bins = scc@genomic_bins,
        cell_names = scc@cell_names,
        cell_to_metacell = cell_to_metacell
    )
    return(res)
})


#' Given metacells (usually from RNA data), project ATAC counts from multiple batches to get a McCounts object
#'
#' @param scc_dirs a vector of paths to directories containing ScCounts objects
#'
#' @inheritParams scc_to_mcc
#'
#' @export
scc_to_mcc_multi_batch <- function(scc_dirs, cell_to_metacell, ignore_metacells = c(-1, -2)) {
    cell_to_metacell <- read_cell_to_metacell(cell_to_metacell, ignore_metacells)

    add_scc <- function(dir, x = NULL) {
        if (!dir.exists(dir)) {
            cli_abort("Directory {.x} does not exist")
        }
        scc <- scc_read(dir)
        cell_to_metacell_batch <- cell_to_metacell %>%
            filter(cell_id %in% scc@cell_names)
        mcc1 <- scc_to_mcc(scc, cell_to_metacell_batch, ignore_metacells)
        if (!is.null(x)) {
            return(x + mcc1)
        }
        return(mcc1)
    }

    mcc <- NULL
    for (dir in scc_dirs) {
        mcc <- add_scc(dir, mcc)
    }
    return(mcc)
}

read_cell_to_metacell <- function(cell_to_metacell, ignore_metacells) {
    if (is.character(cell_to_metacell)) {
        if (!file.exists(cell_to_metacell)) {
            cli_abort("File {cell_to_metacell} does not exist.")
        }
        cell_to_metacell <- anndata::read_h5ad(cell_to_metacell)$obs

        if (!"metacell" %in% colnames(cell_to_metacell)) {
            cli_abort("Column 'metacell' not found in {.field cell_to_metacell}.")
        }

        cell_to_metacell <- cell_to_metacell %>%
            tibble::rownames_to_column("cell_id") %>%
            select(cell_id, metacell)
    }

    if (any(cell_to_metacell$metacell %in% ignore_metacells)) {
        ignored <- cell_to_metacell$metacell[cell_to_metacell$metacell %in% ignore_metacells]
        cli_alert_warning("Ignoring metacells: {.val {ignore_metacells}}")
        cell_to_metacell <- cell_to_metacell %>% filter(!(metacell %in% ignored))
    }

    return(cell_to_metacell)
}

#' Given metacells (usually from RNA data), project ATAC counts to get a McCounts object
#'
#' @description Given cell to metacell association, summarise atac read counts to generate a McCounts object. This can
#' take a while - around 5 minutes using 24 cores on the PBMC dataset.
#'
#' @param sc_counts A ScCounts object
#' @param cell_to_metacell a data.frame with columns \code{cell_id} and \code{metacell} containing the mapping from single cell names to metacell names, or the name of an 'h5ad' file containing this information at the 'obs' slot. In such a case, the 'obs' slot should contain
#' a column named \code{metacell} and the rownames should be the cell names.
#' @param ignore_metacells a vector of metacells to ignore. Default: [-1, -2] (the "outliers" metacell in the metacell2 python package).
#'
#' @return A McCounts object
#'
#' @examples
#' \dontrun{
#' data(cell_to_metacell_pbmc_example)
#' scc_to_mcc(sc_counts, cell_to_metacell_pbmc_example)
#' }
#'
#' @export
scc_to_mcc <- function(sc_counts, cell_to_metacell, ignore_metacells = c(-1, -2)) {
    assert_atac_object(sc_counts, class = "ScCounts")

    cell_to_metacell <- read_cell_to_metacell(cell_to_metacell, ignore_metacells)
    cell_to_metacell <- deframe(cell_to_metacell)
    assert_that(all(names(cell_to_metacell) %in% sc_counts@cell_names))

    removed_cells <- setdiff(sc_counts@cell_names, names(cell_to_metacell))
    if (length(removed_cells) > 0) {
        cli_alert_info("{.val {length(removed_cells)}} cells (out of {.val {length(sc_counts@cell_names)}}) do not have a metacell and have been removed.")
    }

    cell_to_metacell <- cell_to_metacell[intersect(names(cell_to_metacell), sc_counts@cell_names)]

    new_data <- plyr::llply(sc_counts@data, function(sc_mat) {
        cells <- intersect(colnames(sc_mat), names(cell_to_metacell))
        mc_mat <- sparse_matrix_tapply_sum(sc_mat[, cells], cell_to_metacell[cells])
        return(mc_mat)
    }, .parallel = getOption("mcatac.parallel"))

    res <- new("McCounts", data = new_data, cell_names = as.character(sort(unique(cell_to_metacell))), genome = sc_counts@genome, genomic_bins = sc_counts@genomic_bins, id = sc_counts@id, description = sc_counts@description, path = sc_counts@path, cell_to_metacell = enframe(cell_to_metacell, "cell_id", "metacell"))

    return(res)
}


#' @rdname scc_to_mcc
#' @export
scc_project_on_mc <- function(sc_counts, cell_to_metacell, ignore_metacells = c(-1, -2)) {
    lifecycle::deprecate_soft(
        when = "0.0.1",
        what = "scc_project_on_mc()",
        with = "scc_to_mcc()",
    )

    scc_to_mcc(sc_counts, cell_to_metacell, ignore_metacells)
}

#' Read an McCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{scc_from_bam})
#'
#' @return A McCounts object
#'
#' @examples
#' \dontrun{
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' }
#'
#' @inheritParams scc_read
#' @export
mcc_read <- function(path, id = NULL, description = NULL, verbose = TRUE) {
    md_file <- file.path(path, "metadata.yaml")
    if (!file.exists(md_file)) {
        cli_abort("Directory {.file {path}} does not contain a valid McCounts object")
    }

    md <- yaml::read_yaml(md_file)

    if (is.null(md$cell_to_metacell)) {
        cli_abort("Directory {.file {path}} does not contain a valid McCounts object (the metadata file {.file {md_file}} is missing the {.field cell_to_metacell} field.)")
    }

    sc_counts <- read_counts_object(path, "McCounts", id, description, verbose = FALSE)

    mc_counts <- new("McCounts", data = sc_counts@data, cell_names = sc_counts@cell_names, genome = sc_counts@genome, genomic_bins = sc_counts@genomic_bins, id = sc_counts@id, description = sc_counts@description, path = sc_counts@path, cell_to_metacell = as_tibble(md$cell_to_metacell))

    if (verbose) {
        cli_alert_success("Succesfully read a McCounts object from {.file {path}}")
    }

    return(mc_counts)
}

#' Given a sparse matrix of genomic counts and intervals, return a sparse matrix of peaks
#'
#' @param mat A sparse matrix of genomic counts (rows are 1bp, columns are metacells)
#' @param bin a single row data frame with the bin intervals (chrom, start, end)
#' @param intervs an intervals set with another column called "peak_name"
#' @param metacells names of metacells to include. Default: all metacells.
#'
#' @return a sparse matrix where rows are the peaks, columns are the metacells, and values are the counts summed
#' over each peak.
#'
#' @noRd
summarise_bin <- function(mat, bin, intervs, metacells = NULL) {
    metacells <- metacells %||% colnames(mat)
    intervs <- as.data.frame(intervs)
    if (!has_name(intervs, "peak_name")) {
        cli_abort("The {.var intervs} must have a column called {.field peak_name}")
    }
    # find the coordinates that overlap with the intervals (note that we assume mat@i is zero based, as in the standard 'dgCMatrix' implementation)
    mat_intervs <- tibble(chrom = bin$chrom, start = unique(mat@i) + bin$start, end = start + 1) %>%
        gintervals.intersect(intervs)

    # no overlap, return empty matrix
    if (is.null(mat_intervs)) {
        msum <- Matrix::Matrix(
            nrow = length(intervs$peak_name), ncol = length(metacells), data = 0,
            sparse = TRUE
        )
        msum <- as(msum, "dgCMatrix")
        rownames(msum) <- intervs$peak_name
        colnames(msum) <- metacells
        return(msum)
    }

    mat_intervs <- giterator.intervals(intervals = mat_intervs, iterator = 1)

    # make sure we do not take coordinates outside the bin
    mat_intervs <- mat_intervs %>%
        mutate(start = pmax(start, bin$start), end = pmin(end, bin$end))

    mat_intervs <- mat_intervs %>%
        misha.ext::gintervals.neighbors1(intervs) %>%
        # although mat@i is 0-based, R indices are 1-based (hence the +1)
        mutate(ind = start - bin$start + 1)

    group <- factor(mat_intervs$peak_name, levels = intervs$peak_name)
    mat_f <- mat[mat_intervs$ind, metacells, drop=FALSE]
    res <- t(sparse_matrix_tapply_sum(t(mat_f), group))

    return(res)
}

#' Create an McPeaks object from an McCounts object
#'
#' @description given an McCounts object and peaks, summarise the counts over the peaks and return a McPeaks object
#'
#' @param mc_counts a McCounts object
#' @param peaks a data frame with the peak intervals (chrom, start, end) and a column called "peak_name"
#' @param metacells names of metacells to include. Default: all metacells.
#'
#' @inheritParams project_atac_on_mc
#'
#' @return a McPeaks object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' atac_mc <- mcc_to_mcatac(mc_counts, atac_sc@peaks)
#' }
#'
#' @export
mcc_to_mcatac <- function(mc_counts, peaks, metacells = NULL, metadata = NULL, mc_size_eps_q = 0.1) {
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)
    if (!has_name(peaks, "peak_name")) {
        peaks <- peaks %>% mutate(peak_name = peak_names(.))
        cli_alert_warning("The {.var peaks} didn't have a column called {.field peak_name}, so it was added using the {.code peak_names} function.")
    }
    matrices <- plyr::alply(mc_counts@genomic_bins, 1, function(bin) {
        return(
            summarise_bin(mc_counts@data[[bin$name]], bin, peaks, metacells)
        )
    }, .parallel = getOption("mcatac.parallel"))

    mat <- Reduce("+", matrices)

    mc_atac <- new("McPeaks", mat = mat, peaks = peaks, genome = mc_counts@genome, id = mc_counts@id, description = mc_counts@description, metadata = metadata, cell_to_metacell = mc_counts@cell_to_metacell, mc_size_eps_q = mc_size_eps_q, path = mc_counts@path)

    cli_alert_success("Created a new McPeaks object with {.val {ncol(mc_atac@mat)}} metacells and {.val {nrow(mc_atac@mat)}} ATAC peaks.")

    return(mc_atac)
}

#' Extract summary statistics per cell for a set of intervals
#'
#' @param scc a ScCounts object
#' @param intervals an intervals set. Can have a column called "peak_name" with the peak name.
#' @param cells a vector of cell names to include. Default: all cells.
#'
#' @return a sparse matrix where rows are the intervals, columns are the cells, and values are the counts summed. If the intervals have a column called "peak_name", the rows will be the peak names, otherwise the rows will be of the form "{chr}:{start}_{end}"
#'
#' @examples
#' \dontrun{
#' scc_extract(scc, gintervals(1, 0, 100))
#' }
#'
#' @export
scc_extract <- function(scc, intervals, cells = NULL) {
    assert_atac_object(scc, class = "ScCounts")
    cells <- cells %||% scc@cell_names
    cells <- as.character(cells)

    if (!has_name(intervals, "peak_name")) {
        intervals <- intervals %>% mutate(peak_name = peak_names(., tad_based = FALSE))
    }

    matrices <- plyr::alply(scc@genomic_bins, 1, function(bin) {
        return(
            summarise_bin(scc@data[[bin$name]], bin, intervals, cells)
        )
    }, .parallel = getOption("mcatac.parallel"))

    mat <- Reduce("+", matrices)

    return(mat)
}

#' Return the total coverage of each metacell in an McCounts object
#'
#' @param mc_counts a McCounts object
#' @param metacells names of metacells to include. Default: all metacells.
#'
#' @return a named vector with the total coverage for each metacell
#'
#' @examples
#' \dontrun{
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' mc_covs <- mcc_cov(mc_counts)
#' }
#'
#' @export
mcc_metacell_total_cov <- function(mc_counts, metacells = NULL) {
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    sums <- plyr::llply(mc_counts@data, function(m) {
        colSums(m[, metacells])
    }, .parallel = getOption("mcatac.parallel"))

    mc_covs <- Reduce("+", sums)
    return(mc_covs)
}

#' Return the total coverage of each non-zero coordinate in an ScCounts/McCounts object
#'
#' @param mc_counts a McCounts object
#' @param metacells names of metacells to include. Default: all metacells.
#' @param cells names of cells to include. Default: all cells.
#'
#' @return an intervals set with an additional column called "cov" with the total coverage for each coordinate
#'
#' @examples
#' \dontrun{
#' # cell counts
#' sc_counts <- scc_read("pbmc_reads")
#' sc_mars <- scc_marginal(sc_counts)
#'
#' # metacell counts
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' mc_mars <- mcc_marginal(mc_counts)
#' }
#'
#' @export
mcc_marginal <- function(mc_counts, metacells = NULL) {
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    mars <- plyr::adply(mc_counts@genomic_bins, 1, function(bin) {
        mat <- mc_counts@data[[bin$name]]
        mat <- mat[, metacells]
        res <- tibble(chrom = bin$chrom, start = unique(mat@i) + bin$start, end = start + 1) %>%
            mutate(cov = rowSums(mat[end - bin$start, ]))
    }, .parallel = getOption("mcatac.parallel"))

    mars <- mars %>%
        select(chrom, start, end, cov) %>%
        as_tibble()

    return(mars)
}

#' @rdname mcc_marginal
scc_marginal <- function(sc_counts, cells = NULL) {
    mcc_marginal(sc_counts, cells)
}

#' Normalize each metacell in a McCounts object by its total counts
#'
#' @param mc_counts a McCounts object
#' @param metacells metacells to normalize. Default: all metacells.
#'
#' @return an McCounts object with the metacells normalized
#'
#' @examples
#' \dontrun{
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' mc_counts_norm <- mcc_normalize(mc_counts)
#' }
#'
#' @export
mcc_normalize_metacells <- function(mc_counts, metacells = NULL) {
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    mc_covs <- mcc_metacell_total_cov(mc_counts, metacells)

    new_data <- plyr::llply(mc_counts@data, function(m) {
        m_f <- m[, metacells]
        m_f@x <- m_f@x / rep.int(mc_covs, diff(m_f@p))
        return(m_f)
    }, .parallel = getOption("mcatac.parallel"))

    stopifnot(all(names(new_data) == names(mc_counts@data)))

    res <- mc_counts
    res@data <- new_data
    return(res)
}

#' Extract an intervals set from a counts matrix
#'
#' @param mat a counts sparse matrix (e.g. from McCounts object)
#' @param bin a data frame with a single line with the genomic bin coordinates
#' @param metacells names of metacells to include.
#' @param all_metacells names of all metacells in the counts matrix.
#'
#' @return a data frame with "chrom", "start", "metacell" and "value" columns
#'
#' @noRd
extract_bin <- function(mat, bin, metacells, all_metacells) {
    Matrix::summary(mat[, metacells, drop = FALSE]) %>%
        mutate(
            start = i + bin$start - 1, # we remove 1 in order to transform to 0 based coordinates
            end = start + 1,
            chrom = bin$chrom,
            metacell = all_metacells[j]
        ) %>%
        select(chrom, start, end, metacell, value = x) %>%
        as_tibble()
}

#' Extract McCounts to a tidy data frame
#'
#' @param mc_counts A McCounts object
#' @param metacells metacells to extract.
#'
#' @return an intervals set (tibble) with additional "metacell" and "value" columns
#'
#' @examples
#' \dontrun{
#' df <- mcc_extract_to_df(mc_counts, c("1", "2"))
#' }
#'
#' @export
mcc_extract_to_df <- function(mc_counts, metacells = NULL) {
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    mc_data <- plyr::adply(mc_counts@genomic_bins, 1, function(bin) {
        return(
            extract_bin(mc_counts@data[[bin$name]], bin, metacells, mc_counts@cell_names)
        )
    }, .parallel = getOption("mcatac.parallel"))

    mc_data <- mc_data %>%
        select(chrom, start, end, metacell, value) %>%
        as_tibble()

    return(mc_data)
}

#' Create a misha track for each metacell in a McCounts object
#'
#' @description This would create a sparse track of the form "{track_prefix}.mc{metacell}" for each metacell in the McCounts object. The counts are smoothed by a running sum window of size \code{window_size * 2 + 1} (default: 201).
#' The resulted track would be stored in a 'dense' format with a resolution of \code{resolution} bp (default: 10). \cr
#' Each track would have a track attribute named "total_cov" with the total coverage of the metacell. \cr
#' If \code{create_marginal_track} is \code{TRUE} (default), an additional track with the total counts over all the metacells will be created. \cr
#' Counts can be normalized to the total counts of each metacell by setting \code{normalize = TRUE}.
#'
#' @param mc_counts A McCounts object
#' @param track_prefix The prefix of the tracks to create (i.e. misha directory).
#' Track names will be of the form "{track_prefix}.mc{metacell}"
#' @param metacells metacells for which to create tracks. If NULL, all metacells will be used.
#' @param overwrite Whether to overwrite existing tracks.
#' @param resolution The resolution of each track (default: 20bp)
#' @param window_size The size of the window used to smooth the counts. The counts of at position i are smoothed by a sum of the counts of [i - window_size, i + window_size]. If NULL - the counts are not smoothed.
#' @param create_marginal_track Create a track with the total counts from all the metacells. The track would be named "{track_prefix}.marginal"
#' @param normalize Normalize each metacell by the sum of its counts.
#'
#'
#' @return An McTracks object with the new tracks.
#'
#' @examples
#' \dontrun{
#' mcc_to_tracks(mc_counts, "pbmc_mc")
#' }
#'
#' @export
mcc_to_tracks <- function(mc_counts, track_prefix, metacells = NULL, overwrite = FALSE, resolution = 20, window_size = NULL, create_marginal_track = TRUE, normalize = FALSE) {
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)
    gset_genome(mc_counts@genome)

    cli_alert_info("Creating tracks for {.val {length(metacells)}} metacells")
    cli_alert_info("Smoothing over {.val {window_size*2+1}} bp window")
    cli_alert_info("Tracks resolution: {.val {resolution}} bp")

    misha.ext::gtrack.create_dirs(paste0(track_prefix, ".mc"), showWarnings = FALSE)
    withr::local_options(list(gmax.data.size = 1e9))

    marginal_track <- NULL
    if (create_marginal_track) {
        marginal_track <- paste0(track_prefix, ".marginal")
        cli_alert("Creating marginal track")
        mcc_to_marginal_track(
            mc_counts,
            track = marginal_track,
            metacells = metacells,
            window_size = window_size,
            resolution = resolution,
            overwrite = overwrite
        )
    }

    if (normalize) {
        cli_alert("Normalizing each metacell by its total counts")
        mc_counts <- mcc_normalize_metacells(mc_counts, metacells)
    }

    withr::local_options(list(gmultitasking = !getOption("mcatac.parallel")))

    plyr::l_ply(mc_counts@cell_names, function(metacell) {
        track <- glue("{track_prefix}.mc{metacell}")
        cli_alert("Creating {track} track")

        withr::with_options(list(mcatac.parallel = FALSE), {
            x <- mcc_extract_to_df(mc_counts, metacell)
        })

        if (is.null(window_size)) {
            description <- glue("Counts from metacell {metacell} in {resolution}bp resolution")
        } else {
            description <- glue("Smoothed counts over {window_size*2} bp from metacell {metacell}")
        }

        create_smoothed_track_from_dataframe(
            x %>% select(chrom, start, end, value),
            track_prefix = track_prefix,
            track = track,
            description = description,
            window_size = window_size,
            resolution = resolution,
            overwrite = overwrite
        )

        gtrack.attr.set(track, "total_cov", sum(x$value))
        num_cells <- mc_counts@cell_to_metacell %>%
            filter(metacell == !!metacell) %>%
            nrow()
        gtrack.attr.set(track, "num_cells", num_cells)
        gc()
    }, .parallel = getOption("mcatac.parallel"))

    gdb.reload()

    cli_alert_success("Created {length(metacells)} tracks at {track_prefix}")

    mct <- mct_create(genome = mc_counts@genome, tracks = glue("{track_prefix}.mc{metacells}"), metacells = metacells, id = mc_counts@id, description = mc_counts@description, path = mc_counts@path, metadata = mc_counts@metadata, resolution = resolution, window_size = window_size, marginal_track = marginal_track)
    return(mct)
}

#' Create a track with smoothed marginal counts from an ScCounts/McCounts object
#'
#' @param mc_counts An McCounts object
#' @param sc_counts An ScCounts object
#' @param track The name of the track to create.
#' @param metacells The metacells for which to create the track. If NULL, all metacells will be used.
#' @param cells The cells for which to create the track.
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' mcc_to_marginal_track(mc_counts, "pbmc_mc.marginal")
#' }
#'
#' @inheritParams mcc_to_tracks
#'
#' @export
mcc_to_marginal_track <- function(mc_counts, track, metacells = NULL, resolution = 10, window_size = 100, overwrite = FALSE) {
    assert_atac_object(mc_counts, class = "ATAC")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)
    gset_genome(mc_counts@genome)

    if (is.null(window_size)) {
        description <- glue("Marginal counts from {length(metacells)} metacells")
    } else {
        description <- glue("Smoothed marginal counts over {window_size*2} bp from {length(metacells)} metacells")
    }

    mcc_mars <- mcc_marginal(mc_counts, metacells)
    gdir.create("temp", showWarnings = FALSE)
    create_smoothed_track_from_dataframe(
        mcc_mars,
        track_prefix = "temp",
        track = track,
        description = description,
        window_size = window_size,
        resolution = resolution,
        overwrite = overwrite
    )
}

#' @rdname mcc_to_marginal_track
#' @export
scc_to_marginal_track <- function(sc_counts, track, cells = NULL, resolution = 10, window_size = 100, overwrite = FALSE) {
    mcc_to_marginal_track(sc_counts, track, cells, resolution, window_size, overwrite)
}

#' Create a smoothed track from a data frame
#'
#' @description This function first creates a temporary track with the data in the data frame, then calls \code{create_smoothed_track} to create the smoothed track.
#'
#' @param x an intervals set with an additional value column
#' @param track_prefix The prefix of the track to create (i.e. misha directory).
#'
#' @return None.
#'
#' @inheritDotParams create_smoothed_track
#' @noRd
create_smoothed_track_from_dataframe <- function(x, track_prefix, ...) {
    raw_track <- temp_track_name(paste0(track_prefix, "."))
    x <- as_tibble(x)
    gtrack.create_sparse(
        track = raw_track,
        description = "",
        intervals = x[, 1:3],
        values = x[[4]]
    )

    create_smoothed_track(raw_track = raw_track, ...)
}

#' Create a smoothed track from a sparse track
#'
#' @description Create a track smoothed by a sum of [i - window_size, i + window_size], where i is the middle of the \code{resolution}. For
#' examples, if \code{resolution} is 5, and \code{window_size} is 10, when the iterator is at position 20 the smoothed track will be created by summing the values of positions 10-29.
#'
#' @param raw_track name of the track to smooth.
#' @param track The name of the track to create.
#' @param description The description of the track.
#' @param resolution The resolution of the track (in bp)
#' @param window_size The size of the window to smooth (in bp). If NULL, no smoothing will be performed.
#' @param overwrite Overwrite the track if it already exists.
#'
#' @return None
#'
#' @noRd
create_smoothed_track <- function(raw_track, track, description, resolution, window_size = NULL, overwrite = FALSE) {
    if (gtrack.exists(track)) {
        if (overwrite) {
            cli_alert_warning("Removing previous track {.val {track}}")
            gtrack.rm(track, force = TRUE)
        } else {
            cli_abort("{track} already exists. Use 'overwrite = TRUE' to overwrite.")
        }
    }
    vt <- glue("vt_{raw_track}")
    gvtrack.create(vt, raw_track, func = "sum")
    withr::defer(gvtrack.rm(vt))

    if (!is.null(window_size)) {
        shift <- window_size - floor(resolution / 2)
        gvtrack.iterator(vt, sshift = -shift, eshift = shift)
    }

    gtrack.create(
        track,
        description = description,
        expr = vt,
        iterator = resolution
    )
}
