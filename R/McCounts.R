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

#' Given metacells (usually from RNA data), project ATAC counts to get a McCounts object
#'
#' @description Given cell to metacell association, summarise atac read counts to generate a McCounts object. This can
#' take a while - around 5 minutes using 24 cores on the PBMC dataset.
#'
#' @param sc_counts A ScCounts object
#' @param cell_to_metacell a data frame with a column named "cell_id" with cell id and another column named "metacell" with the metacell the cell is part of.
#' @param ignore_metacells a vector of metacells to ignore. Default: [-1] (the "outliers" metacell in the metacell2 python package).
#'
#' @return A McCounts object
#'
#' @examples
#' \dontrun{
#' data(cell_to_metacell_pbmc_example)
#' scc_project_on_mc(sc_counts, cell_to_metacell_pbmc_example)
#' }
#'
#' @export
scc_project_on_mc <- function(sc_counts, cell_to_metacell, ignore_metacells = -1) {
    assert_atac_object(sc_counts, class = "ScCounts")
    if (any(cell_to_metacell$metacell %in% ignore_metacells)) {
        ignored <- cell_to_metacell$metacell[cell_to_metacell$metacell %in% ignore_metacells]
        cli_alert_warning("Ignoring metacells: {.val {ignore_metacells}}")
        cell_to_metacell <- cell_to_metacell %>% filter(!(metacell %in% ignored))
    }
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

#' Read an McCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{write_sc_counts_from_bam})
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
    sc_counts <- scc_read(path, id, description, verbose = FALSE)
    md_file <- file.path(path, "metadata.yaml")
    md <- yaml::read_yaml(md_file)

    if (is.null(md$cell_to_metacell)) {
        cli_abort("Directory {.file {path}} does not contain a valid McCounts object (the metadata file {.file {md_file}} is missing the {.field cell_to_metacell} field.)")
    }

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

    # find the coordinates that overlap with the intervals (note that we mat@i is zero based)
    mat_intervs <- tibble(chrom = bin$chrom, start = mat@i + bin$start, end = start + 1) %>%
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
    mat_f <- mat[mat_intervs$ind, metacells]
    res <- t(sparse_matrix_tapply_sum(t(mat_f), group))

    return(res)
}

#' Create an McATAC object from an McCounts object
#'
#' @description given an McCounts object and peaks, summarise the counts over the peaks and return a McATAC object
#'
#' @param mc_counts a McCounts object
#' @param peaks a data frame with the peak intervals (chrom, start, end) and a column called "peak_name"
#' @param metacells names of metacells to include. Default: all metacells.
#'
#' @inheritParams project_atac_on_mc
#'
#' @return a McATAC object
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

    matrices <- plyr::alply(mc_counts@genomic_bins, 1, function(bin) {
        return(
            summarise_bin(mc_counts@data[[bin$name]], bin, peaks, metacells)
        )
    }, .parallel = getOption("mcatac.parallel"))

    mat <- Reduce("+", matrices)

    mc_atac <- new("McATAC", mat = mat, peaks = peaks, genome = mc_counts@genome, id = mc_counts@id, description = mc_counts@description, metadata = metadata, cell_to_metacell = mc_counts@cell_to_metacell, mc_size_eps_q = mc_size_eps_q, path = mc_counts@path)

    cli_alert_success("Created a new McATAC object with {.val {ncol(mc_atac@mat)}} metacells and {.val {nrow(mc_atac@mat)}} ATAC peaks.")

    return(mc_atac)
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
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    sums <- plyr::llply(mc_counts@data, function(m) {
        colSums(m[, metacells])
    }, .parallel = getOption("mcatac.parallel"))

    mc_covs <- Reduce("+", sums)
    return(mc_covs)
}


#' Normalize each metacell in a McCounts object by its total counts
#'
#' @param mc_counts a McCounts object
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
mcc_normalize_metacells <- function(mc_counts) {
    assert_atac_object(mc_counts, class = "McCounts")

    mc_covs <- mcc_metacell_total_cov(mc_counts)

    new_data <- plyr::llply(mc_counts@data, function(m) {
        t(t(m) / mc_covs)
    }, .parallel = TRUE)

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
#' @description This would create a sparse track of the form "{track_prefix}.mc{metacell}" for each metacell in the McCounts object. By default, the counts are normalized to the total counts of each metacell, set \code{normalize = FALSE}
#' to disable this.
#'
#' @param mc_counts A McCounts object
#' @param track_prefix The prefix of the tracks to create. Track names will be of the form "{track_prefix}.mc{metacell}"
#' @param metacells metacells for which to create tracks. If NULL, all metacells will be used.
#' @param normalize Normalize each metacell by the sum of its counts.
#' @param overwrite Whether to overwrite existing tracks.
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' mcc_to_tracks(mc_counts, "pbmc_mc")
#' }
#'
#' @export
mcc_to_tracks <- function(mc_counts, track_prefix, metacells = NULL, overwrite = FALSE, normalize = TRUE) {
    assert_atac_object(mc_counts, class = "McCounts")
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)
    gset_genome(mc_counts@genome)

    if (normalize) {
        cli_alert("Normalizing each metacell by its total counts")
        mc_counts <- mcc_normalize_metacells(mc_counts)
    }

    cli_alert("Extracting metacell data")
    d <- mcc_extract_to_df(mc_counts, metacells)

    misha.ext::gtrack.create_dirs(paste0(track_prefix, ".mc"), showWarnings = FALSE)

    plyr::l_ply(mc_counts@cell_names, function(metacell) {
        track <- glue("{track_prefix}.mc{metacell}")
        cli_alert("Creating {track} track")
        x <- d %>% filter(metacell == !!metacell)
        if (gtrack.exists(track)) {
            if (overwrite) {
                cli_alert_warning("Removing previous track {.val {track}}")
                gtrack.rm(track, force = TRUE)
            } else {
                cli_abort("{track} already exists. Use 'overwrite = TRUE' to overwrite.")
            }
        }
        gtrack.create_sparse(
            track = track,
            description = "",
            intervals = x %>% select(chrom, start, end),
            values = x$value
        )
    }, .parallel = getOption("mcatac.parallel"))

    gdb.reload()

    cli_alert_success("Created {length(metacells)} tracks at {track_prefix}")
}
