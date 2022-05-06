#' McCounts object
#'
#' @description a McCounts object is a collection of sparse matrices where rows are genomic coordinates and columns are metacells. It is derived from ScCounts using \code{project_counts_on_mc}. The metacell names are store in the \code{cell_names} slot, and the only additional slot is \code{cell_to_metacell} which contains the mapping from single cell names to metacell names.
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
#' @param num_cores The number of cores to use.
#'
#' @return A McCounts object
#'
#' @examples
#' \dontrun{
#' data(cell_to_metacell_pbmc_example)
#' project_counts_on_mc(sc_counts, cell_to_metacell_pbmc_example)
#' }
#'
#' @export
project_counts_on_mc <- function(sc_counts, cell_to_metacell, num_cores = parallel::detectCores()) {
    cell_to_metacell <- deframe(cell_to_metacell)
    cell_to_metacell <- cell_to_metacell[cell_to_metacell %!in% c(-2, -1, 0)] # make sure we don't have any outlier metacells
    assert_that(all(names(cell_to_metacell) %in% sc_counts@cell_names))


    removed_cells <- setdiff(sc_counts@cell_names, names(cell_to_metacell))
    if (length(removed_cells) > 0) {
        cli_alert_info("{.val {length(removed_cells)}} cells (out of {.val {length(sc_counts@cell_names)}}) do not have a metacell and have been removed.")
    }

    cell_to_metacell <- cell_to_metacell[intersect(names(cell_to_metacell), sc_counts@cell_names)]

    doMC::registerDoMC(num_cores)
    new_data <- plyr::llply(sc_counts@data, function(sc_mat) {
        cells <- intersect(colnames(sc_mat), names(cell_to_metacell))
        mc_mat <- sparse_matrix_tapply_sum(sc_mat[, cells], cell_to_metacell[cells])
        return(mc_mat)
    }, .parallel = TRUE)

    res <- new("McCounts", data = new_data, cell_names = as.character(sort(unique(cell_to_metacell))), genome = sc_counts@genome, genomic_bins = sc_counts@genomic_bins, id = sc_counts@id, description = sc_counts@description, path = sc_counts@path, cell_to_metacell = enframe(cell_to_metacell, "cell_id", "metacell"))

    return(res)
}

#' Read an McCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{write_sc_counts_from_bam})
#'
#' @examples
#' \dontrun{
#' mc_counts <- mcc_read("pbmc_reads_mc")
#' }
#'
#' @inheritParams scc_read
#' @export
mcc_read <- function(path, id = NULL, description = NULL) {
    sc_counts <- scc_read(path, id, description)
    md_file <- file.path(path, "metadata.yaml")
    md <- yaml::read_yaml(md_file)

    if (is.null(md$cell_to_metacell)) {
        cli_abort("Directory {.file {path}} does not contain a valid McCounts object (the metadata file {.file {md_file}} is missing the {.field cell_to_metacell} field.)")
    }

    mc_counts <- new("McCounts", data = sc_counts@data, cell_names = sc_counts@cell_names, genome = sc_counts@genome, genomic_bins = sc_counts@genomic_bins, id = sc_counts@id, description = sc_counts@description, path = sc_counts@path, cell_to_metacell = as_tibble(md$cell_to_metacell))
    return(mc_counts)
}

#' Extract McCounts to a tidy data frame
#'
#' @param mc_counts A McCounts object
#' @param metacells metacells to extract.
#' @param num_cores The number of cores to use.
#'
#' @return an intervals set (tibble) with additional "metacell" and "value" columns
#'
#' @examples
#' \dontrun{
#' df <- mcc_extract_to_df(mc_counts, c("1", "2"))
#' }
#'
#' @export
mcc_extract_to_df <- function(mc_counts, metacells = NULL, num_cores = parallel::detectCores()) {
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)
    extract_bin <- function(mat, bin, metacells) {
        Matrix::summary(mat[, metacells, drop = FALSE]) %>%
            mutate(
                start = i + bin$start,
                end = start + 1,
                chrom = bin$chrom,
                metacell = mc_counts@cell_names[j]
            ) %>%
            select(chrom, start, end, metacell, value = x) %>%
            as_tibble()
    }

    doMC::registerDoMC(num_cores)
    mc_data <- plyr::adply(mc_counts@genomic_bins, 1, function(bin) {
        return(
            extract_bin(mc_counts@data[[bin$name]], bin, metacells)
        )
    }, .parallel = TRUE)

    mc_data <- mc_data %>%
        select(chrom, start, end, metacell, value) %>%
        as_tibble()

    return(mc_data)
}

#' Create a misha track for each metacell in a McCounts object
#'
#' @description This would create a sparse track of the form "{track_prefix}.mc{metacell}" for each metacell in the McCounts object.
#'
#' @param mc_counts A McCounts object
#' @param track_prefix The prefix of the tracks to create. Track names will be of the form "{track_prefix}.mc{metacell}"
#' @param metacells metacells for which to create tracks. If NULL, all metacells will be used.
#' @param num_cores The number of cores to use.
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' mcc_to_tracks(mc_counts, "pbmc_mc")
#' }
#'
#' @export
mcc_to_tracks <- function(mc_counts, track_prefix, metacells = NULL, num_cores = parallel::detectCores()) {
    metacells <- metacells %||% mc_counts@cell_names
    metacells <- as.character(metacells)

    cli_alert("Extracting metacell data")
    d <- mcc_extract_to_df(mc_counts, metacells, num_cores)

    misha.ext::gtrack.create_dirs(track_prefix, showWarnings = FALSE)

    doMC::registerDoMC(num_cores)
    plyr::l_ply(mc_counts@cell_names, function(metacell) {
        track <- glue("{track_prefix}.mc{metacell}")
        cli_alert("Creating {track} track")
        x <- d %>% filter(metacell == !!metacell)
        gtrack.create_sparse(
            track = track,
            description = "",
            intervals = x %>% select(chrom, start, end),
            values = x$value
        )
    }, .parallel = TRUE)

    gdb.reload()

    cli_alert_success("Created {length(metacells)} tracks at {track_prefix}")
}
