#' ScCounts object
#'
#'
#' @description ScCounts is a data structure that contains the counts of reads for a set of cells in a set of genomic bins.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot data a named list with sparse matrices for each chromosome. Each sparse matrix at position i,j contains the number of reads that map to the cell j in the genomic coordinate i (**in 1 based format**).
#' @slot path original path from which the object was loaded (optional)
#' @slot genomic_bins an intervals set with the genomic bins. Should have columns named "chrom", "start", "end" and "name" with the chromosome, start and end positions and name of the bin. The names of the data matrices are the "name" field. Note that even though the data is stored in
#' 1-based format, the genomic bins are in 0-based format (misha format), hence the filenames are in 0-based format.
#'
#' @exportClass ScCounts
ScCounts <- setClass(
    "ScCounts",
    slots = c(
        data = "list",
        cell_names = "character",
        genomic_bins = "data.frame"
    ),
    contains = "ATAC"
)

setMethod(
    "initialize",
    signature = "ScCounts",
    definition = function(.Object, data, cell_names, genome, genomic_bins, id = "", description = "", path = "") {
        gset_genome(genome)
        .Object@id <- id
        .Object@description <- description
        .Object@data <- data
        .Object@genome <- genome
        .Object@genomic_bins <- genomic_bins
        .Object@cell_names <- cell_names
        if (path != "") {
            path <- as.character(normalizePath(path))
        }
        .Object@path <- path
        validate_ScCounts(.Object)
        return(.Object)
    }
)

validate_ScCounts <- function(obj) {
    gset_genome(obj@genome)
    if (!all(has_name(obj@genomic_bins, c("chrom", "start", "end", "name")))) {
        cli_abort("genomic bins should have columns named 'chrom', 'start', 'end' and 'name'")
    }

    if (!intervals_in_genome(obj@genomic_bins)) {
        cli_abort("some genomic bins are not within the genome")
    }

    withr::with_options(list(scipen = 1e5), {
        if (!all(obj@genomic_bins$name == paste(obj@genomic_bins$chrom, obj@genomic_bins$start, obj@genomic_bins$end, sep = "_"))) {
            cli_abort("genomic bins should have a name that is the concatenation of the chromosome, start and end position")
        }
    })

    if (!all(names(obj@data) == obj@genomic_bins$name)) {
        cli_abort("data be a list with the same names as the genomic bins")
    }

    purrr::walk2(obj@data, 1:length(obj@data), ~ {
        if (!is_sparse_matrix(.x)) {
            cli_abort("data should be a list with sparse matrices")
        }
        if (ncol(.x) != length(obj@cell_names)) {
            cli_abort("data should have the same number of columns as the number of cell names")
        }
        if (!all(colnames(.x) == obj@cell_names)) {
            cli_abort("data matrices should have the same columns as the cell names")
        }
        intervs <- obj@genomic_bins %>% slice(.y)
        if (nrow(.x) != (intervs$end - intervs$start)) {
            cli_abort("data should have the same number of rows as coordintes in the genomic bin")
        }
    })
}

#' @export
#' @noRd
setMethod(
    "show",
    signature = "ScCounts",
    definition = function(object) {
        print_counts_object(object, "ScCounts", "cell")
    }
)

print_counts_object <- function(object, object_type, column_type) {
    cli::cli_text("{.cls {object_type}} object with {.val {nrow(object@genomic_bins)}} genomic bins (of {.val {object@genome}}) and {.val {length(object@cell_names)}} {column_type}s.")
    if (length(object@id) > 0 && object@id != "") {
        cli::cli_text(c("id: {.val {object@id}}"))
    }
    if (length(object@description) > 0 && object@description != "") {
        cli::cli_text(c("description: {.val {object@description}}"))
    }
    if (length(object@path) > 0 && object@path != "") {
        cli::cli_text(c("Loaded from: {.file {object@path}}"))
    }

    cli::cli_text("Slots include:")
    cli_ul(c("{.code @data}: a named list of sparse matrices where rows are genomic coordinates and columns are {column_type}s."))
    cli_ul(c("{.code @cell_names}: names of the {column_type}s."))
    cli_ul(c("{.code @genomic_bins}: intervals set with the genomic bins."))
    cli_ul(c("{.code @genome}: genome assembly."))
}

#' Write a counts object to a file
#'
#' @param object a counts object (ScCounts or McCounts)
#' @param out_dir the output directory
#' @param overwrite whether to overwrite existing directory
#'
#' @examples
#' \dontrun{
#' scc_write(sc_counts, "sc_counts")
#' mcc_write(mc_counts, "mc_counts")
#' }
#'
#' @export
scc_write <- function(object, out_dir, overwrite = FALSE) {
    data_dir <- file.path(out_dir, "data")

    if (dir.exists(out_dir)) {
        if (!overwrite) {
            cli_abort("Output directory {.file {out_dir}} already exists. Use 'overwrite = TRUE' to overwrite.")
        } else {
            cli_warn("Output directory {.file {out_dir}} already exists. Overwriting.")
            unlink(out_dir, recursive = TRUE)
        }
    }
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    cli_alert("Writing to {.file {out_dir}}")

    plyr::l_ply(names(object@data), function(region) {
        out_file <- file.path(data_dir, glue("{region}.mtx"))
        tgutil::fwrite_mm(object@data[[region]], out_file)
        system(glue("gzip {out_file} && rm -f {out_file}"))
    }, .parallel = getOption("mcatac.parallel"))

    cli_alert_info("Written sparse matrices for {.val {nrow(object@genomic_bins)}} genomic bins")

    counts_md <- list(
        cell_names = object@cell_names,
        genome = object@genome,
        id = object@id,
        description = object@description,
        data_dir = "data",
        genomic_bins = object@genomic_bins$name
    )

    if (class(object) == "McCounts") {
        counts_md$cell_to_metacell <- object@cell_to_metacell
    }

    yaml::write_yaml(counts_md, file = file.path(out_dir, "metadata.yaml"))

    cli_alert_success("Written object to {.file {out_dir}}")
}

#' @rdname scc_write
#' @export
mcc_write <- function(object, out_dir, overwrite = FALSE) {
    scc_write(object, out_dir, overwrite = overwrite)
}

#' Read an ScCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{write_sc_counts_from_fragments/bam})
#' @param id an identifier for the object (optional)
#' @param description description of the object (optional)
#' @param verbose print success message (default: TRUE)
#'
#' @return a ScCounts object
#'
#' @examples
#' \dontrun{
#' counts <- scc_read("pbmc_reads")
#' }
#'
#' @export
scc_read <- function(path, id = NULL, description = NULL, verbose = TRUE) {
    read_counts_object(path, "ScCounts", id = id, description = description, verbose = verbose)
}

read_counts_object <- function(path, type, id = NULL, description = NULL, verbose = TRUE) {
    md_file <- file.path(path, "metadata.yaml")
    if (!file.exists(md_file)) {
        cli_abort("Directory {.file {path}} does not contain a valid ScCounts object")
    }
    md <- yaml::read_yaml(md_file)
    if (type == "ScCounts" && !is.null(md$cell_to_metacell)) {
        cli_abort("Directory {.file {path}} contains a McCounts object. Please use {.code mcc_read} instead.")
    }

    data_dir <- file.path(path, md$data_dir)
    required_fields <- c("cell_names", "genome", "id", "description", "data_dir", "genomic_bins")
    purrr::walk(required_fields, ~ {
        if (.x %!in% names(md)) {
            cli_abort("Directory {.file {path}} does not contain a valid ScCounts object (the metadata file {.file {md_file}} is missing the {.field {.x}} field.")
        }
    })

    genomic_bins <- parse_genomic_bin_names(md$genomic_bins, md$genome)

    data <- read_sc_counts_data(data_dir, genomic_bins, md$cell_names, md$genome)

    counts <- new("ScCounts", data = data, cell_names = md$cell_names, genome = md$genome, genomic_bins = genomic_bins, id = id %||% md$id, description = description %||% md$description, path = path)

    if (verbose) {
        cli_alert_success("Succesfully read a ScCounts object from {.file {path}}")
    }

    return(counts)
}


#' Covert genomic bin names to intervals
#'
#' @param bin_names genomic bin names in the format of "{chrom}_{start}_{end}"
#' @param genome name of the genome assembly
#'
#' @return an intervals set with an additional "name" field with the genomic bin name. Names that were not parsed correctly (e.g.
#'  because the chromosome name is not recognized) are removed.
#'
#' @examples
#' \dontrun{
#' parse_genomic_bin_names(c("chr1_0_100", "chr1_100_200", "chrUn_KI270747v1_0_100"), "hg38")
#' parse_genomic_bin_names(c("saba_0_100", "chr1_0_100"), "hg38")
#' }
#'
#' @noRd
parse_genomic_bin_names <- function(bin_names, genome) {
    gset_genome(genome)

    chroms <- gintervals.all()$chrom
    pattern <- glue("^({chr})_(.+)_(.+)$", chr = paste0(chroms, collapse = "|"))
    genomic_bins <- stringr::str_match_all(bin_names, pattern) %>%
        purrr::map_dfr(~ tibble(chrom = .x[2], start = .x[3], end = .x[4], name = .x[1])) %>%
        mutate(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end)) %>%
        select(chrom, start, end, name) %>%
        # make sure that the order remained
        mutate(name = factor(name, levels = bin_names)) %>%
        arrange(name) %>%
        mutate(name = as.character(name)) %>%
        na.omit()

    unknown_bins <- setdiff(bin_names, genomic_bins$name)
    if (length(unknown_bins) > 0) {
        cli_warn("The following bins were not found in the genome: {.val {unknown_bins}}. This usually means that the genome ({.val {genome}}) is not compatible with the bins.")
    }

    return(genomic_bins)
}

#' Read the sparse matrices of a ScCounts object from a directory
#'
#' @param data_dir path to the directory containing the sparse matrices
#' @param genomic_bins genomic_bins slot of the ScCounts object
#' @param cell_names cell_names slot of the ScCounts object
#' @param genome genome slot of the ScCounts object
#'
#' @noRd
read_sc_counts_data <- function(data_dir, genomic_bins, cell_names, genome) {
    gset_genome(genome)
    intervals <- gintervals.all()

    if (any(setdiff(unique(genomic_bins$chrom), intervals$chrom))) {
        missing_chroms <- setdiff(unique(genomic_bins$chrom), intervals$chrom)
        cli_abort("The following chromosomes are not present in the genome: {.val {missing_chroms}}")
    }

    # Make sure all genomic bins are within the genome boundries
    if (!intervals_in_genome(genomic_bins)) {
        cli_abort("Some genomic bins are outside the genome boundaries.")
    }

    data <- plyr::llply(genomic_bins$name, function(region) {
        region_file <- file.path(data_dir, paste0(region, ".mtx.gz"))
        if (!file.exists(region_file)) {
            cli_abort("File {.file {region_file}} does not exist")
        }

        m <- tgutil::fread_mm(region_file)
        colnames(m) <- cell_names
        return(m)
    }, .parallel = getOption("mcatac.parallel"))

    names(data) <- genomic_bins$name

    return(data)
}


#' Create an ScPeaks object from an ScCounts object
#'
#' @description given an ScCounts object and peaks, summarise the counts over the peaks and return a ScPeaks object
#'
#' @param sc_counts a ScCounts object
#' @param peaks a data frame with the peak intervals (chrom, start, end) and a column called "peak_name"
#' @param cells names of cells to include. Default: all cells.
#'
#' @inheritParams project_atac_on_mc
#'
#' @return a McPeaks object
#'
#' @examples
#' \dontrun{
#' atac_sc <- import_from_10x("pbmc_data", genome = "hg38")
#' sc_counts <- scc_read("pbmc_reads")
#' atac_sc <- scc_to_peaks(sc_counts, atac_sc@peaks)
#' }
#'
#' @export
scc_to_peaks <- function(sc_counts, peaks, cells = NULL, metadata = NULL, mc_size_eps_q = 0.1) {
    assert_atac_object(sc_counts, class = "ScCounts")
    cells <- cells %||% sc_counts@cell_names
    cells <- as.character(cells)
    if (!has_name(peaks, "peak_name")) {
        peaks <- peaks %>% mutate(peak_name = peak_names(.))
        cli_alert_warning("The {.var peaks} didn't have a column called {.field peak_name}, so it was added using the {.code peak_names} function.")
    }
    matrices <- plyr::alply(sc_counts@genomic_bins, 1, function(bin) {
        return(
            summarise_bin(sc_counts@data[[bin$name]], bin, peaks, cells)
        )
    }, .parallel = getOption("mcatac.parallel"))

    mat <- Reduce("+", matrices)

    sc_atac <- new("ScPeaks", mat = mat, peaks = peaks, genome = sc_counts@genome, id = sc_counts@id, description = sc_counts@description, metadata = metadata, path = sc_counts@path)

    cli_alert_success("Created a new ScPeaks object with {.val {ncol(sc_atac@mat)}} cells and {.val {nrow(sc_atac@mat)}} ATAC peaks.")

    return(sc_atac)
}
