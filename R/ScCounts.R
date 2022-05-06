#' ScCounts object
#'
#'
#' @description ScCounts is a data structure that contains the counts of reads for a set of cells in a set of genomic bins.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot data a named list with sparse matrices for each chromosome.
#' @slot path original path from which the object was loaded (optional)
#' @slot genomic_bins an intervals set with the genomic bins. Should have columns named "chrom", "start", "end" and "name" with the chromosome, start and end positions and name of the bin.
#'
#' @exportClass ScCounts
ScCounts <- setClass(
    "ScCounts",
    slots = c(
        id = "character",
        description = "character",
        data = "list",
        genome = "character",
        cell_names = "character",
        genomic_bins = "data.frame",
        path = "character"
    )
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
        .Object@path <- path
        validate_ScCounts(.Object)
        return(.Object)
    }
)

validate_ScCounts <- function(obj) {
    # make sure genomic bins exist in the genome
    # make sure files of genomic bins exist
    # make sure that the dimensions of the matrices are correct
    # make sure that the cell names are correct
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

#' Write a counts object to a file.
#'
#' @param object a counts object (ScCounts or McCounts)
#' @param out_dir the output directory
#' @param num_cores the number of cores to use
#'
#' @examples
#' \dontrun{
#' scc_write(sc_counts, "sc_counts")
#' }
#'
#' @export
scc_write <- function(object, out_dir, num_cores = parallel::detectCores()) {
    data_dir <- file.path(out_dir, "data")
    if (!dir.exists(data_dir)) {
        dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
        cli_alert("Created directory {.file {out_dir}}")
    }

    doMC::registerDoMC(num_cores)
    plyr::l_ply(names(object@data), function(region) {
        print(region)
        out_file <- file.path(data_dir, glue("{region}.mtx"))
        tgutil::fwrite_mm(object@data[[region]], out_file)
        system(glue("gzip {out_file} && rm -f {out_file}"))
    }, .parallel = TRUE)

    cli_alert_info("Written sparse matrices for {.val {nrow(object@genomic_bins)}} genomic bins")

    counts_md <- list(
        cell_names = object@cell_names,
        genome = object@genome,
        id = object@id,
        description = object@description,
        data_dir = "data",
        genomic_bins = object@genomic_bins
    )

    if (class(object) == "McCounts") {
        counts_md$cell_to_metacell <- object@cell_to_metacell
    }

    yaml::write_yaml(counts_md, file = file.path(out_dir, "metadata.yaml"))

    cli_alert_success("Written object to {.file {out_dir}}")
}

#' Read an ScCounts object from a directory
#'
#' @param path path to the directory containing the object (which was created by \code{write_sc_counts})
#' @param id an identifier for the object (optional)
#' @param description description of the object (optional)
#'
#' @return an ScCounts object
#'
#' @examples
#' \dontrun{
#' counts <- scc_read("pbmc_reads")
#' }
#'
#' @export
scc_read <- function(path, id = NULL, description = NULL) {
    md_file <- file.path(path, "metadata.yaml")
    if (!file.exists(md_file)) {
        cli_abort("Directory {.file {path}} does not contain a valid ScCounts object")
    }
    md <- yaml::read_yaml(md_file)
    data_dir <- file.path(path, md$data_dir)
    required_fields <- c("cell_names", "genome", "id", "description", "data_dir", "genomic_bins")
    purrr::walk(required_fields, ~ {
        if (.x %!in% names(md)) {
            cli_abort("Directory {.file {path}} does not contain a valid ScCounts object (the metadata file {.file {md_file}} is missing the {.field {.x}} field.")
        }
    })

    genomic_bins <- as_tibble(md$genomic_bins)
    data <- read_sc_counts_data(data_dir, genomic_bins, md$cell_names, md$genome)

    counts <- new("ScCounts", data = data, cell_names = md$cell_names, genome = md$genome, genomic_bins = genomic_bins, id = id %||% md$id, description = description %||% md$description, path = path)
    return(counts)
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

read_sc_counts_data <- function(data_dir, genomic_bins, cell_names, genome, num_cores = parallel::detectCores()) {
    gset_genome(genome)
    intervals <- gintervals.all()

    # make sure the scope is OK

    if (any(setdiff(unique(genomic_bins$chrom), intervals$chrom))) {
        missing_chroms <- setdiff(unique(genomic_bins$chrom), intervals$chrom)
        cli_abort("The following chromosomes are not present in the genome: {.val {missing_chroms}}")
    }

    # chromosomes <- intersect(chromosomes, intervals$chrom)

    # again: make sure the scope is OK

    doMC::registerDoMC(num_cores)
    data <- plyr::llply(genomic_bins$name, function(region) {
        region_file <- file.path(data_dir, paste0(region, ".mtx.gz"))
        if (!file.exists(region_file)) {
            cli_abort("File {.file {region_file}} does not exist")
        }

        m <- tgutil::fread_mm(region_file)
        colnames(m) <- cell_names
        return(m)
    }, .parallel = TRUE)

    names(data) <- genomic_bins$name

    return(data)
}
