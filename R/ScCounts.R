#' ScCounts object
#'
#'
#' @description ScCounts is a data structure that contains the counts of reads for a set of cells in a set of chromosomes.
#'
#' @slot id an identifier for the object, e.g. "pbmc".
#' @slot description description of the object, e.g. "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)"
#' @slot genome genome assembly of the peaks. e.g. "hg38", "hg19", "mm9", "mm10"
#' @slot data a named list with sparse matrices for each chromosome.
#' @slot path original path from which the object was loaded (optional)
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
        chromosomes = "character",
        path = "character"
    )
)

setMethod(
    "initialize",
    signature = "ScCounts",
    definition = function(.Object, data, cell_names, genome, chromosomes = NULL, id = "", description = "", path = "") {
        gset_genome(genome)
        if (is.null(chromosomes)) {
            chromosomes <- as.character(unique(gintervals.all()$chrom))
        }
        .Object@id <- id
        .Object@description <- description
        .Object@data <- data
        .Object@genome <- genome
        .Object@chromosomes <- chromosomes
        .Object@cell_names <- cell_names
        .Object@path <- path
        validate_ScCounts(.Object)
        return(.Object)
    }
)

validate_ScCounts <- function(obj) {
    if (any(obj@chromosomes %!in% names(obj@data))) {
        missing_chroms <- obj@chromosomes[obj@chromosomes %!in% names(obj@data)]
        cli_alert_warning("The following chromosome(s) were not found in the {.field ScCounts} object data: {.val {missing_chroms}}")
    }
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
    cli::cli_text("{.cls {object_type}} object with {.val {length(object@chromosomes)}} chromosomes (of {.val {object@genome}}) and {.val {length(object@cell_names)}} {column_type}s.")
    if (object@id != "") {
        cli::cli_text(c("id: {.val {object@id}}"))
    }
    if (object@description != "") {
        cli::cli_text(c("description: {.val {object@description}}"))
    }
    if (object@path != "") {
        cli::cli_text(c("Loaded from: {.file {object@path}}"))
    }

    cli::cli_text("Slots include:")
    cli_ul(c("{.code @data}: a named list of sparse matrices where rows are genomic coordinates and columns are single cells."))
    cli_ul(c("{.code @cell_names}: names of the {column_type}s."))
    cli_ul(c("{.code @chromosomes}: names of the chromosomes."))
    cli_ul(c("{.code @genome}: genome assembly."))
}
