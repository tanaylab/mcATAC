#' Add per-metacell metadata to an McPeaks object
#'
#' @param atac_mc an McPeaks object
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#'
#' @examples
#' \dontrun{
#' data(mcmd)
#' atac_mc <- add_metadata(atac_mc, mcmd)
#' }
#'
#' @export
add_mc_metadata <- function(atac_mc, metadata) {
    add_metadata(atac_mc, metadata, "metacell")
}

add_metadata <- function(obj, metadata, metadata_id_field) {
    if (!is.null(metadata)) {
        if (is.character(metadata)) {
            metadata <- tgutil::fread(metadata)
        }
        if (!is.data.frame(metadata)) {
            cli_abort("{.field metadata} is not a data frame")
        }

        metadata <- as_tibble(metadata)

        if (!has_name(metadata, metadata_id_field)) {
            cli_abort("{.field metadata} doesn't have the required field {.field {metadata_id_field}}")
        }

        metadata <- as_tibble(metadata)

        if (methods::is(obj, "ATACPeaks")) {
            cell_names <- colnames(obj@mat)
        } else if (methods::is(obj, "ScCounts")) {
            cell_names <- obj@cell_names
        } else if (methods::is(obj, "McTracks")) {
            cell_names <- obj@metacells
        }

        # make sure that all cells/metacells exist within the matrix
        missing_cells <- metadata[[metadata_id_field]] %!in% cell_names
        if (any(missing_cells)) {
            missing_cells <- paste(unique(metadata[[metadata_id_field]][missing_cells]), collapse = ", ")
            cli_warn("The following {metadata_id_field}s are missing from {.field mat} colnames: {.val {missing_cells}}")
        }

        metadata[[metadata_id_field]] <- as.character(metadata[[metadata_id_field]])
    }

    obj@metadata <- metadata

    return(obj)
}


get_cell_type_colors <- function(metadata) {
    metadata %>%
        distinct(cell_type, color) %>%
        select(cell_type, color) %>%
        deframe()
}

get_metacell_colors <- function(metadata) {
    metadata %>%
        distinct(metacell, color) %>%
        select(metacell, color) %>%
        deframe()
}

#' Does the McPeaks object contain per-metacell cell type annotation
#'
#' @param atac_mc an McPeaks object
#'
#' @return TRUE if the McPeaks object contains per-metacell cell type annotation, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_cell_type(atac_mc)
#' }
#'
#' @export
has_cell_type <- function(atac_mc) {
    assert_atac_object(atac_mc, class = "ATAC")
    return(!is.null(atac_mc@metadata) && !is.null(atac_mc@metadata$cell_type))
}

#' Does the McPeaks object contain cell type color annotation
#'
#' @param atac_mc an McPeaks object
#'
#' @return TRUE if the McPeaks object contains cell type color annotation, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_cell_type_color(atac_mc)
#' }
#'
#' @export
has_cell_type_colors <- function(atac_mc) {
    assert_atac_object(atac_mc, class = "ATAC")
    return(has_cell_type(atac_mc) && !is.null(atac_mc@metadata$color))
}
