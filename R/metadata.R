
#' Add per-metacell metadata to an McATAC object
#'
#' @param mcatac an McATAC object
#' @param metadata data frame with a column called 'metacell' and additional metacell annotations, or the name of a delimited file which contains such annotations.
#'
#' @examples
#' \dontrun{
#' data(mcmd)
#' mc_atac <- add_metadata(mc_atac, mcmd)
#' }
#'
#' @export
add_mc_metadata <- function(mcatac, metadata) {
    add_metadata(mcatac, metadata, "metacell")
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

        # make sure that all cells/metacells exist within the matrix
        missing_cells <- metadata[[metadata_id_field]] %!in% colnames(obj@mat)
        if (any(missing_cells)) {
            missing_cells <- paste(unique(metadata[[metadata_id_field]][missing_cells]), collapse = ", ")
            cli_abort("The following {metadata_id_field}s are missing from {.field mat} colnames: {.val {missing_cells}}")
        }
    }

    obj@metadata <- metadata

    return(obj)
}


get_cell_type_colors <- function(metadata) {
    metadata %>%
        distinct(cell_type, color) %>%
        deframe()
}

#' Does the McATAC object contain per-metacell cell type annotation
#'
#' @param mc_atac an McATAC object
#'
#' @return TRUE if the McATAC object contains per-metacell cell type annotation, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_cell_type(atac_mc)
#' }
#'
#' @export
has_cell_type <- function(mc_atac) {
    assert_atac_object(mc_atac)
    return(!is.null(mc_atac@metadata) && !is.null(mc_atac@metadata$cell_type))
}

#' Does the McATAC object contain cell type color annotation
#'
#' @param mc_atac an McATAC object
#'
#' @return TRUE if the McATAC object contains cell type color annotation, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' has_cell_type_color(atac_mc)
#' }
#'
#' @export
has_cell_type_colors <- function(mc_atac) {
    assert_atac_object(mc_atac)
    return(has_cell_type(mc_atac) && !is.null(mc_atac@metadata$color))
}
