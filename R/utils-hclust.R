#' Generate an hclust object from a data.frame with cell_type, group and order columns
#'
#' @param df data.frame with cell_type, group and order
#'
#' @return hclust object with cell_type as labels and order as order of the leaves
#'
#' @examples
#' df <- data.frame(
#'     cell_type = c(
#'         "PGC", "Surface ectoderm", "Neural crest",
#'         "Neural tube/Floor plate", "Forebrain/Midbrain/Hindbrain", "Caudal neural plate",
#'         "Rostral neural plate", "Neural plate boundary", "Definitive ectoderm",
#'         "Epiblast", "Primitive streak", "Caudal epiblast", "Tail bud - neural",
#'         "Tail bud - mesoderm", "Early nascent mesoderm"
#'     ), order = 1:15,
#'     group = c(
#'         "ecto", "ecto", "ecto", "ecto", "ecto", "ecto",
#'         "ecto", "ecto", "ecto", "early", "early", "early", "early",
#'         "Meso", "Meso"
#'     )
#' )
#'
#' group_to_hclust(df)
#'
#' @export
group_to_hclust <- function(df) {
    # Create an empty matrix of 1s
    distance_matrix <- matrix(1, nrow = nrow(df), ncol = nrow(df))
    rownames(distance_matrix) <- df$cell_type
    colnames(distance_matrix) <- df$cell_type

    # set the distance between elements the same group to 0 (using df$group)
    for (grp in unique(df$group)) {
        distance_matrix[df$group == grp, df$group == grp] <- 0
    }

    # Hierarchical clustering
    hc <- hclust(as.dist(distance_matrix))

    # Label leaves with cell_type
    hc$labels <- df$cell_type

    hc <- as.hclust(reorder(as.dendrogram(hc), df$order, agglo.FUN = mean))

    return(hc)
}

get_heatmap_idx <- function(x, y, m) {
    n_x <- nrow(m)
    n_y <- ncol(m)
    x_idx <- floor(x * n_x) + 1

    y_idx <- n_y - max(floor(y * n_y), 0)
    return(c(x_idx, y_idx))
}
