#' Negation of the %in% operator
#'
#' @export
#' @noRd
`%!in%` <- Negate(`%in%`)


check_files_exist <- function(files) {
    for (file in files) {
        if (!file.exists(file) && !grepl("^http", file)) {
            cli_abort("{.file {file}} file doesn't exist.", call = parent.frame(1))
        }
    }
}

is_sparse_matrix <- function(mat) {
    return(methods::is(mat, "sparseMatrix"))
}
