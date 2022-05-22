#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Set misha root based on genome name
#'
#' @name gset_genome
#' @importFrom misha.ext gset_genome
#' @inheritParams misha.ext::gset_genome
#' @noRd
#' @export gset_genome
NULL
