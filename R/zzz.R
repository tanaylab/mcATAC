.onLoad <- function(libname, pkgname) {
    set_parallel()
    if (!exists(".mcatac_cache__", envir = globalenv())) {
        .mcatac_cache__ <<- new.env()
    }
}
