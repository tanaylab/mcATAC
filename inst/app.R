if (dir.exists("code")) {
    pkgload::load_all("code", export_all = FALSE)
}
options(gmultitasking = FALSE)
mct <- readr::read_rds("mct.rds")
if (file.exists("hc.rds")){
    hc <- readr::read_rds("hc.rds")
} else {
    hc <- NULL
}
mcATAC::run_app(mct, hc = hc)
