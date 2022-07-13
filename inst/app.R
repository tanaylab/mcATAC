if (dir.exists("code")) {
    pkgload::load_all("code", export_all = FALSE)
}
options(gmultitasking = FALSE)
mct <- readr::read_rds("mct.rds")
mcATAC::run_app(mct)
