set.seed(60427)
library(dplyr)
library(glue)

test_dir <- tempdir()

raw_dir <- fs::path(test_dir, "pbmc_data")

fs::dir_create(raw_dir)

# set the timeout to more than 60 seconds in order to be able to download the files on github
# options(timeout = 1e4)
# download_pbmc_example_data()

system(glue("cp -rf /net/mraid14/export/tgdata/users/aviezerl/src/mcATAC/pbmc_data/* {raw_dir}/"))

data(cell_to_metacell_pbmc_example)
atac_sc <- import_from_10x(raw_dir, "hg38")
atac_mc <- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example)

withr::defer(
    {
        unlink(raw_dir)
    },
    teardown_env()
)
