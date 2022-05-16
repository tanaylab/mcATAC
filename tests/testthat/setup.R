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
atac_sc <<- import_from_10x(raw_dir, "hg38", id = "pbmc", description = "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)")

data(mcmd)
atac_mc <<- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example, metadata = mcmd, min_int_frac = 0.5, mc_size_eps_q = 0.1)

data(rna_mc_mat)

# Make atac_sc and atac_mc const in order to test them independetly
lockBinding("atac_sc", globalenv())
lockBinding("atac_mc", globalenv())


samtools_bin <- "samtools"

withr::defer(
    {
        unlink(raw_dir, recursive = TRUE)
        unlockBinding("atac_sc", globalenv())
        unlockBinding("atac_mc", globalenv())
    },
    teardown_env()
)
