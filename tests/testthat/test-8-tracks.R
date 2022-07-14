mct <- mct_create(
    genome = "hg38",
    track_prefix = "pbmc_mc",
    id = "pbmc",
    description = "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)",
    metadata = mcmd,
    resolution = 20,
    window_size = 10
)


test_that("add rna works", {
    mct_rna <- add_mc_rna(mct, rna_mc_mat)
    egc <- t(t(rna_mc_mat) / colSums(rna_mc_mat))
    both_mcs <- intersect(mct_rna@metacells, colnames(egc))
    expect_equal(mct_rna@rna_egc, egc[, both_mcs])
})

test_that("add metadata works", {
    mct_md <- add_mc_metadata(mct, mcmd)
    expect_true(all(mct_md@metacells == mcmd$metacell))
})
