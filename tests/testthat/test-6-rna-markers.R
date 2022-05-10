test_that("get_rna_egc works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    rna_mat <- get_rna_egc(atac_mc, c("CD4", "CD8A"), epsilon = 0)
    expect_equal(rna_mat, atac_mc@rna_egc[c("CD4", "CD8A"), ])
    rna_mat <- get_rna_egc(atac_mc, c("CD4", "CD8A"), epsilon = 1e-5)
    expect_equal(rna_mat, atac_mc@rna_egc[c("CD4", "CD8A"), ] + 1e-5)
    rna_mat <- get_rna_egc(atac_mc, "CD4", epsilon = 0)
    expect_equal(rna_mat, atac_mc@rna_egc["CD4", , drop = FALSE])
    expect_error(get_rna_egc(atac_mc, "savta"))

    rna_mat <- get_rna_egc(atac_mc, c("CD4", "CD8A", "MIR1302-2HG", "FAM138A"), epsilon = 0)
    expect_equal(rna_mat, atac_mc@rna_egc[c("CD4", "CD8A"), ])

    rna_mat <- get_rna_egc(atac_mc, c("CD4", "CD8A", "MIR1302-2HG", "FAM138A"), rm_zeros = FALSE, epsilon = 0)
    expect_equal(rna_mat, atac_mc@rna_egc[c("CD4", "CD8A", "MIR1302-2HG", "FAM138A"), ])
})

test_that("get_rna_fp works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    rna_fp <- get_rna_fp(atac_mc, c("CD4", "CD8A"))
    rna_egc <- get_rna_egc(atac_mc, c("CD4", "CD8A"))
    rna_fp1 <- rna_egc / apply(rna_egc, 1, median, na.rm = TRUE)
    expect_equal(rna_fp, rna_fp1)
})

test_that("get_rna_markers works", {
    atac_mc <- add_mc_rna(atac_mc, rna_mc_mat)
    expect_equal(get_rna_markers(atac_mc, 2), c("CA6", "CACHD1"))
})
