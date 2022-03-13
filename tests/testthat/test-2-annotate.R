test_that("annotate_peaks works", {
    mc_atac1 <- annotate_peaks(mc_atac)
    expect_equal(mc_atac$peaks, mc_atac1$peaks %>% dplyr::select(-peak_annot))
})
