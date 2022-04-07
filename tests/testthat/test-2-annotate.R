test_that("annotate_peaks works", {
    atac_mc1 <- annotate_peaks(atac_mc)
    expect_equal(atac_mc@peaks, atac_mc1@peaks %>% select(-peak_type, -closest_tss))
})
