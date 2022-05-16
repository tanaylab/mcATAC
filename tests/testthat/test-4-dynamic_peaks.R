test_that("identify_dynamic_peaks works with bmq", {
    fig_file <- tempfile(fileext = ".png")
    png(fig_file)
    dynamic_peaks <- identify_dynamic_peaks(atac_mc, method = "bmq", mean_thresh_q = 0.1, cov_q_thresh = 0.6, num_bins = 100)
    dev.off()
    expect_true(file.exists(fig_file))
    expect_true(all(rlang::has_name(dynamic_peaks, c("chrom", "start", "end", "peak_name"))))
    expect_true(all(dynamic_peaks$peak_name %in% atac_mc@peaks$peak_name))
})

test_that("identify_dynamic_peaks works with gmm", {
    fig_file <- tempfile(fileext = ".png")
    png(fig_file)
    dynamic_peaks <- identify_dynamic_peaks(atac_mc, method = "gmm")
    dev.off()
    expect_true(file.exists(fig_file))
    expect_true(all(rlang::has_name(dynamic_peaks, c("chrom", "start", "end", "peak_name"))))
    expect_true(all(dynamic_peaks$peak_name %in% atac_mc@peaks$peak_name))
})

test_that("identify_dynamic_peaks works when plot = FALSE", {
    fig_file <- tempfile(fileext = ".png")
    png(fig_file)
    dynamic_peaks <- identify_dynamic_peaks(atac_mc, method = "bmq", mean_thresh_q = 0.1, cov_q_thresh = 0.6, num_bins = 100, plot = FALSE)
    dev.off()
    expect_true(!file.exists(fig_file))
    expect_true(all(rlang::has_name(dynamic_peaks, c("chrom", "start", "end", "peak_name"))))
    expect_true(all(dynamic_peaks$peak_name %in% atac_mc@peaks$peak_name))
})

test_that("identify_dynamic_peaks errors with unknown methods", {
    expect_error(identify_dynamic_peaks(atac_mc, mean_thresh_q = 0.05, method = "savta"))
})
