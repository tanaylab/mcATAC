filter_features_ok <- function(atac_sc, atac_sc_f) {
    peak_len <- atac_sc@peaks$end - atac_sc@peaks$start
    rmx <- sparseMatrixStats::rowMaxs(atac_sc@mat)
    low_max_peaks <- atac_sc@peaks$peak_name[rmx < 3]
    too_short_peaks <- atac_sc@peaks$peak_name[peak_len < 200]
    too_long_peaks <- atac_sc@peaks$peak_name[peak_len > 1000]
    peaks_to_remove <- unique(c(low_max_peaks, too_short_peaks, too_long_peaks))
    expect_equal(nrow(atac_sc_f@ignore_pmat), length(peaks_to_remove))
    expect_equal(nrow(atac_sc_f@ignore_peaks), length(peaks_to_remove))
    expect_true(all(atac_sc_f@ignore_peaks$peak_name %in% peaks_to_remove))
    expect_true(!any(atac_sc_f@peaks$peak_name %in% peaks_to_remove))
    expect_equal(nrow(atac_sc_f@peaks), nrow(atac_sc@peaks) - length(peaks_to_remove))
}

test_that("filter_features works", {
    atac_sc_f <- filter_features(atac_sc, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
    filter_features_ok(atac_sc, atac_sc_f)
})

test_that("filter_features doesn't remove additional fields from peaks", {
    atac_sc1 <- atac_sc
    atac_sc1@peaks <- atac_sc1@peaks %>% mutate(savta = 3)
    atac_sc_f <- filter_features(atac_sc1, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
    filter_features_ok(atac_sc, atac_sc_f)
    expect_equal(colnames(atac_sc_f@peaks), colnames(atac_sc1@peaks))
})

test_that("export works when there are ignored peaks", {
    atac_sc_f <- filter_features(scatac = atac_sc, minimal_max_umi = 3, min_peak_length = 200, max_peak_length = 1000)
    export_to_h5ad(atac_sc_f, fs::path(raw_dir, "atac_sc.h5ad"), compression = "gzip")
    atac_sc1 <- import_from_h5ad(fs::path(raw_dir, "atac_sc.h5ad"))
    expect_true(Matrix::mean(abs(atac_sc_f@mat - atac_sc1@mat)) <= 1e-9)
    expect_true(Matrix::mean(abs(atac_sc_f@ignore_pmat - atac_sc1@ignore_pmat)) <= 1e-9)
    expect_equal(atac_sc_f@peaks, atac_sc1@peaks, ignore_attr = TRUE)
    expect_equal(atac_sc_f@ignore_peaks, atac_sc1@ignore_peaks, ignore_attr = TRUE)
})

test_that("find_blacklist_overlaps works", {
    blacklist_overlaps <- find_blacklist_overlaps(atac_sc)
    atac_sc_f <- atac_ignore_peaks(atac_sc, blacklist_overlaps)
    misha.ext::gset_genome("hg38")
    expect_equal(
        atac_sc@peaks %>%
            as.data.frame() %>%
            misha.ext::gintervals.neighbors1("ENCODE.blacklist") %>%
            filter(dist == 0) %>%
            anti_join(atac_sc_f@ignore_peaks, by = c("chrom", "start", "end", "peak_name")) %>%
            nrow(),
        0
    )

    expect_equal(
        atac_sc_f@ignore_peaks,
        atac_sc@peaks %>%
            as.data.frame() %>%
            misha.ext::gintervals.neighbors1("ENCODE.blacklist") %>%
            filter(dist == 0) %>%
            select(any_of(colnames(atac_sc_f@peaks))),
        ignore_attr = TRUE
    )
})

test_that("atac_ignore_peaks works twice when reset=FALSE", {
    set.seed(60427)
    ignore_peaks1 <- atac_sc@peaks %>% sample_n(100)
    ignore_peaks2 <- atac_sc@peaks %>%
        anti_join(ignore_peaks1, by = c("chrom", "start", "end", "peak_name")) %>%
        sample_n(100)

    atac_sc1 <- atac_ignore_peaks(atac_sc, ignore_peaks1, reset = FALSE)
    atac_sc2 <- atac_ignore_peaks(atac_sc1, ignore_peaks2, reset = FALSE)

    expect_setequal(atac_sc2@ignore_peaks$peak_name, unique(c(ignore_peaks1$peak_name, ignore_peaks2$peak_name)))
    expect_setequal(rownames(atac_sc2@ignore_pmat), unique(c(ignore_peaks1$peak_name, ignore_peaks2$peak_name)))
    expect_true(all(atac_sc2@ignore_peaks$peak_name %!in% atac_sc2@peaks$peak_name))
    expect_true(all(atac_sc1@ignore_peaks$peak_name %in% atac_sc2@ignore_peaks$peak_name))
})

test_that("atac_ignore_peaks works twice when reset=TRUE", {
    set.seed(60427)
    ignore_peaks1 <- atac_sc@peaks %>% sample_n(100)
    ignore_peaks2 <- atac_sc@peaks %>%
        anti_join(ignore_peaks1, by = c("chrom", "start", "end", "peak_name")) %>%
        sample_n(100)

    atac_sc1 <- atac_ignore_peaks(atac_sc, ignore_peaks1, reset = TRUE)
    atac_sc2 <- atac_ignore_peaks(atac_sc1, ignore_peaks2, reset = TRUE)

    expect_setequal(atac_sc2@ignore_peaks$peak_name, ignore_peaks2$peak_name)
    expect_true(!any(atac_sc2@ignore_peaks$peak_name %in% atac_sc1@ignore_peaks$peak_name))
})
