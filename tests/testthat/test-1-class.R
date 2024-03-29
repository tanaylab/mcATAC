test_that("import from 10x works", {
    expect_equal(class(atac_sc), "ScPeaks", ignore_attr = TRUE)
    expect_equal(nrow(atac_sc@mat), 107829)
    expect_equal(ncol(atac_sc@mat), 11909)
    expect_equal(colnames(atac_sc@peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_sc@peaks), atac_sc@peaks$peak_name)
    expect_equal(atac_sc@id, "pbmc")
    expect_equal(atac_sc@tad_based, TRUE)
    expect_equal(atac_sc@description, "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)")
    expect_equal(atac_sc@path, normalizePath(file.path(raw_dir, "matrix.mtx")), ignore_attr = TRUE)
})

test_that("import from matrix works", {
    mat <- atac_sc@mat
    peaks <- atac_sc@peaks %>%
        select(chrom, start, end) %>%
        as.data.frame()

    atac_sc1 <- import_from_matrix(mat, peaks, "hg38", id = atac_sc@id, rm_zero_peaks = FALSE)
    expect_true(sum(abs(mat - atac_sc1@mat)) <= 1e-9)
    expect_equal(peaks, atac_sc1@peaks %>% select(chrom, start, end) %>% as.data.frame())

    # mcATAC
    atac_mc1 <- import_from_matrix(atac_mc@mat, atac_mc@peaks, atac_mc@genome, class = "McPeaks", id = atac_mc@id, description = atac_mc@description, metadata = atac_mc@metadata)
    expect_true(mean(abs(atac_mc@mat - atac_mc1@mat)) <= 1e-9)
    expect_true(mean(abs(atac_mc@egc - atac_mc1@egc)) <= 1e-9)
    expect_true(mean(abs(atac_mc@fp - atac_mc1@fp)) <= 1e-9)
    expect_equal(atac_mc@peaks, atac_mc1@peaks, ignore_attr = TRUE)
    expect_equal(atac_mc@mc_size_eps_q, atac_mc1@mc_size_eps_q)
    expect_equal(atac_mc@genome, atac_mc1@genome)
    expect_equal(atac_mc@metadata, atac_mc1@metadata, ignore_attr = TRUE)
    expect_equal(atac_mc@tad_based, atac_mc1@tad_based)
    expect_equal(atac_mc@description, atac_mc1@description)
})

test_that("projection works", {
    expect_equal(class(atac_mc), "McPeaks", ignore_attr = TRUE)
    expect_equal(nrow(atac_mc@mat), 107687)
    expect_equal(ncol(atac_mc@mat), length(unique(cell_to_metacell_pbmc_example$metacell)))
    expect_setequal(colnames(atac_mc@mat), cell_to_metacell_pbmc_example$metacell)
    expect_equal(colnames(atac_mc@peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_mc@peaks), atac_mc@peaks$peak_name)
    # make sure we didn't change the peak names form the scATAC
    expect_equal(
        atac_mc@peaks %>%
            left_join(atac_sc@peaks %>% rename(peak_name_old = peak_name), by = c("chrom", "start", "end")) %>%
            filter(peak_name != peak_name_old) %>%
            nrow(),
        0
    )
    expect_equal(colSums(atac_mc@egc), rep(quantile(colSums(atac_mc@mat), 0.1), ncol(atac_mc@egc)), ignore_attr = TRUE)
    expect_equal(atac_mc@id, "pbmc")
    expect_equal(atac_mc@tad_based, TRUE)
    expect_equal(atac_mc@description, "PBMC from a healthy donor - granulocytes removed through cell sorting (10k)")
    expect_equal(atac_mc@cell_to_metacell, cell_to_metacell_pbmc_example)
})

test_that("projection from a metacell1 object works", {
    atac_mc1 <- project_atac_on_mc_from_metacell1(atac_sc, fs::path(raw_dir, "scdb"), "rna")
    expect_equal(atac_mc@mat, atac_mc1@mat)
    expect_equal(atac_mc@peaks, atac_mc1@peaks)
})

test_that("export works", {
    export_to_h5ad(atac_mc, fs::path(raw_dir, "atac_mc.h5ad"), compression = "gzip")
    atac_mc1 <- import_from_h5ad(fs::path(raw_dir, "atac_mc.h5ad"))
    expect_true(mean(abs(atac_mc@mat - atac_mc1@mat)) <= 1e-9)
    expect_true(mean(abs(atac_mc@egc - atac_mc1@egc)) <= 1e-9)
    expect_true(mean(abs(atac_mc@fp - atac_mc1@fp)) <= 1e-9)
    expect_equal(atac_mc@peaks, atac_mc1@peaks, ignore_attr = TRUE)
    expect_equal(atac_mc@mc_size_eps_q, atac_mc1@mc_size_eps_q)
    expect_equal(atac_mc@genome, atac_mc1@genome)
    expect_equal(atac_mc@metadata, atac_mc1@metadata, ignore_attr = TRUE)
    expect_equal(atac_mc@id, atac_mc1@id)
    expect_equal(atac_mc@tad_based, atac_mc1@tad_based)
    expect_equal(atac_mc@description, atac_mc1@description)
    expect_equal(atac_mc1@path, normalizePath(fs::path(raw_dir, "atac_mc.h5ad")))
    expect_equal(atac_mc1@cell_to_metacell, atac_mc@cell_to_metacell, ignore_attr = TRUE)
})

test_that("add_metadata works", {
    data(mcmd)
    atac_mc1 <- add_mc_metadata(atac_mc, mcmd)
    expect_true(all(colnames(atac_mc1@mat) == atac_mc1@metadata$metacell))
})
