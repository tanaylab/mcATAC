test_that("import from 10x works", {
    expect_equal(class(atac_sc), "ScATAC", ignore_attr = TRUE)
    expect_equal(nrow(atac_sc@mat), 107829)
    expect_equal(ncol(atac_sc@mat), 11909)
    expect_equal(colnames(atac_sc@peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_sc@peaks), atac_sc@peaks$peak_name)
})

test_that("projection works", {
    expect_equal(class(atac_mc), "McATAC", ignore_attr = TRUE)
    expect_equal(nrow(atac_mc@mat), 107687)
    expect_equal(ncol(atac_mc@mat), length(unique(cell_to_metacell_pbmc_example$metacell)))
    expect_setequal(colnames(atac_mc@mat), cell_to_metacell_pbmc_example$metacell)
    expect_equal(colnames(atac_mc@peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_mc@peaks), atac_mc@peaks$peak_name)
})

test_that("projection from a metacell1 object works", {
    atac_mc1 <- project_atac_on_mc_from_metacell1(atac_sc, fs::path(raw_dir, "scdb"), "rna")
    expect_equal(atac_mc@mat, atac_mc1@mat)
    expect_equal(atac_mc@peaks, atac_mc1@peaks)
})

test_that("export works", {
    export_to_h5ad(atac_mc, fs::path(raw_dir, "atac_mc.h5ad"), compression = "gzip")
    atac_mc1 <- import_from_h5ad(fs::path(raw_dir, "atac_mc.h5ad"))
    expect_true(mean(atac_mc@mat - atac_mc1@mat) <= 1e-9)
    expect_equal(atac_mc@peaks, atac_mc1@peaks, ignore_attr = TRUE)
})

test_that("add_metadata works", {
    data(mcmd)
    atac_mc1 <- add_mc_metadata(atac_mc, mcmd)
    expect_true(all(colnames(atac_mc1@mat) == atac_mc1@metadata$metacell))
})
