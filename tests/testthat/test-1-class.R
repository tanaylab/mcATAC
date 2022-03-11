test_that("import from 10x works", {
    atac_sc <<- import_from_10x(raw_dir)
    expect_equal(class(atac_sc), "ScATAC")
    expect_true(all(has_name(atac_sc, c("mat", "peaks", "metadata"))))
    expect_equal(nrow(atac_sc$mat), 108377)
    expect_equal(ncol(atac_sc$mat), 11909)
    expect_equal(colnames(atac_sc$peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_sc$peaks), atac_sc$peaks$peak_name)
})

test_that("projection works", {
    atac_mc <- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example)
    expect_equal(class(atac_mc), "McATAC")
    expect_equal(nrow(atac_mc$mat), 108377)
    expect_equal(ncol(atac_mc$mat), length(unique(cell_to_metacell_pbmc_example$metacell)))
    expect_setequal(colnames(atac_mc$mat), cell_to_metacell_pbmc_example$metacell)
    expect_true(all(has_name(atac_mc, c("mat", "peaks", "metadata"))))
    expect_equal(colnames(atac_mc$peaks), c("chrom", "start", "end", "peak_name"))
    expect_equal(peak_names(atac_mc$peaks), atac_mc$peaks$peak_name)
})

test_that("projection from a metacell1 object works", {
    atac_mc <- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example)
    atac_mc1 <- project_atac_on_mc_from_metacell1(atac_sc, fs::path(raw_dir, "scdb"), "rna")
    expect_equal(atac_mc, atac_mc1)
})

test_that("export works", {
    atac_mc <- project_atac_on_mc(atac_sc, cell_to_metacell_pbmc_example)
    export_to_h5ad(atac_mc, fs::path(raw_dir, "atac_mc.h5ad"), compression = "gzip")
    atac_mc1 <- import_from_h5ad(fs::path(raw_dir, "atac_mc.h5ad"))
    expect_equal(atac_mc, atac_mc1, ignore_attr = TRUE)
})
