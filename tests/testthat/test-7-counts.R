
test_that("write_sparse_matrix_from_bam works", {
    gset_genome("hg38")
    reg <- giterator.intervals(intervals = gintervals.all() %>%
        filter(chrom == "chrY"), iterator = 5e7) %>%
        slice(1)
    bam_file <- file.path(raw_dir, "possorted_bam.bam")
    mm_file <- tempfile()
    write_sparse_matrix_from_bam(bam_file,
        out_file = mm_file, cell_names = atac_sc,
        region = reg
    )
    expect_true(file.exists(paste0(mm_file, ".gz")))
    mat <- tgutil::fread_mm(paste0(mm_file, ".gz"))
    expect_equal(dim(mat), c(50000000L, 11909))
    expect_equal(sum(mat), 99340)

    # Extract the reads directly from the BAM file and compare
    cell_names <- colnames(atac_sc@mat)
    region_str <- paste0(reg$chrom, ":", reg$start, "-", reg$end)
    cmd <- glue("{samtools_bin} view --keep-tag CB {bam_file} {region_str} | grep CB | {awk_cmd}", awk_cmd = "awk '{print $4,substr($12, 6)}'")
    df <- tgutil::fread(cmd = cmd, col.names = c("start", "cell_name"), header = FALSE) %>% as_tibble()
    df_count <- df %>%
        filter(cell_name %in% cell_names) %>%
        count(start, cell_name)

    mm_df <- Matrix::summary(mat) %>%
        as_tibble() %>%
        rlang::set_names(c("start", "cell_name", "n_mm")) %>%
        mutate(cell_name = cell_names[cell_name])

    expect_equal(
        mm_df %>%
            anti_join(df_count) %>%
            nrow(),
        0
    )

    expect_equal(
        df_count %>%
            anti_join(mm_df) %>%
            nrow(),
        0
    )

    expect_equal(
        mm_df %>%
            left_join(df_count) %>%
            filter(n_mm != n) %>%
            nrow(),
        0
    )
})
