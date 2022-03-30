expect_ggplot_ok <- function(gg, fn = tempfile()) {
    png(fn)
    print(gg)
    dev.off()
    expect_true(file.exists(fn))
    expect_gt(file.size(fn), 0)
}
